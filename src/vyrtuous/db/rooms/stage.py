"""stage.py The purpose of this program is to provide the Stage utility class.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from datetime import datetime
from typing import Optional
import time

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.fields.duration import DurationObject
from vyrtuous.db.actions.voice_mute import VoiceMute
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    clean_guild_dictionary,
    flush_page,
)
from vyrtuous.utils.check import (
    check,
)
from vyrtuous.utils.logger import logger
from vyrtuous.db.roles.administrator import is_administrator
from vyrtuous.db.roles.coordinator import is_coordinator
from vyrtuous.db.roles.developer import is_developer
from vyrtuous.db.roles.guild_owner import is_guild_owner
from vyrtuous.db.roles.moderator import (
    is_moderator,
)
from vyrtuous.db.roles.sysadmin import is_sysadmin
from vyrtuous.utils.check import has_equal_or_lower_role


class Stage(DatabaseFactory):

    ACT = "stage"
    CATEGORY = "stage"
    PLURAL = "Stages"
    SCOPES = ["channels"]
    SINGULAR = "Stage"
    UNDO = "stage"

    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "expires_in",
        "guild_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "expired",
        "updated_at",
    ]

    TABLE_NAME = "active_stages"

    def __init__(
        self,
        channel_snowflake: int,
        expires_in: Optional[datetime],
        guild_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.expires_in = expires_in
        self.guild_snowflake = guild_snowflake
        self.updated_at = updated_at

    async def send_stage_ask_to_speak_message(
        self, join_log: dict[int, discord.Member], member: discord.Member
    ):
        bot = DiscordBot.get_instance()
        now = time.time()
        join_log[member.id] = [t for t in join_log[member.id] if now - t < 300]
        if len(join_log[member.id]) < 1:
            join_log[member.id].append(now)
            embed = discord.Embed(
                title=f"{get_random_emoji()} {self.channel_snowflake} — Stage Mode",
                description=f"Ends <t:{int(self.expires_in.timestamp())}:R>",
                color=discord.Color.green(),
            )
            embed.add_field(name="\u200b", value="**Ask to speak!**", inline=False)
            await bot.get_channel(self.channel_snowflake).send(embed=embed)

    @classmethod
    async def build_dictionary(cls, where_kwargs):
        dictionary = {}
        stages = await Stage.select(**where_kwargs)
        for stage in stages:
            dictionary.setdefault(stage.guild_snowflake, {"channels": {}})
            dictionary[stage.guild_snowflake]["channels"].setdefault(
                stage.channel_snowflake, {}
            )
            dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ].setdefault("stages", {})
            dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ]["stages"].update(
                {"expires_in": DurationObject.from_expires_in(stage.expires_in)}
            )
        return dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} Stages"
        where_kwargs = object_dict.get("columns", None)

        guild_dictionary = await Stage.build_dictionary(where_kwargs=where_kwargs)

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, stage_dictionary in guild_data.get(
                "channels"
            ).items():
                channel = guild.get_channel(channel_snowflake)
                lines.append(
                    f"**Expires in:** {stage_dictionary.get("expires_in", None)}"
                )
                field_count += 1
                if field_count == chunk_size:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
                if lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
            pages.append(embed)

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
        return pages

    @classmethod
    async def toggle_stage(cls, channel_dict, duration, snowflake_kwargs):
        bot = DiscordBot.get_instance()
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        member_snowflake = snowflake_kwargs.get("member_snowflake", None)
        guild = bot.get_guild(guild_snowflake)
        failed, pages, skipped, succeeded = [], [], [], []
        kwargs = channel_dict.get("columns", None)

        stage = await Stage.select(**kwargs, singular=True)
        if stage:
            title = f"{get_random_emoji()} Stage Ended in {channel_dict.get("mention", None)}"
            await Stage.delete(**kwargs)
            for member in channel_dict.get("object", None).members:
                await VoiceMute.delete(
                    **kwargs,
                    member_snowflake=member.id,
                    target="room",
                )
                voice_mute = await VoiceMute.select(
                    **kwargs, member_snowflake=member.id, target="user", singular=True
                )
                if not voice_mute and member.voice and member.voice.mute:
                    try:
                        await member.edit(
                            mute=False,
                            reason="Stage ended — no user-specific mute found",
                        )
                        succeeded.append(member)
                    except discord.Forbidden as e:
                        logger.warning(
                            f"Unable to undo voice-mute "
                            f"for member {member.display_name} ({member.id}) in "
                            f"channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) in "
                            f"guild {guild.name} ({guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict.get("mention", None)}",
                f"**Unmuted:** {len(succeeded)} users",
            ]
            if failed:
                description_lines.append(f"**Failed:** {len(failed)}")
            embed = discord.Embed(
                description="\n".join(description_lines),
                title=title,
                color=discord.Color.blurple(),
            )
            pages.append(embed)
        else:
            stage = Stage(
                **kwargs,
                expires_in=duration.expires_in,
            )
            await stage.create()
            for member in channel_dict.get("object", None).members:
                if (
                    await check(
                        snowflake_kwargs=snowflake_kwargs,
                        lowest_role="Coordinator",
                    )
                    or member.id == member_snowflake
                ):
                    skipped.append(member)
                    continue
                voice_mute = await VoiceMute(
                    **kwargs,
                    expires_in=duration.expires_in,
                    member_snowflake=member.id,
                    target="room",
                    reason="Stage mute",
                )
                await voice_mute.create()
                try:
                    if member.voice and member.voice.channel.id == channel_dict.get(
                        "id", None
                    ):
                        await member.edit(mute=True)
                    succeeded.append(member)
                except Exception as e:
                    logger.warning(
                        f"Unable to voice-mute "
                        f"member {member.display_name} ({member.id}) "
                        f"in channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) "
                        f"in guild {guild.name} ({guild.id}). "
                        f"{str(e).capitalize()}"
                    )
                    failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict.get("mention", None)}",
                f"**Expires:** {duration}",
                f"**Muted:** {len(succeeded)} users",
                f"**Skipped:** {len(skipped)}",
            ]
            if failed:
                description_lines.append(f"**Failed:** {len(failed)}")
            embed = discord.Embed(
                description="\n".join(description_lines),
                title=f"{get_random_emoji()} Stage Created in {channel_dict.get('name', None)}",
                color=discord.Color.blurple(),
            )
            pages.append(embed)
        return pages

    @classmethod
    async def survey(cls, channel_dict, guild_snowflake):

        chunk_size, pages = 7, []
        (
            sysadmins,
            developers,
            guild_owners,
            administrators,
            coordinators,
            moderators,
        ) = ([], [], [], [], [], [])

        for member in channel_dict.get("object", None).members:
            try:
                if await is_sysadmin(member.id):
                    sysadmins.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_developer(member.id):
                    developers.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_guild_owner(guild_snowflake, member.id):
                    guild_owners.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_administrator(guild_snowflake, member.id):
                    administrators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_coordinator(
                    channel_dict.get("id", None), guild_snowflake, member.id
                ):
                    coordinators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_moderator(
                    channel_dict.get("id", None), guild_snowflake, member.id
                ):
                    moderators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
        sysadmins_chunks = [
            sysadmins[i : i + chunk_size] for i in range(0, len(sysadmins), chunk_size)
        ]
        guild_owners_chunks = [
            guild_owners[i : i + chunk_size]
            for i in range(0, len(guild_owners), chunk_size)
        ]
        developers_chunks = [
            developers[i : i + chunk_size]
            for i in range(0, len(developers), chunk_size)
        ]
        administrators_chunks = [
            administrators[i : i + chunk_size]
            for i in range(0, len(administrators), chunk_size)
        ]
        coordinators_chunks = [
            coordinators[i : i + chunk_size]
            for i in range(0, len(coordinators), chunk_size)
        ]
        moderators_chunks = [
            moderators[i : i + chunk_size]
            for i in range(0, len(moderators), chunk_size)
        ]
        roles_chunks = [
            ("Sysadmins", sysadmins, sysadmins_chunks),
            ("Developers", developers, developers_chunks),
            ("Guild Owners", guild_owners, guild_owners_chunks),
            ("Administrators", administrators, administrators_chunks),
            ("Coordinators", coordinators, coordinators_chunks),
            ("Moderators", moderators, moderators_chunks),
        ]
        max_pages = max(len(c[2]) for c in roles_chunks)
        for page in range(max_pages):
            embed = discord.Embed(
                title=f"{get_random_emoji()} Survey results for {channel_dict.get('name', None)}",
                description=f"Total surveyed: {len(channel_dict.get('object', None).members)}",
                color=discord.Color.blurple(),
            )
            for role_name, role_list, chunks in roles_chunks:
                chunk = chunks[page] if page < len(chunks) else []
                embed.add_field(
                    name=f"{role_name} ({len(chunk)}/{len(role_list)})",
                    value=", ".join(u.mention for u in chunk) if chunk else "*None*",
                    inline=False,
                )
            pages.append(embed)
        return pages

    @classmethod
    async def toggle_stage_mute(cls, channel_dict, member_dict, snowflake_kwargs):
        await has_equal_or_lower_role(
            snowflake_kwargs=snowflake_kwargs,
            member_snowflake=member_dict.get("id", None),
        )
        where_kwargs = channel_dict.get("columns", None)
        stage = await Stage.select(**where_kwargs, singular=True)
        if stage:
            await member_dict.get("object", None).edit(
                mute=not member_dict.get("object", None).voice.mute
            )
            return f"Successfully toggled the mute for {member_dict.get("mention", None)} in {channel_dict.get("mention", None)}."
