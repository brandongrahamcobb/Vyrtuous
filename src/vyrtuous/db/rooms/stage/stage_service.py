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

import time

import discord
from discord.ext import commands

from vyrtuous.base.service import Service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.fields.duration import DurationObject
from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.infractions.vmute.voice_mute import VoiceMute
from vyrtuous.db.roles.admin.administrator_service import is_administrator
from vyrtuous.db.roles.coord.coordinator_service import is_coordinator
from vyrtuous.db.roles.dev.developer_service import is_developer
from vyrtuous.db.roles.mod.moderator_service import is_moderator
from vyrtuous.db.roles.owner.guild_owner_service import is_guild_owner
from vyrtuous.db.roles.sysadmin.sysadmin_service import is_sysadmin
from vyrtuous.db.rooms.stage.stage import Stage
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_channels,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger


class StageService(Service):

    lines, pages = [], []

    @classmethod
    async def send_stage_ask_to_speak_message(
        cls, join_log: dict[int, discord.Member], member: discord.Member, stage: Stage
    ):
        bot = DiscordBot.get_instance()
        now = time.time()
        join_log[member.id] = [t for t in join_log[member.id] if now - t < 300]
        if len(join_log[member.id]) < 1:
            join_log[member.id].append(now)
            embed = discord.Embed(
                title=f"{get_random_emoji()} {stage.channel_snowflake} — Stage Mode",
                description=f"Ends <t:{int(stage.expires_in.timestamp())}:R>",
                color=discord.Color.green(),
            )
            embed.add_field(name="\u200b", value="**Ask to speak!**", inline=False)
            await bot.get_channel(stage.channel_snowflake).send(embed=embed)

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        stages = await Stage.select(singular=False, **where_kwargs)
        for stage in stages:
            dictionary.setdefault(stage.guild_snowflake, {"channels": {}})
            dictionary[stage.guild_snowflake]["channels"].setdefault(
                stage.channel_snowflake, {}
            )
            dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ].setdefault("stages", {})
            dictionary[stage.guild_snowflake]["channels"][stage.channel_snowflake][
                "stages"
            ].update({"expires_in": DurationObject.from_expires_in(stage.expires_in)})
        skipped_channels = generate_skipped_channels(dictionary)
        skipped_guilds = generate_skipped_guilds(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_guilds:
                StageService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_channels:
                StageService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Stages"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await StageService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            stage_n = 0
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, stage_dictionary in guild_data.get(
                "channels"
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                StageService.lines.append(
                    f"**Expires in:** {stage_dictionary.get("stages", {}).get("expires_in", None)}"
                )
                stage_n += 1
                field_count += 1
                if field_count == CHUNK_SIZE:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(StageService.lines),
                        inline=False,
                    )
                    embed = flush_page(embed, StageService.pages, title, guild.name)
                    StageService.lines = []
                    field_count = 0
                if StageService.lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(StageService.lines),
                        inline=False,
                    )
            StageService.pages.append(embed)
            StageService.pages[0].description = f'{guild.name} **({stage_n})**'
        return StageService.pages

    @classmethod
    async def toggle_stage(cls, channel_dict, default_kwargs, duration):
        bot = DiscordBot.get_instance()
        updated_kwargs = default_kwargs.copy()
        guild = bot.get_guild(updated_kwargs.get("guild_snowflake", None))
        failed, pages, skipped, succeeded = [], [], [], []
        updated_kwargs.update(channel_dict.get("columns", None))
        stage_kwargs = updated_kwargs.copy()
        del stage_kwargs["member_snowflake"]
        stage = await Stage.select(**stage_kwargs, singular=True)
        if stage:
            title = f"{get_random_emoji()} Stage Ended in {channel_dict.get("mention", None)}"
            await Stage.delete(**updated_kwargs)
            for member in channel_dict.get("object", None).members:
                await VoiceMute.delete(
                    **updated_kwargs,
                    member_snowflake=member.id,
                    target="room",
                )
                voice_mute = await VoiceMute.select(
                    **updated_kwargs,
                    member_snowflake=member.id,
                    target="user",
                    singular=True,
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
                **stage_kwargs,
                expires_in=duration.expires_in,
            )
            await stage.create()
            for member in channel_dict.get("object", None).members:
                if await PermissionService.check(
                    updated_kwargs=updated_kwargs,
                    lowest_role="Coordinator",
                ) or member.id == updated_kwargs.get("member_snowflake", None):
                    skipped.append(member)
                    continue
                voice_mute = await VoiceMute(
                    **updated_kwargs,
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
    async def toggle_stage_mute(cls, channel_dict, default_kwargs, member_dict):
        await PermissionService.has_equal_or_lower_role(
            updated_kwargs=default_kwargs,
            member_snowflake=member_dict.get("id", None),
        )
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(channel_dict.get("columns", None))
        stage = await Stage.select(singular=True, **updated_kwargs)
        if stage:
            await member_dict.get("object", None).edit(
                mute=not member_dict.get("object", None).voice.mute
            )
            return f"Successfully toggled the mute for {member_dict.get("mention", None)} in {channel_dict.get("mention", None)}."
