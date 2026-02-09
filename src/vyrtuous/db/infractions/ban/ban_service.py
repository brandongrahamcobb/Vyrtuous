"""!/bin/python3
ban_service.py The purpose of this program is to extend AliasService to service ban infractions.

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

from datetime import datetime, timezone

import discord

from vyrtuous.base.record_service import RecordService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.fields.duration import DurationObject
from vyrtuous.db.infractions.ban.ban import Ban
from vyrtuous.db.mgmt.stream.stream_service import StreamService
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger


class BanService(RecordService):
    lines, pages = [], []
    model = Ban

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        bans = await Ban.select(singular=False, **where_kwargs)
        for ban in bans:
            dictionary.setdefault(ban.guild_snowflake, {"members": {}})
            dictionary[ban.guild_snowflake]["members"].setdefault(
                ban.member_snowflake, {"bans": {}}
            )
            dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                "bans"
            ].setdefault(ban.channel_snowflake, {})
            dictionary[ban.guild_snowflake]["members"][ban.member_snowflake]["bans"][
                ban.channel_snowflake
            ] = {
                "reason": ban.reason,
                "expires_in": DurationObject.from_expires_in(ban.expires_in),
            }
        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_members = generate_skipped_members(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        if is_at_home:
            if skipped_guilds:
                BanService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                BanService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_members,
                        title="Skipped Members in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        thumbnail = False
        where_kwargs = object_dict.get("columns", None)
        title = f"{get_random_emoji()} Bans {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        dictionary = await BanService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        ban_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, ban_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    BanService.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in ban_dictionary.get(
                    "bans"
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(
                        object_dict.get("object"), discord.abc.GuildChannel
                    ):
                        BanService.lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(object_dict.get("object"), discord.Member):
                        BanService.lines.append(
                            f"**Expires in:** {channel_dictionary['expires_in']}"
                        )
                        BanService.lines.append(
                            f"**Reason:** {channel_dictionary['reason']}"
                        )
                    ban_n += 1
                    field_count += 1
                    if field_count >= CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(BanService.lines),
                            inline=False,
                        )
                        embed = flush_page(embed, BanService.pages, title, guild.name)
                        BanService.lines = []
                        field_count = 0
            if BanService.lines:
                embed.add_field(
                    name="Information", value="\n".join(BanService.lines), inline=False
                )
            BanService.pages.append(embed)
        if BanService.pages:
            BanService.pages[0].description = f"**({ban_n})**"
        return BanService.pages

    @classmethod
    async def ban_overwrite(cls, channel, member):
        if channel:
            kwargs = {
                "channel_snowflake": channel.id,
                "guild_snowflake": channel.guild.id,
                "member_snowflake": member.id,
            }
            ban = await Ban.select(**kwargs, singular=True)
            if ban:
                targets = []
                for target, overwrite in channel.overwrites.items():
                    if any(value is not None for value in overwrite._values.values()):
                        if isinstance(target, discord.Member):
                            targets.append(target)
                if member not in targets:
                    try:
                        await channel.set_permissions(
                            member, view_channel=False, reason="Reinstating active ban."
                        )
                    except discord.Forbidden as e:
                        logger.warning(e)
                    if (
                        member.voice
                        and member.voice.channel
                        and member.voice.channel.id == channel.id
                    ):
                        try:
                            await member.move_to(None, reason="Reinstating active ban.")
                            set_kwargs = {
                                "last_kicked": datetime.now(timezone.utc),
                                "reset": False,
                            }
                            await Ban.update(set_kwargs=set_kwargs, where_kwargs=kwargs)
                        except discord.Forbidden as e:
                            logger.warning(e)

    @classmethod
    async def enforce(cls, ctx, source, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        ban = Ban(
            channel_snowflake=ctx.target_channel_snowflake,
            expires_in=ctx.expires_in,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
            reason=ctx.reason,
        )
        await ban.create()
        is_channel_scope = False
        channel = guild.get_channel(ctx.target_channel_snowflake)
        if channel:
            try:
                await channel.set_permissions(
                    member,
                    view_channel=False,
                    reason=ctx.reason,
                )
                if (
                    member.voice
                    and member.voice.channel
                    and member.voice.channel.id == channel.id
                ):
                    is_channel_scope = True
                    await member.move_to(None, reason=ctx.reason)
                    where_kwargs = {
                        "channel_snowflake": ctx.target_channel_snowflake,
                        "guild_snowflake": ctx.source_guild_snowflake,
                        "member_snowflake": ctx.target_member_snowflake,
                    }
                    set_kwargs = {"last_kicked": datetime.now(timezone.utc)}
                    await Ban.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
            except discord.Forbidden as e:
                logger.error(str(e).capitalize())

                return await state.end(error=str(e).capitalize())
        await StreamService.send_entry(
            channel_snowflake=ctx.channel_snowflake,
            duration=DurationObject.from_expires_in(ctx.expires_in),
            identifier="ban",
            is_channel_scope=is_channel_scope,
            member=member,
            source=source,
            reason=ctx.reason,
        )
        embed = await BanService.act_embed(ctx=ctx)
        return await state.end(success=embed)

    @classmethod
    async def undo(cls, ctx, source, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        await Ban.delete(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
        )
        channel = guild.get_channel(ctx.target_channel_snowflake)
        if channel:
            try:
                await channel.set_permissions(member, view_channel=None)
            except discord.Forbidden as e:
                logger.error(str(e).capitalize())
                return await state.end(error=str(e).capitalize())
        await StreamService.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            identifier="unban",
            is_modification=True,
            member=member,
            source=source,
        )
        embed = await BanService.undo_embed(ctx=ctx)
        return await state.end(success=embed)

    @classmethod
    async def act_embed(cls, ctx):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(ctx.target_channel_snowflake)
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        embed = discord.Embed(
            title=f"{get_random_emoji()} {member.display_name} has been banned",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Expires:** {DurationObject.from_expires_in(ctx.expires_in)}\n"
                f"**Reason:** {ctx.reason}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, ctx):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(ctx.target_channel_snowflake)
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        embed = discord.Embed(
            title=f"{get_random_emoji()} {member.display_name} has been unbanned",
            description=(f"**User:** {member.mention}\n**Channel:** {channel.mention}"),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed
