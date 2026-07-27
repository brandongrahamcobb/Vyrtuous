"""!/bin/python3
flag_service.py The purpose of this program is to extend AliasService to service flag infractions.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from dataclasses import dataclass

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.flag import Flag
from vyrtuous.models.duration import DurationObject
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.tracking import data_builder, stream_service

MODEL = Flag


@dataclass(frozen=True)
class UnflagMessageContext:
    author_snowflake: int
    channel_snowflake: int
    guild_snowflake: int
    member_snowflake: int
    message_snowflake: int
    message_channel_snowflake: int
    reason: str = "No reason provided."


async def unflag_by_message(
    ctx: UnflagMessageContext, display: bool = True
) -> discord.Embed:
    await unflag(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
    )
    await log_unflag(
        author_snowflake=ctx.author_snowflake,
        channel_snowflake=ctx.channel_snowflake,
        display=display,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        message_snowflake=ctx.message_snowflake,
        message_channel_snowflake=ctx.message_channel_snowflake,
        reason=ctx.reason,
    )
    embed = build_unflag_embed(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
    )
    return embed


async def log_unflag(
    author_snowflake: int | None,
    channel_snowflake: int,
    display: bool,
    guild_snowflake: int,
    member_snowflake: int,
    message_snowflake: int | None,
    message_channel_snowflake: int | None,
    reason: str,
) -> None:
    duration = DurationObject(number=0, prefix="", sign=1, unit="")
    is_channel_scope = None
    role_snowflake = None
    target = None
    await data_builder.save_data(
        author_snowflake=author_snowflake or None,
        channel_snowflake=channel_snowflake,
        duration=duration,
        guild_snowflake=guild_snowflake,
        identifier="unflag",
        member_snowflake=member_snowflake,
        reason=reason,
        role_snowflake=role_snowflake or None,
        target=target or None,
    )
    if display:
        await stream_service.send_log(
            author_snowflake=author_snowflake or None,
            channel_snowflake=channel_snowflake,
            duration=duration,
            guild_snowflake=guild_snowflake,
            identifier="unflag",
            is_channel_scope=is_channel_scope or None,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake or None,
            message_channel_snowflake=message_channel_snowflake or None,
            reason=reason,
            role_snowflake=role_snowflake or None,
            target=target or None,
        )


async def unflag(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    member = guild.get_member(member_snowflake)
    if member is None:
        simplified_member = bot.registry.get(MemberState).active.get(member_snowflake)
        if not simplified_member:
            raise commands.MemberNotFound(str(member_snowflake))
    await database_factory.delete(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
    )
    original_set = bot.registry.get(MemberState).flagged
    original_set[guild_snowflake].remove(member_snowflake)


def build_unflag_embed(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(str(channel_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        display_name = member.display_name
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(member_snowflake)
        if not simplified_member:
            raise commands.MemberNotFound(str(member_snowflake))
        display_name = simplified_member[0]
        member_str = display_name
    embed = discord.Embed(
        title=f"{emojis.get_random_emoji()} {display_name} has been unflagged",
        description=(f"**User:** {member_str}\n**Channel:** {channel.mention}"),
        color=discord.Color.yellow(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed
