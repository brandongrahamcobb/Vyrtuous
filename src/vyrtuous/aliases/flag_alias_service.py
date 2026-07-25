"""!/bin/python3
flag_service.py The purpose of this program is to extend AliasService to service flag infractions.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from dataclasses import dataclass

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot, TargetIsBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.flag import Flag
from vyrtuous.models.duration import DurationObject
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.tracking import data_builder, stream_service

MODEL = Flag


@dataclass(frozen=True)
class FlagMessageContext:
    author_snowflake: int
    channel_snowflake: int
    guild_snowflake: int
    member_snowflake: int
    message_snowflake: int
    message_channel_snowflake: int
    reason: str


async def log_flag(
    author_snowflake: int | None,
    channel_snowflake: int,
    display: bool,
    guild_snowflake: int,
    member_snowflake: int,
    message_snowflake: int | None,
    message_channel_snowflake: int | None,
    reason: str | None,
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
        identifier="flag",
        member_snowflake=member_snowflake,
        reason=reason or "No reason provided.",
        role_snowflake=role_snowflake or None,
        target=target or None,
    )
    if display:
        await stream_service.send_log(
            author_snowflake=author_snowflake or None,
            channel_snowflake=channel_snowflake,
            duration=duration,
            guild_snowflake=guild_snowflake,
            identifier="flag",
            is_channel_scope=is_channel_scope or None,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake or None,
            message_channel_snowflake=message_channel_snowflake or None,
            reason=reason or "No reason provided.",
            role_snowflake=role_snowflake or None,
            target=target or None,
        )


async def flag_by_message(
    ctx: FlagMessageContext, display: bool = True
) -> discord.Embed:
    await flag(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        reason=ctx.reason,
    )
    await log_flag(
        author_snowflake=ctx.author_snowflake,
        channel_snowflake=ctx.channel_snowflake,
        display=display,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        message_snowflake=ctx.message_snowflake,
        message_channel_snowflake=ctx.message_channel_snowflake,
        reason=ctx.reason,
    )
    embed = build_flag_embed(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        reason=ctx.reason,
    )
    return embed


async def flag(
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
    reason: str,
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    if guild.me.id == member_snowflake:
        raise TargetIsBot
    member = guild.get_member(member_snowflake)
    if member is None:
        simplified_member = bot.registry.get(MemberState).active.get(member_snowflake)
        if not simplified_member:
            raise commands.MemberNotFound(str(member_snowflake))
    flag = MODEL(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
        reason=reason,
    )
    await database_factory.create(flag)
    original_set = bot.registry.get(MemberState).flagged
    original_set[guild_snowflake].add(member_snowflake)


def build_flag_embed(
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
    reason: str,
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
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
        title=f"{emojis.get_random_emoji()} {display_name} has been flagged",
        description=(
            f"**User:** {member_str}\n"
            f"**Channel:** {channel.mention}\n"
            f"**Reason:** {reason}"
        ),
        color=discord.Color.blue(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed
