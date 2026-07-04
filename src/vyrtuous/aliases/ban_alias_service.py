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

from dataclasses import dataclass
from datetime import datetime, timezone

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.ban import Ban
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.moderation import cap_service
from vyrtuous.utils.tracking import data_builder, stream_service
from vyrtuous.utils.users import moderator_service

MODEL = Ban


@dataclass(frozen=True)
class BanMessageContext:
    author_snowflake: int
    channel_snowflake: int
    duration_value: str
    guild_snowflake: int
    member_snowflake: int
    message_snowflake: int
    message_channel_snowflake: int
    reason: str


async def log_ban(
    author_snowflake: int | None,
    channel_snowflake: int,
    display: bool,
    duration_value: str | None,
    guild_snowflake: int,
    is_channel_scope: bool,
    member_snowflake: int,
    message_snowflake: int | None,
    message_channel_snowflake: int | None,
    reason: str,
):
    role_snowflake = None
    target = None
    await data_builder.save_data(
        author_snowflake=author_snowflake or None,
        channel_snowflake=channel_snowflake,
        duration_value=duration_value or None,
        guild_snowflake=guild_snowflake,
        identifier="ban",
        member_snowflake=member_snowflake,
        reason=reason or "No reason provided.",
        role_snowflake=role_snowflake or None,
        target=target or None,
    )
    if display:
        await stream_service.send_log(
            author_snowflake=author_snowflake or None,
            channel_snowflake=channel_snowflake,
            identifier="ban",
            duration_value=duration_value or None,
            guild_snowflake=guild_snowflake,
            is_channel_scope=is_channel_scope,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake or None,
            message_channel_snowflake=message_channel_snowflake or None,
            reason=reason or "No reason provided.",
            role_snowflake=role_snowflake or None,
            target=target or None,
        )


async def ban_by_message(
    ctx: BanMessageContext,
    display: bool = True,
) -> discord.Embed:
    if await cap_service.exceeds_cap(
        category="ban",
        channel_snowflake=ctx.channel_snowflake,
        duration_value=ctx.duration_value,
        guild_snowflake=ctx.guild_snowflake,
    ):
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel_snowflake,
            guild_snowflake=ctx.guild_snowflake,
            member_snowflake=ctx.author_snowflake,
            lowest_role="Coordinator",
        )
    is_channel_scope = await ban(
        channel_snowflake=ctx.channel_snowflake,
        duration_value=ctx.duration_value,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        reason=ctx.reason,
    )
    await log_ban(
        author_snowflake=ctx.author_snowflake,
        channel_snowflake=ctx.channel_snowflake,
        display=display,
        duration_value=ctx.duration_value,
        guild_snowflake=ctx.guild_snowflake,
        is_channel_scope=is_channel_scope,
        member_snowflake=ctx.member_snowflake,
        message_snowflake=ctx.message_snowflake,
        message_channel_snowflake=ctx.message_channel_snowflake or None,
        reason=ctx.reason,
    )
    embed = await build_ban_embed(
        channel_snowflake=ctx.channel_snowflake,
        duration_value=ctx.duration_value,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        reason=ctx.reason,
    )
    return embed


async def set_ban_overwrite(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int, reason: str
):
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(str(channel_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        try:
            await channel.set_permissions(
                member,
                view_channel=False,
                reason=reason,
            )
        except discord.Forbidden:
            raise


async def ban(
    channel_snowflake: int,
    duration_value: str,
    guild_snowflake: int,
    member_snowflake: int,
    reason: str,
):
    is_channel_scope = False
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(str(channel_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        await set_ban_overwrite(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            reason=reason,
        )
        if (
            member.voice
            and member.voice.channel
            and member.voice.channel.id == channel_snowflake
        ):
            is_channel_scope = True
            try:
                await member.move_to(None, reason=reason)
            except discord.Forbidden:
                raise
            where_kwargs = {
                "channel_snowflake": channel_snowflake,
                "guild_snowflake": guild_snowflake,
                "member_snowflake": member_snowflake,
            }
            set_kwargs = {"last_kicked": datetime.now(timezone.utc)}
            await database_factory.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
    else:
        simplified_member = bot.registry.get(MemberState).active.get(member_snowflake)
        if not simplified_member:
            raise commands.MemberNotFound(str(member_snowflake))
    duration_builder = DurationBuilder()
    expires_in = duration_builder.parse(duration_value).to_expires_in()
    ban = MODEL(
        channel_snowflake=channel_snowflake,
        expires_in=expires_in,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
        reason=reason,
    )
    await database_factory.create(ban)
    return is_channel_scope


async def build_ban_embed(
    channel_snowflake: int,
    duration_value: str,
    guild_snowflake: int,
    member_snowflake: int,
    reason: str,
):
    bot = DiscordBot.get_instance()
    duration_builder = DurationBuilder()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    channel = bot.get_channel(channel_snowflake)
    if channel is None or not isinstance(
        channel, (discord.VoiceChannel, discord.StageChannel)
    ):
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
        title=f"{emojis.get_random_emoji()} {display_name} has been banned",
        description=(
            f"**User:** {member_str}\n"
            f"**Channel:** {channel.mention}\n"
            f"**Expires:** {duration_builder.parse(duration_value).to_unix_ts()}\n"
            f"**Reason:** {reason}"
        ),
        color=discord.Color.blue(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed
