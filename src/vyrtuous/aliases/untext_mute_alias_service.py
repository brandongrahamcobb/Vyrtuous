"""!/bin/python3
text_mute_service.py The purpose of this program is to extend AliasService to service the text mute infraction.

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

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.text_mute import TextMute
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.moderation import cap_service
from vyrtuous.utils.tracking import data_builder, stream_service
from vyrtuous.utils.users import moderator_service

MODEL = TextMute


def build_untext_mute_embed(
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
) -> discord.Embed:
    bot = DiscordBot.get_instance()
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
        title=f"{emojis.get_random_emoji()} {display_name} has been untext-muted",
        description=(f"**User:** {member_str}\n**Channel:** {channel.mention}"),
        color=discord.Color.yellow(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed


@dataclass(frozen=True)
class UntextMuteMessageContext:
    author_snowflake: int
    channel_snowflake: int
    guild_snowflake: int
    member_snowflake: int
    message_snowflake: int
    message_channel_snowflake: int


async def log_untext_mute(
    author_snowflake: int | None,
    channel_snowflake: int,
    display: bool,
    guild_snowflake: int,
    member_snowflake: int,
    message_snowflake: int | None,
    message_channel_snowflake: int | None,
) -> None:
    duration_value = None
    is_channel_scope = None
    reason = None
    role_snowflake = None
    target = None
    await data_builder.save_data(
        author_snowflake=author_snowflake or None,
        channel_snowflake=channel_snowflake,
        duration_value=duration_value or None,
        guild_snowflake=guild_snowflake,
        identifier="untmute",
        member_snowflake=member_snowflake,
        reason=reason or "No reason provided.",
        role_snowflake=role_snowflake or None,
        target=target or None,
    )
    if display:
        await stream_service.send_log(
            author_snowflake=author_snowflake or None,
            channel_snowflake=channel_snowflake,
            identifier="untmute",
            duration_value=duration_value or None,
            guild_snowflake=guild_snowflake,
            is_channel_scope=is_channel_scope or None,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake or None,
            message_channel_snowflake=message_channel_snowflake or None,
            reason=reason or "No reason provided.",
            role_snowflake=role_snowflake or None,
            target=target or None,
        )


async def set_untext_mute_overwrite(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> None:
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    channel = bot.get_channel(channel_snowflake)
    if channel is None or not isinstance(
        channel, (discord.VoiceChannel, discord.StageChannel)
    ):
        raise commands.ChannelNotFound(str(channel_snowflake))
    member = guild.get_member(member_snowflake)
    if member is None:
        simplified_member = bot.registry.get(MemberState).active.get(member_snowflake)
        if not simplified_member:
            raise commands.MemberNotFound(str(member_snowflake))
    else:
        try:
            await channel.set_permissions(
                target=member,
                send_messages=None,
                add_reactions=None,
                reason="Undoing text-mute",
            )
        except discord.Forbidden:
            raise


async def untext_mute(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> None:
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(str(channel_snowflake))
    member = guild.get_member(member_snowflake)
    if member is None:
        simplified_member = bot.registry.get(MemberState).active.get(member_snowflake)
        if not simplified_member:
            raise commands.MemberNotFound(str(member_snowflake))
    else:
        await set_untext_mute_overwrite(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
    await database_factory.delete(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
    )


async def untext_mute_by_message(
    ctx: UntextMuteMessageContext, display: bool = True
) -> discord.Embed:
    database_factory = DatabaseFactory(MODEL)
    duration_builder = DurationBuilder()
    text_mute = await database_factory.select(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        singular=True,
    )
    duration_value = duration_builder.from_timestamp(text_mute.expires_in).build(
        as_str=True
    )
    exceeds_cap = await cap_service.exceeds_cap(
        category="tmute",
        channel_snowflake=ctx.channel_snowflake,
        duration_value=str(duration_value),
        guild_snowflake=ctx.guild_snowflake,
    )
    if exceeds_cap:
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel_snowflake,
            guild_snowflake=ctx.guild_snowflake,
            member_snowflake=ctx.author_snowflake,
            lowest_role="Coordinator",
        )
    await untext_mute(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
    )
    await log_untext_mute(
        author_snowflake=ctx.author_snowflake,
        channel_snowflake=ctx.channel_snowflake,
        display=display,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        message_snowflake=ctx.message_snowflake,
        message_channel_snowflake=ctx.message_channel_snowflake,
    )
    embed = build_untext_mute_embed(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
    )
    return embed
