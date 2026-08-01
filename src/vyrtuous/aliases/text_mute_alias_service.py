"""!/bin/python3
text_mute_service.py The purpose of this program is to extend AliasService to service the text mute infraction.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from vyrtuous.bot.discord_bot import DiscordBot, TargetIsBot
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.text_mute import TextMute
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import ChannelNotFound, GuildNotFound, MemberNotFound
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.moderation import cap_service
from vyrtuous.utils.tracking import stream_service

MODEL = TextMute


@dataclass(frozen=True)
class TextMuteMessageContext:
    author_snowflake: int
    channel_snowflake: int
    duration: DurationObject
    guild_snowflake: int
    member_snowflake: int
    message_snowflake: int
    message_channel_snowflake: int
    reason: str


async def set_text_mute_overwrite(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int, reason: str
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise ChannelNotFound(str(channel_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        try:
            await channel.set_permissions(
                target=member,
                send_messages=False,
                add_reactions=False,
                reason=reason,
            )
        except discord.Forbidden:
            raise


async def text_mute_by_message(
    ctx: TextMuteMessageContext, display: bool = True
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    permission_state: PermissionState = bot.registry.get(PermissionState)
    if await cap_service.exceeds_cap(
        category="tmute",
        channel_snowflake=ctx.channel_snowflake,
        duration=ctx.duration,
        guild_snowflake=ctx.guild_snowflake,
    ):
        await permission_service.has_permissions(
            permission_state=permission_state,
            channel_snowflake=ctx.channel_snowflake,
            guild_snowflake=ctx.guild_snowflake,
            member_snowflake=ctx.author_snowflake,
            requested=["command.moderation.uncapped"],
        )
    is_channel_scope = await enable(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        reason=ctx.reason,
    )
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    duration_builder = DurationBuilder()
    expires_in = duration_builder.load(ctx.duration).to_expires_in()
    text_mute = MODEL(
        channel_snowflake=ctx.channel_snowflake,
        expires_in=expires_in,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        reason=ctx.reason,
    )
    await database_factory.create(text_mute)
    await stream_service.log(
        author_snowflake=ctx.author_snowflake,
        channel_snowflake=ctx.channel_snowflake,
        display=display,
        duration=ctx.duration,
        guild_snowflake=ctx.guild_snowflake,
        identifier="tmute",
        is_channel_scope=is_channel_scope,
        member_snowflake=ctx.member_snowflake,
        message_snowflake=ctx.message_snowflake,
        message_channel_snowflake=ctx.message_channel_snowflake,
        reason=ctx.reason,
        target=None,
        role_snowflake=None,
    )
    embed = build_text_mute_embed(
        channel_snowflake=ctx.channel_snowflake,
        duration=ctx.duration,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        reason=ctx.reason,
    )
    return embed


async def enable(
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
    reason: str,
) -> bool:
    is_channel_scope = False
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    if guild.me.id == member_snowflake:
        raise TargetIsBot
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise ChannelNotFound(str(channel_snowflake))
    member = guild.get_member(member_snowflake)
    if member is None:
        simplified_member = bot.registry.get(MemberState).active.get(member_snowflake)
        if not simplified_member:
            raise MemberNotFound(str(member_snowflake))
    else:
        if (
            member.voice
            and member.voice.channel
            and member.voice.channel.id == channel_snowflake
        ):
            is_channel_scope = True
        await channel.set_permissions(
            target=member,
            send_messages=False,
            add_reactions=False,
            reason=reason,
        )
    return is_channel_scope


def build_text_mute_embed(
    channel_snowflake: int,
    duration: DurationObject,
    guild_snowflake: int,
    member_snowflake: int,
    reason: str,
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    duration_builder = DurationBuilder()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise ChannelNotFound(str(channel_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        display_name = member.display_name
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(member_snowflake)
        if not simplified_member:
            raise MemberNotFound(str(member_snowflake))
        display_name = simplified_member[0]
        member_str = display_name
    embed = discord.Embed(
        title=f"{emojis.get_random_emoji()} {display_name} has been text-muted",
        description=(
            f"**User:** {member_str}\n"
            f"**Channel:** {channel.mention}\n"
            f"**Expires:** {duration_builder.load(duration=duration).to_unix_ts()}\n"
            f"**Reason:** {reason}"
        ),
        color=discord.Color.blue(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed
