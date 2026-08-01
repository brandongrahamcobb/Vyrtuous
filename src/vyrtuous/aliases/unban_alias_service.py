"""!/bin/python3
ban_service.py The purpose of this program is to extend AliasService to service ban infractions.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.db.ban import Ban
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import ChannelNotFound, GuildNotFound, MemberNotFound
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.moderation import cap_service
from vyrtuous.utils.tracking import stream_service

MODEL = Ban


@dataclass(frozen=True)
class UnbanMessageContext:
    author_snowflake: int
    channel_snowflake: int
    guild_snowflake: int
    member_snowflake: int
    message_snowflake: int
    message_channel_snowflake: int
    reason: str = "No reason provided."


async def unban_by_message(
    ctx: UnbanMessageContext,
    display: bool = True,
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    duration_builder = DurationBuilder()
    permission_state: PermissionState = bot.registry.get(PermissionState)
    ban = await database_factory.select(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        singular=True,
    )
    if not ban.blacklisted:
        duration = duration_builder.from_timestamp(ban.expires_in).build()
        exceeds_cap = await cap_service.exceeds_cap(
            category="ban",
            channel_snowflake=ctx.channel_snowflake,
            duration=duration,
            guild_snowflake=ctx.guild_snowflake,
        )
        if exceeds_cap:
            await permission_service.has_permissions(
                permission_state=permission_state,
                channel_snowflake=ctx.channel_snowflake,
                guild_snowflake=ctx.guild_snowflake,
                member_snowflake=ctx.author_snowflake,
                requested=["command.moderation.uncapped"],
            )
    else:
        return await build_blacklisted_block_embed(
            channel_snowflake=ctx.channel_snowflake,
            guild_snowflake=ctx.guild_snowflake,
            member_snowflake=ctx.member_snowflake,
        )
    is_channel_scope = await disable(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        reason=ctx.reason,
    )
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    await database_factory.delete(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
    )
    await stream_service.log(
        author_snowflake=ctx.author_snowflake,
        channel_snowflake=ctx.channel_snowflake,
        display=display,
        duration=DurationObject(number=0, prefix="", sign=1, unit=""),
        guild_snowflake=ctx.guild_snowflake,
        identifier="unban",
        is_channel_scope=is_channel_scope,
        member_snowflake=ctx.member_snowflake,
        message_snowflake=ctx.message_snowflake,
        message_channel_snowflake=ctx.message_channel_snowflake,
        reason=ctx.reason,
        role_snowflake=None,
        target=None,
    )
    return build_unban_embed(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
    )


async def disable(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int, reason: str
) -> bool:
    is_channel_scope = False
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise ChannelNotFound(str(channel_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        await channel.set_permissions(member, view_channel=None, reason=reason)
    else:
        simplified_member = bot.registry.get(MemberState).active.get(member_snowflake)
        if not simplified_member:
            raise MemberNotFound(str(member_snowflake))
    return is_channel_scope


def build_unban_embed(
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
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
        title=f"{emojis.get_random_emoji()} {display_name} has been unbanned",
        description=(f"**User:** {member_str}\n**Channel:** {channel.mention}"),
        color=discord.Color.yellow(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed


async def build_blacklisted_block_embed(
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
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
        title=f"{emojis.get_random_emoji()} {display_name} is blacklisted",
        description=f"**User:** {member_str}\n**Channel:** {channel.mention}\nUse {bot.config['discord_command_prefix']}blacklist to unblock.",
        color=discord.Color.yellow(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed
