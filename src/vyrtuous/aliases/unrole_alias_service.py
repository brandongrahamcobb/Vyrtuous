"""!/bin/python3
role_service.py The purpose of this program is to extend AliasService to service the role class.

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
from vyrtuous.cache.registry import MemberState
from vyrtuous.models.duration import DurationObject
from vyrtuous.utils.errors.error import GuildNotFound, MemberNotFound, RoleNotFound
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.tracking import stream_service


@dataclass(frozen=True)
class UnroleMessageContext:
    author_snowflake: int
    channel_snowflake: int
    guild_snowflake: int
    member_snowflake: int
    message_snowflake: int
    message_channel_snowflake: int
    role_snowflake: int


async def disable(
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
    role_snowflake: int,
    reason: str = "No reason provided.",
) -> bool:
    is_channel_scope = False
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    role = guild.get_role(role_snowflake)
    if role is None:
        raise RoleNotFound(str(role_snowflake))
    member = guild.get_member(member_snowflake)
    if member is None:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member is None:
            raise MemberNotFound(str(member_snowflake))
    else:
        if (
            member.voice
            and member.voice.channel
            and member.voice.channel.id == channel_snowflake
        ):
            is_channel_scope = True
        try:
            await member.remove_roles(role, reason=reason)
        except discord.Forbidden:
            raise
    return is_channel_scope


def build_unrole_embed(
    guild_snowflake: int,
    member_snowflake: int,
    role_snowflake: int,
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    role = guild.get_role(role_snowflake)
    if role is None:
        raise RoleNotFound(str(role_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        display_name = member.display_name
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member is None:
            raise MemberNotFound(str(member_snowflake))
        else:
            display_name = simplified_member[0]
            member_str = display_name
    embed = discord.Embed(
        title=f"{emojis.get_random_emoji()} {display_name}'s role has been revoked",
        description=(f"**User:** {member_str}\n" f"**Role:** {role.mention}"),
        color=discord.Color.yellow(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed


async def unrole_by_message(
    ctx: UnroleMessageContext, display: bool = True
) -> discord.Embed:
    is_channel_scope = await disable(
        channel_snowflake=ctx.channel_snowflake,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        role_snowflake=ctx.role_snowflake,
    )
    await stream_service.log(
        author_snowflake=ctx.author_snowflake,
        channel_snowflake=ctx.channel_snowflake,
        display=display,
        duration=DurationObject(number=0, prefix="", sign=1, unit=""),
        guild_snowflake=ctx.guild_snowflake,
        identifier="unrole",
        is_channel_scope=is_channel_scope,
        member_snowflake=ctx.member_snowflake,
        message_snowflake=ctx.message_snowflake,
        message_channel_snowflake=ctx.message_channel_snowflake,
        reason="No reason provided.",
        role_snowflake=ctx.role_snowflake,
        target=None,
    )
    embed = build_unrole_embed(
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        role_snowflake=ctx.role_snowflake,
    )
    return embed
