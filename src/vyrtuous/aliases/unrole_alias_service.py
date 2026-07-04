"""!/bin/python3
role_service.py The purpose of this program is to extend AliasService to service the role class.

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
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.tracking import data_builder, stream_service


@dataclass(frozen=True)
class UnroleMessageContext:
    author_snowflake: int
    guild_snowflake: int
    member_snowflake: int
    message_snowflake: int
    message_channel_snowflake: int
    role_snowflake: int


async def set_unrole_overwrite(
    guild_snowflake: int, member_snowflake: int, role_snowflake: int
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    role = guild.get_role(role_snowflake)
    if role is None:
        raise commands.RoleNotFound(str(role_snowflake))
    member = guild.get_member(member_snowflake)
    if member is None:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member is None:
            raise commands.MemberNotFound(str(member_snowflake))
    else:
        try:
            await member.remove_roles(role, reason="Revoking role.")
        except discord.Forbidden:
            raise


def build_unrole_embed(
    guild_snowflake: int,
    member_snowflake: int,
    role_snowflake: int,
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    role = guild.get_role(role_snowflake)
    if role is None:
        raise commands.RoleNotFound(str(role_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        display_name = member.display_name
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member is None:
            raise commands.MemberNotFound(str(member_snowflake))
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


async def log_unrole(
    author_snowflake: int | None,
    display: bool,
    guild_snowflake: int,
    member_snowflake: int,
    message_snowflake: int | None,
    message_channel_snowflake: int | None,
    role_snowflake: int,
) -> None:
    channel_snowflake = None
    duration_value = None
    is_channel_scope = None
    reason = None
    target = None
    await data_builder.save_data(
        author_snowflake=author_snowflake or None,
        channel_snowflake=channel_snowflake or None,
        duration_value=duration_value or None,
        guild_snowflake=guild_snowflake,
        identifier="unrole",
        member_snowflake=member_snowflake,
        reason=reason or "No reason provided.",
        role_snowflake=role_snowflake or None,
        target=target or None,
    )
    if display:
        await stream_service.send_log(
            author_snowflake=author_snowflake or None,
            channel_snowflake=channel_snowflake or None,
            identifier="unrole",
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


async def unrole_by_message(
    ctx: UnroleMessageContext, display: bool = True
) -> discord.Embed:
    await set_unrole_overwrite(
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        role_snowflake=ctx.role_snowflake,
    )
    await log_unrole(
        author_snowflake=ctx.author_snowflake,
        display=display,
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        message_snowflake=ctx.message_snowflake,
        message_channel_snowflake=ctx.message_channel_snowflake,
        role_snowflake=ctx.role_snowflake,
    )
    embed = build_unrole_embed(
        guild_snowflake=ctx.guild_snowflake,
        member_snowflake=ctx.member_snowflake,
        role_snowflake=ctx.role_snowflake,
    )
    return embed
