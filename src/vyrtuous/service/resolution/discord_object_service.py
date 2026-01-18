"""discord_object_service.py The purpose of this program is to provide the DiscordObject module.

Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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

import re
from typing import Optional, Union

from discord.ext import commands
import discord

from vyrtuous.service.check_service import resolve_highest_permission_role
from vyrtuous.service.logging_service import logger
from vyrtuous.utils.permission import PERMISSION_TYPES


def determine_from_target(source, target) -> dict:
    if isinstance(source, commands.Context):
        do = DiscordObject(ctx=source)
        author = source.author
    elif isinstance(source, discord.Interaction):
        do = DiscordObject(interaction=source)
        author = source.user
    try:
        channel = do.resolve_channel(target)
    except GuildChannelNotFound as e:
        logger.warning(e)
    else:
        return {
            "columns": {
                "channel_snowflake": channel.id,
                "guild_snowflake": channel.guild.id,
            },
            "id": channel.id,
            "mention": channel.mention,
            "name": channel.name,
            "type": type(channel),
            "object": channel,
        }
    try:
        member = do.resolve_member(target)
    except GuildMemberNotFound as e:
        logger.warning(e)
    else:
        if member == source.guild.me:
            raise TargetIsBot()
        await has_equal_or_lower_role(
            source=source, member_snowflake=member.id, sender_snowflake=author.id
        )
        return {
            "columns": {
                "member_snowflake": member.id,
                "guild_snowflake": member.guild.id,
            },
            "id": member.id,
            "mention": member.mention,
            "name": member.display_name,
            "type": type(member),
            "object": member,
        }
    guild = source.bot.get_guild(target)
    if guild:
        return {
            "columns": {"guild_snowflake": guild.id},
            "id": guild.id,
            "name": guild.name,
            "type": type(guild),
            "object": guild,
        }
    try:
        role = do.resolve_role(target)
    except GuildRoleNotFound as e:
        logger.warning(e)
    else:
        return {
            "columns": {"guild_snowflake": role.guild.id, "role_snowflake": role.id},
            "id": role.id,
            "mention": role.mention,
            "name": role.name,
            "type": type(role),
            "object": role,
        }
    raise DiscordObjectNotFound()


async def has_equal_or_lower_role(
    source, member_snowflake: int, sender_snowflake: int
) -> bool:
    kwargs = {
        "channel_snowflake": source.channel.id,
        "guild_snowflake": source.guild.id,
        "member_snowflake": sender_snowflake,
    }
    sender_rank = await resolve_highest_permission_role(**kwargs)

    kwargs.update({"member_snowflake": member_snowflake})
    target_kwargs = kwargs
    target_rank = await resolve_highest_permission_role(**target_kwargs)

    if sender_rank <= target_rank:
        raise HasEqualOrLowerRole(PERMISSION_TYPES[target_rank])

    return False


class DiscordObjectNotFound(Exception):
    "Returns an error of a channel, guild, member or role is not found."


class DiscordSourceNotFound(Exception):

    def __init__(self):
        super().__init__(
            "Unable to resolve a valid Discord object due to missing ctx (commands.Context) or interaction (discord.Interaction)."
        )


class GuildChannelNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            f"Unable to resolve a valid Discord guild channel with the provided context and target ({target})."
        )


class GuildMemberNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            f"Unable to resolve a valid Discord guild member with the provided context and target ({target})."
        )


class GuildRoleNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            f"Unable to resolve a valid Discord guild role with the provided context and target ({target})."
        )


class TargetIsBot(commands.CheckFailure):
    def __init__(
        self, source: Union[commands.Context, discord.Interaction, discord.Message]
    ):
        super().__init__(
            message="You cannot execute actions on {source.guild.me.mention}."
        )


class HasEqualOrLowerRole(commands.CheckFailure):
    def __init__(self, target_rank=str):
        super().__init__(
            message=f"You may not execute this command on this `{PERMISSION_TYPES[target_rank]}` because they have equal or higher role than you in this channel/server."
        )


class DiscordObject:

    def __init__(
        self, *, ctx: commands.Context = None, interaction: discord.Interaction = None
    ):
        if (ctx is None) == (interaction is None):
            raise DiscordSourceNotFound()
        self._source = ctx or interaction

    def resolve_channel(
        self,
        target: Optional[Union[int, str]],
    ) -> Union[discord.TextChannel, discord.VoiceChannel]:
        c_id = None
        if isinstance(target, int):
            c_id = target
        elif isinstance(target, str):
            if target.isdigit():
                c_id = int(target)
            elif re_match := re.match(r"^<#(\d+)>$", target):
                c_id = int(re_match.group(1))
        if c_id:
            c = self._source.guild.get_channel(c_id)
            if isinstance(c, (discord.TextChannel, discord.VoiceChannel)):
                return c
        raise GuildChannelNotFound(target=target)

    def resolve_member(self, target) -> discord.Member:
        m_id = None
        if isinstance(target, int):
            m_id = target
        elif isinstance(target, str):
            if target.isdigit():
                m_id = int(target)
            elif re_match := re.match(r"^<@(\d+)>$", target):
                m_id = int(re_match.group(1))
        if m_id:
            m = self._source.guild.get_member(m_id)
            if isinstance(m, discord.Member):
                return m
        raise GuildMemberNotFound(target=target)

    def resolve_role(self, target) -> discord.Role:
        r_id = None
        if isinstance(target, int):
            r_id = target
        elif isinstance(target, str):
            if target.isdigit():
                r_id = int(target)
            elif re_match := re.match(r"^<@&(\d+)>$", target):
                r_id = int(re_match.group(1))
        if r_id:
            r = self._source.guild.get_role(r_id)
            print(r)
            if isinstance(r, discord.Role):
                return r
        raise GuildRoleNotFound(target=target)
