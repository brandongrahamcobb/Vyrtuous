"""discord_object_service.py The purpose of this program is to provide the DiscordObject module.

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

import re
from typing import Any, Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.discord.source import DiscordSourceNotFound
from vyrtuous.utils.logger import logger


class DiscordObjectNotFound(commands.CheckFailure):
    "Returns an error if a channel, guild, member or role is not found."

    def __init__(self, target: str, *, message: str | None = None):
        super().__init__(
            message=message
            or f"Unable to resolve a valid channel, guild, member or role for target (`{target}`)."
        )


class GuildChannelNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            message=f"Unable to resolve a valid Discord guild channel with the provided context and target (`{target}`).",
            target=target,
        )


class GuildNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            message=f"Unable to resolve a valid Discord guild with the provided context and target (`{target}`).",
            target=target,
        )


class GuildMemberNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            message=f"Unable to resolve a valid Discord guild member with the provided context and target (`{target}`).",
            target=target,
        )


class GuildRoleNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            message=f"Unable to resolve a valid Discord guild role with the provided context and target (`{target}`).",
            target=target,
        )


class TargetIsBot(commands.CheckFailure):
    def __init__(
        self,
        *,
        ctx: commands.Context | None = None,
        interaction: discord.Interaction | None = None,
        message: discord.Message | None = None,
    ):
        if (ctx is None) == (interaction is None) == (message is None):
            raise DiscordSourceNotFound()
        self._source = ctx or interaction or message
        super().__init__(
            message=f"You cannot execute actions on {self._source.guild.me.mention}."
        )


class DiscordObject:

    def __init__(
        self,
        *,
        ctx: commands.Context | None = None,
        interaction: discord.Interaction | None = None,
        message: discord.Message | None = None,
    ):
        if (ctx is None) == (interaction is None) == (message is None):
            raise DiscordSourceNotFound()
        self.bot = DiscordBot.get_instance()
        self._source = ctx or interaction or message

    async def determine_from_target(self, target) -> dict[str, Any]:
        if target is None:
            raise TypeError()
        if isinstance(self._source, (commands.Context, discord.Message)):
            do = DiscordObject(ctx=self._source)
            author = self._source.author
        elif isinstance(self._source, discord.Interaction):
            do = DiscordObject(interaction=self._source)
            author = self._source.user
        if target and str(target).lower() == "all":
            return {"columns": {}}
        try:
            channel = do.resolve_channel(target)
        except GuildChannelNotFound as e:
            logger.warning(e)
        else:
            return {
                "author_snowflake": author.id,
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
            if member == self._source.guild.me:
                raise TargetIsBot()
            return {
                "author_snowflake": author.id,
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
        try:
            guild = do.resolve_guild(target)
        except GuildNotFound as e:
            logger.warning(e)
        else:
            return {
                "author_snowflake": author.id,
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
                "author_snowflake": author.id,
                "columns": {
                    "guild_snowflake": role.guild.id,
                    "role_snowflake": role.id,
                },
                "id": role.id,
                "mention": role.mention,
                "name": role.name,
                "type": type(role),
                "object": role,
            }
        raise DiscordObjectNotFound(target=target)

    def resolve_channel(
        self,
        target: Union[int, str],
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
        raise GuildChannelNotFound(target=str(target))

    def resolve_guild(
        self,
        target: Union[int, str],
    ) -> discord.Guild:
        g_id = None
        if isinstance(target, int):
            g_id = target
        elif isinstance(target, str):
            if target.isdigit():
                g_id = int(target)
            elif re_match := re.match(r"^<#(\d+)>$", target):
                g_id = int(re_match.group(1))
        if g_id:
            g = self.bot.get_guild(g_id)
            if isinstance(g, (discord.Guild)):
                return g
        raise GuildNotFound(target=str(target))

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
        raise GuildMemberNotFound(target=str(target))

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
            if isinstance(r, discord.Role):
                return r
        raise GuildRoleNotFound(target=str(target))
