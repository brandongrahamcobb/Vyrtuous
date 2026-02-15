"""!/bin/python3
guild_owner_service.py The purpose of this program is to extend Service to service the guild owner class.

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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.owner.guild_owner import GuildOwner
from vyrtuous.sysadmin.sysadmin_service import SysadminService


class NotGuildOwner(commands.CommandError):
    def __init__(
        self,
        message="Member is not a guild owner in this server.",
    ):
        super().__init__(message)


class GuildOwnerService:
    MODEL = GuildOwner

    def __init__(self, *, author_service=None, bot=None, database_factory=None):
        self.__author_service = author_service
        self.__bot = bot
        self.__sysadmin_service = SysadminService(author_service=author_service)
        self.__developer_service = DeveloperService(
            author_service=author_service, database_factory=database_factory
        )
        self.__database_factory = database_factory

    async def update_guild_owners(self):
        for guild in self.__bot.guilds:
            guild_owner = await self.__database_factory.select(
                guild_snowflake=guild.id, singular=True
            )
            if guild_owner and guild_owner.member_snowflake == guild.owner_id:
                self.__bot.logger.info(
                    f"Guild owner ({guild_owner.member_snowflake}) already in the db."
                )
                continue
            else:
                guild_owner = self.MODEL(
                    guild_snowflake=guild.id, member_snowflake=guild.owner_id
                )
                await self.__database_factory.create(guild_owner)
                self.__bot.logger.info(
                    f"Guild owner ({guild_owner.member_snowflake}) added to the db."
                )

    async def is_guild_owner_wrapper(
        self,
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        member = self.__author_service.resolve_author(source=source)
        member_snowflake = member.id
        return await self.is_guild_owner(
            guild_snowflake=source.guild.id, member_snowflake=int(member_snowflake)
        )

    async def is_guild_owner_at_all(
        self,
        member_snowflake: int,
    ):
        for guild in self.__bot.guilds:
            if guild and guild.owner_id == member_snowflake:
                return True
        raise NotGuildOwner

    async def is_guild_owner(self, guild_snowflake: int, member_snowflake: int) -> bool:
        guild = self.__bot.get_guild(guild_snowflake)
        if guild and guild.owner_id == member_snowflake:
            return True
        raise NotGuildOwner
