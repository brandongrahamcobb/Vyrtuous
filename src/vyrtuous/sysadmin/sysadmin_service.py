from copy import copy

"""!/bin/python3
sysadmin_service.py The purpose of this program is to extend Service to service the sysadmin class.

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

from vyrtuous.sysadmin.sysadmin import Sysadmin


class NotSysadmin(commands.CheckFailure):
    def __init__(self, message="Member is not a sysadmin."):
        super().__init__(message)


class SysadminService:
    MODEL = Sysadmin

    def __init__(self, *, author_service=None, bot=None, database_factory=None):
        self.__author_service = author_service
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL

    async def update_sysadmin(self):
        member_snowflake = self.__bot.config.get("discord_owner_id", None)
        sysadmin = await self.__database_factory.select(
            member_snowflake=int(member_snowflake), singular=True
        )
        if not sysadmin:
            sysadmin = Sysadmin(member_snowflake=int(member_snowflake))
            await self.__database_factory.create(sysadmin)
            self.__bot.logger.info(f"Sysadmin ({member_snowflake}) added to the db.")
        else:
            self.__bot.logger.info(f"Sysadmin ({member_snowflake}) already in the db.")

    async def is_sysadmin_wrapper(
        self,
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        member = self.__author_service.resolve_author(source=source)
        member_snowflake = member.id
        return await self.is_sysadmin(member_snowflake)

    async def is_sysadmin(self, member_snowflake: int) -> bool:
        if int(self.__bot.config["discord_owner_id"]) == member_snowflake:
            return True
        raise NotSysadmin
