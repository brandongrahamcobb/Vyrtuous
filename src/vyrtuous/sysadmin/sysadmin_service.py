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


class NotSysadmin(commands.CommandError):
    def __init__(self, message="Member is not a sysadmin."):
        super().__init__(message)


class SysadminService:
    MODEL = Sysadmin

    def __init__(self, *, author_service=None, bot=None):
        self.__author_service = author_service
        self.__bot = bot

    async def is_sysadmin_wrapper(
        self,
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        member = self.__author_service.resolve_author(source=source)
        member_snowflake = member.id
        return await self.is_sysadmin(member_snowflake)

    def sysadmin_predicator(self):
        async def predicate(
            source: Union[commands.Context, discord.Interaction, discord.Message],
        ):
            if await self.is_sysadmin_wrapper(source):
                return True
            raise NotSysadmin

        predicate._permission_level = "Sysadmin"
        return commands.check(predicate)

    async def is_sysadmin(self, member_snowflake: int) -> bool:
        if int(self.__bot.config["discord_owner_id"]) == member_snowflake:
            return True
        raise NotSysadmin
