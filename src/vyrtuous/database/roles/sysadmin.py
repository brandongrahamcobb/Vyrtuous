"""sysadmin.py The purpose of this program is to manage the sysadmin.

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

from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.member_snowflake import get_member_snowflake


class NotSysAdmin(commands.CheckFailure):
    def __init__(self, message="Member is not a sysadmin."):
        super().__init__(message)


async def is_sysadmin_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member_snowflake = get_member_snowflake(source=source)
    return is_sysadmin(member_snowflake)


def sysadmin_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        if await is_sysadmin_wrapper(source):
            return True
        raise commands.CheckFailure("Member is not a sysadmin.")

    predicate._permission_level = "SysAdmin"
    return commands.check(predicate)


def is_sysadmin(member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    if int(bot.config["discord_owner_id"]) == member_snowflake:
        return True
    raise NotSysAdmin
