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


from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.author import resolve_author
from vyrtuous.commands.errors import NotSysadmin

from vyrtuous.base.record_service import RecordService


async def is_sysadmin_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    return await is_sysadmin(member_snowflake)


def sysadmin_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        if await is_sysadmin_wrapper(source):
            return True
        raise NotSysadmin

    predicate._permission_level = "Sysadmin"
    return commands.check(predicate)


async def is_sysadmin(member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    if int(bot.config["discord_owner_id"]) == member_snowflake:
        return True
    raise NotSysadmin


class SysadminService(RecordService):
    pass
