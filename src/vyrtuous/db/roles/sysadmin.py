"""sysadmin.py The purpose of this program is to inherit from the DatabaseFactory to provide the sysadmin role.

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

from datetime import datetime
from typing import Optional, Union

from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.dir_to_classes import skip_db_discovery

@skip_db_discovery
class NotSysadmin(commands.CheckFailure):
    def __init__(self, message="Member is not a sysadmin."):
        super().__init__(message)


async def is_sysadmin_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    return is_sysadmin(member_snowflake)


def sysadmin_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        if await is_sysadmin_wrapper(source):
            return True
        raise commands.CheckFailure("Member is not a sysadmin.")

    predicate._permission_level = "Sysadmin"
    return commands.check(predicate)


def is_sysadmin(member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    if int(bot.config["discord_owner_id"]) == member_snowflake:
        return True
    raise NotSysadmin


class Sysadmin(DatabaseFactory):

    ACT = None
    CATEGORY = None
    PLURAL = "Sysadmin"
    SCOPES = [None]
    SINGULAR = "Sysadmin"
    UNDO = None

    REQUIRED_INSTANTIATION_ARGS = ["member_snowflake"]
    OPTIONAL_ARGS = ["created_at", "updated_at"]
    
    TABLE_NAME = "sysadmin"

    def __init__(
        self,
        member_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.created_at = created_at
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.updated_at = updated_at    