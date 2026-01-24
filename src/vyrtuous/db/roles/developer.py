"""developer.py The purpose of this program is to inherit from the DatabaseFactory to provide the developer role.

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

from vyrtuous.db.roles.sysadmin import is_sysadmin_wrapper
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.dir_to_classes import skip_db_discovery

@skip_db_discovery
class NotDeveloper(commands.CheckFailure):
    def __init__(self, message="Member is not a sysadmin or developer."):
        super().__init__(message)


async def is_developer_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    return await is_developer(member_snowflake)


def developer_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (is_sysadmin_wrapper, is_developer_wrapper):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure("Member is not a sysadmin or developer.")

    predicate._permission_level = "Developer"
    return commands.check(predicate)


async def is_developer(member_snowflake: int) -> bool:
    developer = await Developer.select(member_snowflake=member_snowflake)
    if not developer:
        raise NotDeveloper
    return True


class Developer(DatabaseFactory):

    ACT = "dev"
    CATEGORY = "dev"
    PLURAL = "Developers"
    SCOPES = ["member"]
    SINGULAR = "Developer"
    UNDO = "dev"

    REQUIRED_INSTANTIATION_ARGS = ["guild_snowflake", "member_snowflake"]
    OPTIONAL_ARGS = ["created_at", "updated_at"]

    TABLE_NAME = "developers"

    def __init__(
        self,
        guild_snowflake: int,
        member_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.updated_at = updated_at
