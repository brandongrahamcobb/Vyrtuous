"""administrator.py The purpose of this program is to inherit from the DatabaseFactory to provide the administrator role.

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
from vyrtuous.db.roles.developer import is_developer_wrapper
from vyrtuous.db.roles.guild_owner import is_guild_owner_wrapper
from vyrtuous.db.roles.sysadmin import is_sysadmin_wrapper
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.dir_to_classes import skip_db_discovery


@skip_db_discovery
class NotAdministrator(commands.CheckFailure):
    def __init__(
        self,
        message="Member is not a sysadmin, developer, guild owner or administrator in this server.",
    ):
        super().__init__(message)


async def is_administrator_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    return await is_administrator(
        guild_snowflake=source.guild.id,
        member_snowflake=member_snowflake,
    )


async def is_administrator(guild_snowflake: int, member_snowflake: int) -> bool:
    administrator = await Administrator.select(
        guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
    )
    if not administrator:
        raise NotAdministrator
    return True


def administrator_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (
            is_sysadmin_wrapper,
            is_developer_wrapper,
            is_guild_owner_wrapper,
            is_administrator_wrapper,
        ):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure(
            "Member is not a sysadmin, developer, guild owner or administrator in this server."
        )

    predicate._permission_level = "Administrator"
    return commands.check(predicate)


class Administrator(DatabaseFactory):

    ACT = None
    CATEGORY = None
    PLURAL = "Administrators"
    SCOPES = [None]
    SINGULAR = "Administrator"
    UNDO = None

    REQUIRED_INSTANTIATION_ARGS = [
        "guild_snowflake",
        "member_snowflake",
        "role_snowflakes",
    ]
    OPTIONAL_ARGS = ["created_at", "updated_at"]

    TABLE_NAME = "administrators"

    def __init__(
        self,
        guild_snowflake: int,
        member_snowflake: int,
        role_snowflakes: list[int | None],
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
    ):
        super().__init__()
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.role_snowflakes = role_snowflakes
        self.updated_at = updated_at


class AdministratorRole(DatabaseFactory):

    ACT = "arole"
    CATEGORY = "arole"
    PLURAL = "Administrator Roles"
    SCOPES = ["guild"]
    SINGULAR = "Administrator Role"
    UNDO = "arole"

    REQUIRED_INSTANTIATION_ARGS = [
        "guild_snowflake",
        "role_snowflake",
    ]
    OPTIONAL_ARGS = ["created_at", "updated_at"]

    TABLE_NAME = "administrator_roles"

    def __init__(
        self,
        guild_snowflake: list[int | None],
        role_snowflake: list[int | None],
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.role_snowflake = role_snowflake
        self.updated_at = updated_at

    @classmethod
    async def grant_role(cls, guild_snowflake, role_snowflake):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        role = guild.get_role(role_snowflake)
        for member in role.members:
            try:
                await member.add_roles(role, reason="Granting administrator role.")
            except discord.Forbidden:
                continue

    @classmethod
    async def revoke_role(cls, guild_snowflake, role_snowflake):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        role = guild.get_role(role_snowflake)
        for member in role.members:
            try:
                await member.remove_roles(role, reason="Revoking a administrator role.")
            except discord.Forbidden:
                continue
