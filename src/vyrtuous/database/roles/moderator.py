"""moderator.py The purpose of this program is to inherit from the DatabaseFactory to provide the moderator role.

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

from vyrtuous.database.database_factory import DatabaseFactory
from vyrtuous.database.roles.administrator import is_administrator_wrapper
from vyrtuous.database.roles.coordinator import is_coordinator_wrapper
from vyrtuous.database.roles.developer import is_developer_wrapper
from vyrtuous.database.roles.guild_owner import is_guild_owner_wrapper
from vyrtuous.database.roles.sysadmin import is_sysadmin_wrapper
from vyrtuous.service.member_snowflake import get_member_snowflake


class NotModerator(commands.CheckFailure):
    def __init__(
        self,
        message="Member is not a sysadmin, developer, guild owner, administrator, coordinator or moderator in this channel.",
    ):
        super().__init__(message)


async def is_moderator_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member_snowflake = get_member_snowflake(source=source)
    return await is_moderator(
        channel_snowflake=source.channel.id,
        guild_snowflake=source.guild.id,
        member_snowflake=member_snowflake,
    )


async def is_moderator(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> bool:
    moderator = await Moderator.select(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
    )
    if not moderator:
        raise NotModerator
    return True


def moderator_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (
            is_sysadmin_wrapper,
            is_developer_wrapper,
            is_guild_owner_wrapper,
            is_administrator_wrapper,
            is_coordinator_wrapper,
            is_moderator_wrapper,
        ):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure(
            "Member is not a sysadmin, developer, guild owner, administrator, coordinator or moderator in this channel."
        )

    predicate._permission_level = "Moderator"
    return commands.check(predicate)


class Moderator(DatabaseFactory):

    ACT = "mod"
    PLURAL = "Moderators"
    SINGULAR = "Moderator"
    UNDO = "mod"
    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "member_snowflake",
    ]
    OPTIONAL_ARGS = ["created_at", "updated_at"]
    TABLE_NAME = "moderators"

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>" if channel_snowflake else None
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.updated_at = updated_at
