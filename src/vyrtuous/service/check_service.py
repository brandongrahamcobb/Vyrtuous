"""check_service.py The purpose of this program is to provide the check_service module.

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

from vyrtuous.database.roles.administrator import is_administrator_wrapper
from vyrtuous.database.roles.coordinator import is_coordinator_wrapper
from vyrtuous.database.roles.developer import is_developer_wrapper
from vyrtuous.database.roles.guild_owner import is_guild_owner_wrapper
from vyrtuous.database.roles.moderator import is_moderator_wrapper
from vyrtuous.database.roles.sysadmin import is_sysadmin_wrapper
from vyrtuous.utils.permission import PERMISSION_TYPES
from vyrtuous.database.roles.role import resolve_highest_role


class HasEqualOrLowerRole(commands.CheckFailure):
    def __init__(self, target_rank=str):
        super().__init__(
            message=f"You may not execute this command on this `{target_rank}` because they have equal or higher role than you in this channel/server."
        )


async def check(
    source: Union[commands.Context, discord.Interaction, discord.Message],
    lowest_role: str = None,
    member_snowflake: int = None,
) -> str:
    verifications = (
        ("SysAdmin", is_sysadmin_wrapper),
        ("Developer", is_developer_wrapper),
        ("Guild Owner", is_guild_owner_wrapper),
        ("Administrator", is_administrator_wrapper),
        ("Coordinator", is_coordinator_wrapper),
        ("Moderator", is_moderator_wrapper),
    )
    passed_lowest = False
    for role_name, verify in verifications:
        try:
            if await verify(source):
                return role_name
        except commands.CheckFailure:
            if lowest_role is not None and passed_lowest:
                raise
        if role_name == lowest_role:
            passed_lowest = True
    return "Everyone"


async def has_equal_or_lower_role(
    source: Union[commands.Context, discord.Interaction, discord.Message],
    member_snowflake: int,
    sender_snowflake: int,
) -> bool:
    kwargs = {
        "channel_snowflake": source.channel.id,
        "guild_snowflake": source.guild.id,
        "member_snowflake": sender_snowflake,
    }
    sender_name = await resolve_highest_role(**kwargs)
    sender_rank = PERMISSION_TYPES.index(sender_name)
    kwargs.update({"member_snowflake": member_snowflake})
    target_kwargs = kwargs
    target_name = await resolve_highest_role(**target_kwargs)
    target_rank = PERMISSION_TYPES.index(target_name)
    if sender_rank <= target_rank:
        raise HasEqualOrLowerRole(PERMISSION_TYPES[target_rank])

    return sender_name
