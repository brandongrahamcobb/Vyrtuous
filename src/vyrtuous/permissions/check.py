"""check.py The purpose of this program is to provide the check utility module.

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

from vyrtuous.db.roles.admin.administrator_service import is_administrator
from vyrtuous.db.roles.coord.coordinator_service import is_coordinator
from vyrtuous.db.roles.dev.developer_service import is_developer
from vyrtuous.db.roles.mod.moderator_service import is_moderator
from vyrtuous.db.roles.owner.guild_owner_service import is_guild_owner
from vyrtuous.db.roles.sysadmin.sysadmin_service import is_sysadmin
from vyrtuous.inc.helpers import PERMISSION_TYPES
from vyrtuous.permissions.highest_role import resolve_highest_role


class HasEqualOrLowerRole(commands.CheckFailure):
    def __init__(self, target_rank=str):
        super().__init__(
            message=f"You may not execute this command on this `{target_rank}` because they have equal or higher role than you in this channel/server."
        )


async def check(
    snowflake_kwargs,
    lowest_role: str,
) -> str:
    verifications = (
        ("Sysadmin", is_sysadmin),
        ("Developer", is_developer),
        ("Guild Owner", is_guild_owner),
        ("Administrator", is_administrator),
        ("Coordinator", is_coordinator),
        ("Moderator", is_moderator),
    )
    channel_snowflake = snowflake_kwargs.get("channel_snowflake", None)
    guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
    member_snowflake = snowflake_kwargs.get("member_snowflake", None)
    passed_lowest = False
    for role_name, verify in verifications:
        try:
            if role_name in ("Sysadmin", "Developer"):
                if await verify(member_snowflake=int(member_snowflake)):
                    return role_name
            elif role_name in ("Guild Owner", "Administrator"):
                if await verify(
                    guild_snowflake=int(guild_snowflake),
                    member_snowflake=int(member_snowflake),
                ):
                    return role_name
            else:
                if await verify(
                    channel_snowflake=int(channel_snowflake),
                    guild_snowflake=int(guild_snowflake),
                    member_snowflake=int(member_snowflake),
                ):
                    return role_name
        except commands.CheckFailure:
            if lowest_role is not None and passed_lowest:
                raise
        if role_name == lowest_role:
            passed_lowest = True
    return "Everyone"


async def has_equal_or_lower_role_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
    member_snowflake: int,
    sender_snowflake: int,
) -> bool:
    snowflake_kwargs = {
        "channel_snowflake": source.channel.id,
        "guild_snowflake": source.guild.id,
        "member_snowflake": sender_snowflake,
    }
    return await has_equal_or_lower_role(
        snowflake_kwargs=snowflake_kwargs, member_snowflake=int(member_snowflake)
    )


async def has_equal_or_lower_role(
    snowflake_kwargs,
    member_snowflake: int,
) -> bool:
    sender_name = await resolve_highest_role(**snowflake_kwargs)
    sender_rank = PERMISSION_TYPES.index(sender_name)
    snowflake_kwargs.update({"member_snowflake": member_snowflake})
    where_kwargs = snowflake_kwargs
    target_name = await resolve_highest_role(**where_kwargs)
    target_rank = PERMISSION_TYPES.index(target_name)
    if sender_rank <= target_rank:
        raise HasEqualOrLowerRole(PERMISSION_TYPES[target_rank])
    return sender_name
