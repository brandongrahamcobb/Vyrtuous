"""role.py The purpose of this program is to be a child of DatabaseFactory and the parent to all permission roles.

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

from vyrtuous.database.roles.administrator import NotAdministrator, is_administrator
from vyrtuous.database.roles.coordinator import NotCoordinator, is_coordinator
from vyrtuous.database.roles.developer import NotDeveloper, is_developer
from vyrtuous.database.roles.guild_owner import NotGuildOwner, is_guild_owner
from vyrtuous.database.roles.moderator import NotModerator, is_moderator
from vyrtuous.database.roles.sysadmin import NotSysAdmin, is_sysadmin
from vyrtuous.service.logging_service import logger


async def resolve_highest_role(
    channel_snowflake: int = None,
    member_snowflake: int = None,
    guild_snowflake: int = None,
):
    try:
        if is_sysadmin(member_snowflake=member_snowflake):
            return "SysAdmin"
    except NotSysAdmin as e:
        logger.warning(str(e).capitalize())
    try:
        if await is_developer(member_snowflake=member_snowflake):
            return "Developer"
    except NotDeveloper as e:
        logger.warning(str(e).capitalize())
    try:
        if await is_guild_owner(
            guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
        ):
            return "Guild Owner"
    except NotGuildOwner as e:
        logger.warning(str(e).capitalize())
    try:
        if await is_administrator(
            guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
        ):
            return "Administrator"
    except NotAdministrator as e:
        logger.warning(str(e).capitalize())
    if channel_snowflake:
        try:
            if await is_coordinator(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                member_snowflake=member_snowflake,
            ):
                return "Coordinator"
        except NotCoordinator as e:
            logger.warning(str(e).capitalize())
        try:
            if await is_moderator(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                member_snowflake=member_snowflake,
            ):
                return "Moderator"
        except NotModerator as e:
            logger.warning(str(e).capitalize())
    return "Everyone"
