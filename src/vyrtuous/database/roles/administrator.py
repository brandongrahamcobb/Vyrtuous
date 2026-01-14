"""administrator.py The purpose of this program is to inherit from the PermissionRole to provide the administrator role.

Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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
from typing import Optional
from vyrtuous.database.roles.permission_role import PermissionRole

class Administrator(PermissionRole):

    ACT = None
    PLURAL = "Administrators"
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
        guild_snowflake: Optional[int],
        member_snowflake: Optional[int],
        role_snowflakes: list[int | None],
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None
    ):
        super().__init__()
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.role_snowflakes = role_snowflakes
        self.updated_at = updated_at

class AdministratorRole(PermissionRole):

    REQUIRED_INSTANTIATION_ARGS = [
        "guild_snowflake",
        "role_snowflake",
    ]
    OPTIONAL_ARGS = ["created_at", "updated_at"]

    def __init__(
        self,
        guild_snowflake: list[int | None],
        role_snowflake: list[int | None],
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None
    ):
        super().__init__()
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.role_snowflake = role_snowflake
        self.updated_at = updated_at