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

from dataclasses import dataclass, field
from datetime import datetime, timezone

from vyrtuous.base.database_factory import DatabaseFactory


@dataclass(frozen=True)
class Administrator(DatabaseFactory):

    __tablename__ = "administrators"
    identifier = "admin"
    guild_snowflake: int
    member_snowflake: int
    role_snowflakes: list[int | None] = field(default_factory=list)
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))


@dataclass(frozen=True)
class AdministratorRole(DatabaseFactory):

    __tablename__ = "administrator_roles"
    identifier = "arole"
    guild_snowflake: int
    role_snowflake: int
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
