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

from datetime import datetime, timezone

from vyrtuous.db.database_factory import DatabaseFactory


class Sysadmin(DatabaseFactory):

    __tablename__ = "sysadmin"
    category = "sysadmin"
    created_at: datetime
    member_snowflake: int
    updated_at: datetime

    def __init__(
        self,
        member_snowflake: int,
        created_at: datetime | None = None,
        updated_at: datetime | None = None,
    ):
        self.created_at = created_at or datetime.now(timezone.utc)
        self.member_snowflake = member_snowflake
        self.updated_at = updated_at or datetime.now(timezone.utc)
