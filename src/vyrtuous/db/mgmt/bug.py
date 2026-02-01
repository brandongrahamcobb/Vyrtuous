"""developer.py The purpose of this program is to inherit from the user class as a developer.

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


class Bug(DatabaseFactory):

    __tablename__ = "bug_tracking"
    channel_snowflake: int
    guild_snowflake: int
    id: str
    member_snowflakes: list[int]
    message_snowflake: int
    notes: str
    resolved: bool
    created_at: datetime
    updated_at: datetime

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        id: str,
        member_snowflakes: list[int],
        message_snowflake: int,
        notes: str | None = None,
        resolved: bool = False,
        created_at: datetime | None = None,
        updated_at: datetime | None = None,
    ):
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at or datetime.now(timezone.utc)
        self.member_snowflakes = member_snowflakes
        self.id = id
        self.guild_snowflake = guild_snowflake
        self.message_snowflake = message_snowflake
        self.notes = notes
        self.resolved = resolved
        self.updated_at = updated_at or datetime.now(timezone.utc)
