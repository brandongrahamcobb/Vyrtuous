"""developer.py The purpose of this program is to inherit from the user class as a developer.

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

from datetime import datetime, timezone
from typing import Optional

from vyrtuous.database.database_factory import DatabaseFactory


class DeveloperLog(DatabaseFactory):

    ACT = "dlog"
    PLURAL = "Developer Logs"
    SINGULAR = "Developer Log"
    UNDO = "dlog"
    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "developer_snowflakes",
        "guild_snowflake",
        "id",
        "message_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "notes",
        "resolved",
        "updated_at",
    ]
    TABLE_NAME = "developer_logs"

    def __init__(
        self,
        channel_snowflake: Optional[int],
        developer_snowflakes: list[int],
        guild_snowflake: Optional[int],
        id: Optional[str],
        message_snowflake: Optional[int],
        created_at: Optional[datetime] = None,
        notes: Optional[str] = None,
        resolved: bool = False,
        updated_at: Optional[datetime] = None,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.created_at: datetime = created_at or datetime.now(timezone.utc)
        self.developer_snowflakes = developer_snowflakes
        self.id = id
        self.guild_snowflake = guild_snowflake
        self.message_snowflake = message_snowflake
        self.notes = notes
        self.resolved = resolved
        self.updated_at: datetime = updated_at or datetime.now(timezone.utc)
