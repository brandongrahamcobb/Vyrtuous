"""cap.py The purpose of this program is to provide the Cap utility class.

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
from typing import Optional

from vyrtuous.db.database_factory import DatabaseFactory


class Cap(DatabaseFactory):

    ACT = None
    CATEGORY = "cap"
    PLURAL = "Caps"
    SCOPES = ["channels"]
    SINGULAR = "Cap"
    UNDO = None

    REQUIRED_INSTANTIATION_ARGS = [
        "category",
        "channel_snowflake",
        "duration_seconds",
        "guild_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "updated_at",
    ]

    TABLE_NAME = "active_caps"

    def __init__(
        self,
        category: str,
        channel_snowflake: int,
        duration_seconds: int,
        guild_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
    ):
        self.category = category
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.duration_seconds = duration_seconds
        self.guild_snowflake = guild_snowflake
        self.updated_at = updated_at
