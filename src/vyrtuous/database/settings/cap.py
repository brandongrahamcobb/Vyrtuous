"""cap.py The purpose of this program is to provide the Cap utility class.

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

from vyrtuous.database.database_factory import DatabaseFactory


class Cap(DatabaseFactory):

    ACT = "cap"
    PLURAL = "Caps"
    SINGULAR = "Cap"
    UNDO = "cap"
    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "duration_seconds",
        "guild_snowflake",
        "moderation_type",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "updated_at",
    ]
    TABLE_NAME = "active_caps"

    def __init__(
        self,
        channel_snowflake: Optional[int],
        duration_seconds: Optional[int],
        guild_snowflake: Optional[int],
        moderation_type: Optional[str],
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
    ):
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.duration = duration_seconds
        self.guild_snowflake = guild_snowflake
        self.moderation_type = moderation_type
        self.updated_at = updated_at
