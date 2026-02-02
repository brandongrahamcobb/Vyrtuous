"""streaming.py A utility module for managing and streaming of messages in the Vyrtuous Discord bot.

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


class Stream(DatabaseFactory):

    __tablename__ = "streaming"
    category = "stream"
    channel_snowflake: int
    enabled: bool
    entry_type: str
    guild_snowflake: int
    created_at: datetime
    snowflakes: list[int | None]
    updated_at: datetime

    TABLE_NAME = "streaming"
    lines, pages = [], []

    def __init__(
        self,
        channel_snowflake: int,
        enabled: bool,
        entry_type: str,
        guild_snowflake: int,
        created_at: datetime | None = None,
        snowflakes: list[int | None] = list[None],
        updated_at: datetime | None = None,
    ):
        self._action: str
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at or datetime.now(timezone.utc)
        self.enabled = enabled
        self.entry_type = entry_type
        self.guild_snowflake = guild_snowflake
        self.snowflakes = snowflakes
        self.updated_at = updated_at or datetime.now(timezone.utc)

    @property
    def action(self):
        return self._action

    @action.setter
    def action(self, action: str):
        if action not in ("create", "modify", "delete"):
            raise ValueError("Invalid action.")
        self._action = action

    @property
    def entry_type(self):
        return self._entry_type

    @entry_type.setter
    def entry_type(self, entry_type: str):
        if entry_type not in ("all", "channel"):
            raise ValueError("Invalid entry type.")
        self._entry_type = entry_type
