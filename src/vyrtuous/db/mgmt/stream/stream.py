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

from dataclasses import dataclass, field
from datetime import datetime, timezone

from vyrtuous.base.database_factory import DatabaseFactory


@dataclass(frozen=True)
class Stream(DatabaseFactory):

    __tablename__ = "streaming"
    identifier = "stream"
    channel_snowflake: int
    enabled: bool
    entry_type: str
    guild_snowflake: int
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    snowflakes: list[int | None] = field(default_factory=list)
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
