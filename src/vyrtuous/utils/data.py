"""!/bin/python3

data.py The purpose of this program is to manage statistics of Vyrtuous.

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

from dataclasses import dataclass
from datetime import datetime, timezone


@dataclass(frozen=True)
class Data:
    __table_name__ = "moderation_logs"
    current_channel_members: int = 0
    total_guild_members: int = 0
    online_members: int = 0
    total_voice_members: int = 0
    author_snowflake: int | None = None
    channel_snowflake: int | None = None
    guild_snowflake: int | None = None
    target_snowflake: int | None = None
    expires_at: datetime = datetime.now(timezone.utc)
    identiifer: str = ""
    reason: str = "No reason provided."
    is_modification: bool = False
    target_highest_role: str = "Role undetermined"
    executor_highest_role: str = "Role undetermined"
