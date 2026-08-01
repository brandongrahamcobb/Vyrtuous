"""!/bin/python3
ban.py The purpose of this program is to provide the ban database model.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone


@dataclass(frozen=True)
class Ban:
    __tablename__ = "active_bans"
    identifier = "ban"
    name = "Ban"
    channel_snowflake: int
    guild_snowflake: int
    member_snowflake: int
    blacklisted: bool = False
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    expires_in: datetime | None = None
    last_kicked: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    reason: str = "No reason provided."
    reset: bool = False
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
