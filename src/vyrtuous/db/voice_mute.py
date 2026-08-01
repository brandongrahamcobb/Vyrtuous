"""!/bin/python3
voice_mute.py The purpose of this program is to provide the voice-mute database model.

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
class VoiceMute:
    __tablename__ = "active_voice_mutes"
    identifier = "vmute"
    name = "Voice-Mute"
    guild_snowflake: int
    member_snowflake: int
    channel_snowflake: int | None = None
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    expires_in: datetime | None = None
    reason: str = "No reason provided."
    target: str = "click"
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
