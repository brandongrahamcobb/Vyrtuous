"""voice_mute.py The purpose of this program is to inherit from DatabaseFactory to provide the voice mute moderation.

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


class VoiceMute(DatabaseFactory):

    __tablename__ = "active_voice_mutes"
    category = "vmute"
    channel_snowflake: int
    created_at: datetime
    expires_in: datetime
    guild_snowflake: int
    member_snowflake: int
    reason: str
    updated_at: datetime

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
        created_at: datetime | None = None,
        expires_in: datetime | None = None,
        reason: str = "No reason provided.",
        target: str = "user",
        updated_at: datetime | None = None,
    ):
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at or datetime.now(timezone.utc)
        self.expires_in = expires_in or datetime.now(timezone.utc)
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.reason = reason
        self.target = target
        self.updated_at = updated_at or datetime.now(timezone.utc)

    @property
    def target(self):
        return self._target

    @target.setter
    def target(self, target):
        if target not in ["room", "user"]:
            raise ValueError("Invalid target.")
        self._target = target

    @property
    def expired(self) -> bool:
        return datetime.now(timezone.utc) > self.expires_in
