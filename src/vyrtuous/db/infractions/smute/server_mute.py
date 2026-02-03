"""server_mute.py The purpose of this program is to inherit from DatabaseFactory to provide the server mute moderation.

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

from vyrtuous.db.base.database_factory import DatabaseFactory
from vyrtuous.utils.dir_to_classes import skip_db_discovery


@skip_db_discovery
@dataclass(frozen=True)
class ServerMute(DatabaseFactory):

    __tablename__ = "active_server_voice_mutes"
    identifier = "smute"
    guild_snowflake: int
    member_snowflake: int
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    reason: str = "No reason provided."
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
