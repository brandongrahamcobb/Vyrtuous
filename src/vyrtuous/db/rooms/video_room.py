"""video_rooms.py A utility module for managing video rooms in the Vyrtuous Discord bot.

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


class VideoRoom(DatabaseFactory):

    __tablename__ = "video_rooms"
    category = "vroom"
    channel_snowflake: int
    guild_snowflake: int
    created_at: datetime
    updated_at: datetime

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        created_at: datetime | None = None,
        updated_at: datetime | None = None,
    ):
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at or datetime.now(timezone.utc)
        self.guild_snowflake = guild_snowflake
        self.is_video_room: bool = True
        self.updated_at = updated_at or datetime.now(timezone.utc)
