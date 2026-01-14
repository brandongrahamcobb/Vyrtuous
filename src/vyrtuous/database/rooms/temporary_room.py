"""temporary_rooms.py A utility module for managing temporary rooms in the Vyrtuous Discord bot.

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
from vyrtuous.database.rooms.room import Room

class TemporaryRoom(Room):

    ACT = "temp"
    PLURAL = "Temporary Rooms"
    SINGULAR = "Temporary Rooms"
    UNDO = "temp"
    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "member_snowflake",
        "room_name",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "updated_at",
    ]
    TABLE_NAME = "temporary_rooms"

    def __init__(
        self,
        channel_snowflake: Optional[int],
        guild_snowflake: Optional[int],
        member_snowflake: Optional[int],
        room_name: Optional[str],
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None
    ):
        super().__init__()
        self.channel_mention = f"<#{channel_snowflake}>"
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.is_temp_room: Optional[bool] = True
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.member_snowflake = member_snowflake
        self.room_name = room_name
        self.updated_at = updated_at
