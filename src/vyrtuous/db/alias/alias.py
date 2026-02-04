"""action.py The purpose of this program is to be a child of DatabaseFactory and the parent to all moderation actions.

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
class Alias(DatabaseFactory):

    __tablename__ = "command_aliases"
    identifier = "alias"
    alias_name: str
    category: str
    channel_snowflake: int
    guild_snowflake: int
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    role_snowflake: int | None = None
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))

    CATEGORY_TO_HELP = {
        "ban": [
            "**member**: Tag a member or include their ID",
            "**duration**: m/h/d",
            "**reason**: Reason for ban",
        ],
        "flag": [
            "**member**: Tag a member or include their ID",
            "**reason**: Reason for flag",
        ],
        "role": [
            "**member**: Tag a member or include their ID",
            "**role**: Tag a role or include its ID",
        ],
        "tmute": [
            "**member**: Tag a member or include their ID",
            "**duration**: m/h/d",
            "**reason**: Reason for text-mute",
        ],
        "vegan": ["**member**: Tag a member or include their ID"],
        "vmute": [
            "**member**: Tag a member or include their ID",
            "**duration**: m/h/d",
            "**reason**: Reason for voice-mute",
        ],
    }
    CATEGORY_TO_DESCRIPTION = {
        "ban": "Toggles a ban.",
        "flag": "Toggles a moderation flag.",
        "role": "Toggles a role to a user.",
        "tmute": "Toggles a mute in text channels.",
        "vegan": "Toggles a going vegan flag.",
        "vmute": "Toggles a mute in voice channels.",
    }
    CATEGORY_TO_PERMISSION_LEVEL = {
        "ban": "Moderator",
        "flag": "Moderator",
        "role": "Coordinator",
        "tmute": "Moderator",
        "vegan": "Moderator",
        "vmute": "Moderator",
    }
