"""flag.py The purpose of this program is to inherit from Action to provide the flag moderation.

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
from vyrtuous.database.actions.action import Action


class Flag(Action):

    ACT = "flag"
    PLURAL = "Flags"
    SINGULAR = "Flag"
    UNDO = "unflag"
    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "member_snowflake",
    ]
    OPTIONAL_ARGS = ["created_at", "reason", "updated_at"]
    TABLE_NAME = "active_flags"

    def __init__(
        self,
        channel_snowflake: Optional[int],
        guild_snowflake: Optional[int],
        member_snowflake: Optional[int],
        created_at: Optional[datetime] = None,
        reason: Optional[str] = "No reason provided.",
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>" if channel_snowflake else None
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.reason = reason
        self.updated_at = updated_at
