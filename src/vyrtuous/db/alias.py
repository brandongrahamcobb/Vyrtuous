"""!/bin/python3
alias.py The purpose of this program is to extend DatabaseFactory for text-based specific aliases.

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

from discord.ext import commands


@dataclass(frozen=True)
class Alias:
    __tablename__ = "command_aliases"
    identifier = "alias"
    alias_name: str
    category: str
    channel_snowflake: int
    guild_snowflake: int
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    role_snowflake: int | None = None
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))


class NotAlias(commands.CheckFailure):
    def __init__(self):
        super().__init__(
            message=("Invalid alias."),
        )
