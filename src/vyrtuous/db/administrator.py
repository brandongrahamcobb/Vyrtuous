"""!/bin/python3
administrator.py The purpose of this program is to extend DatabaseFactory provide the administrator and administrator role classes.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot


@dataclass(frozen=True)
class Administrator:

    __tablename__ = "administrators"
    identifier = "admin"
    guild_snowflake: int
    member_snowflake: int
    role_snowflakes: list[int | None] = field(default_factory=list)
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))


@dataclass(frozen=True)
class AdministratorRole:

    __tablename__ = "administrator_roles"
    identifier = "arole"
    guild_snowflake: int
    role_snowflake: int
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))


class NotAdministrator(commands.CheckFailure):
    def __init__(self, guild_snowflake: int):
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(guild_snowflake))
        message = f"You lack sufficient permissions of an administrator in the requested server ({guild.name}."
        super().__init__(message)
