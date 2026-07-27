"""!/bin/python3
guild_owner.py The purpose of this program is to extend DatabaseFactory to provide the guild owner class.

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

from dataclasses import dataclass

from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot


@dataclass(frozen=True)
class GuildOwner:
    guild_snowflake: int
    member_snowflake: int


class NotGuildOwner(app_commands.CheckFailure, commands.CheckFailure):
    def __init__(self, guild_snowflake: int | None):
        name = None
        bot: DiscordBot = DiscordBot.get_instance()
        if guild_snowflake:
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                raise GuildNotFound(str(guild_snowflake))
            name = guild.name
        message: str = (
            f"You lack sufficient permissions of a guild owner in the requested server{f" ({name})" if name else ""}."
        )
        super().__init__(message)
