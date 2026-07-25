"""!/bin/python3
guild_owner.py The purpose of this program is provide the guild owner model.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from dataclasses import dataclass

from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot


@dataclass(frozen=True)
class GuildOwner:
    guild_snowflake: int
    member_snowflake: int


class NotGuildOwner(commands.CheckFailure):
    def __init__(self, guild_snowflake: int | None):
        name = None
        bot: DiscordBot = DiscordBot.get_instance()
        if guild_snowflake:
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                raise commands.GuildNotFound(str(guild_snowflake))
            name = guild.name
        message: str = (
            f"You lack sufficient permissions of a guild owner in the requested server{f" ({name})" if name else ""}."
        )
        super().__init__(message)
