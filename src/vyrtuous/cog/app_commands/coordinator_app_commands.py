"""!/bin/python3

coordinator_app_commands.py A discord.py cog containing coordinator commands for the Vyrtuous bot.

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

from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.coordinator.coordinator import Coordinator


class CoordinatorAppCommands(commands.Cog):
    ROLE = Coordinator

    def __init__(self, *, bot: DiscordBot | None = None):
        self.__bot = bot


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorAppCommands(bot))
