"""mock_discord_bot.py The purpose of this program is to support integration testing for Vyrtuous.

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

from unittest.mock import AsyncMock, Mock

from discord.ext import commands
import asyncpg
import discord

from vyrtuous.config import Config
from vyrtuous.inc.helpers import DISCORD_COGS


class MockBot(commands.Bot):

    def __init__(self, config: Config, db_pool: asyncpg.Pool):
        intents = discord.Intents.all()
        self.config = config
        self.db_pool = db_pool
        self._guilds = []
        self._tree = AsyncMock()
        self._tree.sync = AsyncMock()
        self._tree.add_command = Mock()
        self._tree.remove_command = AsyncMock()
        self._tree.copy_global_to = AsyncMock()
        super().__init__(command_prefix="!", help_command=None, intents=intents)

    @classmethod
    def get_instance(cls):
        return cls()

    def get_guild(self, target: int):
        for guild in self._guilds:
            if target == guild.id:
                return guild

    async def setup_hook(self):
        for cog in DISCORD_COGS:
            if cog != "vyrtuous.cogs.scheduled_tasks":
                await self.load_extension(cog)

    @property
    def tree(self):
        return self._tree
