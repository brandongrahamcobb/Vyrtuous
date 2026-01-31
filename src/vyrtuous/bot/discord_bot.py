"""!/bin/python3

discord_bot.py This is essentially a stripped version of Rapptz advanced_startup.py.

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
import asyncpg
import discord

from vyrtuous.inc.helpers import DISCORD_COGS


class DiscordBot(commands.Bot):

    _instance = None

    def __init__(self, *, config, db_pool: asyncpg.Pool, **kwargs):
        DiscordBot._instance = self
        intents = discord.Intents.all()
        super().__init__(
            command_prefix=config["discord_command_prefix"],
            help_command=None,
            intents=intents,
            **kwargs,
        )
        self.config = config
        self.db_pool = db_pool
        self.testing_guild_snowflake = self.config["discord_testing_guild_snowflake"]

    async def setup_hook(self):
        for cog in DISCORD_COGS:
            await self.load_extension(cog)

    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            raise RuntimeError("DiscordBot instance has not been created yet")
        return cls._instance
