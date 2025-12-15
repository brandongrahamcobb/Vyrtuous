''' discord.py This is essentially a stripped version of Rapptz advanced_startup.py.

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
'''
from collections import defaultdict
from discord.ext import commands
from vyrtuous.inc.helpers import DISCORD_COGS
from vyrtuous.utils.setup_logging import logger

import asyncio
import asyncpg
import discord

class DiscordBot(commands.Bot):
    
    _instance = None
        
    def __init__(self, *, config, db_pool: asyncpg.Pool, **kwargs): # oauth_token,
        try:
            DiscordBot._instance = self
            self.lock = asyncio.Lock()
            intents = discord.Intents.all()
            super().__init__(command_prefix=config['discord_command_prefix'],
                             help_command=None,
                             intents=intents,
                             **kwargs)
            self.config = config
            self.db_pool = db_pool
            self.testing_guild_id = self.config['discord_testing_guild_id']
        except Exception as e:
            logger.error(f'Error during Discord bot initialization: {e}')

    async def setup_hook(self) -> None:
        for cog in DISCORD_COGS:
            await self.load_extension(cog)

    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            raise RuntimeError("DiscordBot instance has not been created yet")
        return cls._instance
