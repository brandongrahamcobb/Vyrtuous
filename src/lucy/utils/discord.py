''' discord.py  This is essentially a stripped version of Rapptz advanced_startup.py.
    Copyright (C) 2024  github.com/brandongrahamcobb

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
from discord.ext import commands
from typing import List, Optional
from utils.setup_logging import logger
from utils.create_https_completion import Conversations
from utils.twitch import Vyrtuous

import asyncpg
import discord

conversations = Conversations()

class Lucy(commands.Bot):

    def __init__(
        self,
        *args,
        initial_extensions: List[str],
        db_pool: asyncpg.Pool,
        testing_guild_id: Optional[int] = None,
        conversations: Conversations = None,
        **kwargs
    ):
        logger.info('Initializing Lucy bot...')
        try:
            super().__init__(*args, **kwargs)
            self.db_pool = db_pool
            self.testing_guild_id = testing_guild_id
            self.initial_extensions = initial_extensions
            self.conversations = conversations
            logger.info('Lucy bot initialized successfully.')
        except Exception as e:
            logger.error(f'Error during initialization: {e}')

    async def setup_hook(self) -> None:
        try:
            logger.info('Starting setup_hook...')
            for cog in self.initial_extensions:
                logger.debug(f'Loading extension: {cog}')
                await self.load_extension(cog)
                logger.info(f'Extension {cog} loaded successfully.')

            if self.testing_guild_id:
                logger.debug(f'Setting up testing guild with ID: {self.testing_guild_id}')
                guild = discord.Object(self.testing_guild_id)
                self.tree.copy_global_to(guild=guild)
                await self.tree.sync(guild=guild)
                logger.info(f'Commands synced with testing guild ID: {self.testing_guild_id}')
            else:
                logger.info('No testing guild ID provided; skipping guild sync.')
        except Exception as e:
            logger.error(f'Error during setup_hook: {e}')
