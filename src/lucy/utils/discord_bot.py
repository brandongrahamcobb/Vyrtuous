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
import asyncpg
import discord
import logging

from discord.ext import commands
from typing import List, Optional
from .setup_logging import logger
from .create_https_completion import Conversations


class DiscordBot(commands.Bot):
    def __init__(self, *, config, db_pool: asyncpg.Pool, conversations, lock, oauth_token,**kwargs):
        logger.info('Initializing Lucy Discord bot...')
        try:
            intents = discord.Intents.all()
            super().__init__(command_prefix=config['discord_command_prefix'], intents=intents, **kwargs)

            self.config = config
            self.conversations = conversations
            self.db_pool = db_pool
            self.lock = lock
            self.oauth_token = oauth_token
            self.api_key = self.config['api_keys']['Discord']['api_key']
            self.testing_guild_id = self.config['discord_testing_guild_id']

            logger.info('Lucy Discord bot initialized successfully.')
        except Exception as e:
            logger.error(f'Error during Discord bot initialization: {e}')

    async def setup_hook(self) -> None:
        try:
            logger.info('Starting Discord bot setup_hook...')
            cogs = [
                'lucy.cogs.hybrid',
                'lucy.cogs.indica',
                'lucy.cogs.sativa',
            ]
            for cog in cogs:
                logger.debug(f'Loading extension: {cog}')
                await self.load_extension(cog)
                logger.info(f'Extension {cog} loaded successfully.')
            if self.testing_guild_id:
                logger.debug(f'Setting up testing guild with ID: {self.testing_guild_id}')
                guild = discord.Object(id=self.testing_guild_id)
                self.tree.copy_global_to(guild=guild)
                await self.tree.sync(guild=guild)
                logger.info(f'Commands synced with testing guild ID: {self.testing_guild_id}')
            else:
                logger.info('No testing guild ID provided; skipping guild sync.')
        except Exception as e:
            logger.error(f'Error during Discord bot setup_hook: {e}')
