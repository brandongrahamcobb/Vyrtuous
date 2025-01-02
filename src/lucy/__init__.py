''' __init__.py  The purpose of this program is to provide the advanced_startup.py logic provided by Rapptz; from cd ../.
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
from datetime import datetime, timedelta
from quart import Quart, Response, redirect, request, session
from lucy.utils.config import Config
from lucy.utils.create_https_completion import Conversations
from lucy.utils.discord_utils import Lucy
from lucy.utils.increment_version import increment_version
from lucy.utils.setup_logging import setup_logging
from lucy.utils.twitch import app, startup
from lucy.utils.helpers import *

import aiohttp
import asyncio
import asyncpg
import discord
import os
import yaml

async def main():
    config = Config().get_config()
    conversations = Conversations()
    async with asyncpg.create_pool(
        database='lucy',
        user='postgres',
        command_timeout=30
    ) as pool:
        setup_logging(config, PATH_LOG)
        tasks = []
        async with Lucy(
            command_prefix=config['discord_command_prefix'],
            db_pool=pool,
            initial_extensions=['lucy.cogs.hybrid', 'lucy.cogs.indica', 'lucy.cogs.sativa'],
            intents=discord.Intents.all(),
            testing_guild_id=config['discord_testing_guild_id'],
            conversations=conversations
        ) as bot:
            bot.config = config
            increment_version(config, PATH_CONFIG_YAML)
            tasks.append(asyncio.create_task(startup()))
            tasks.append(asyncio.create_task(app.run_task(port=5000)))
            tasks.append(asyncio.create_task(bot.start(config['discord_token'])))
            await asyncio.gather(*tasks)

def run():
    import asyncio
    asyncio.run(main())
