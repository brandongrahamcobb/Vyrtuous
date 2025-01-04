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

from lucy.utils.config import Config
from lucy.utils.create_https_completion import Conversations
from lucy.utils.discord_utils import Lucy
from lucy.utils.helpers import *
from lucy.utils.increment_version import increment_version
from lucy.utils.linkedin import LinkedIn
from lucy.utils.setup_logging import setup_logging
from lucy.utils.twitch import app, ensure_token, Vyrtuous
from pathlib import Path

import aiohttp
import asyncio
import asyncpg
import discord
import os
import sys
import yaml

package_root = Path(__file__).resolve().parent
sys.path.insert(0, str(package_root))

config = Config().get_config()
conversations = Conversations()

start_discord_event = asyncio.Event()
start_linkedin_event = asyncio.Event()
start_twitch_event = asyncio.Event()

async def discord_init(db_pool, lock):
    discord_bot = Lucy(
        config=config,
        command_prefix=config['discord_command_prefix'],
        db_pool=db_pool,
        initial_extensions=DISCORD_COGS,
        intents=DISCORD_INTENTS,
        testing_guild_id=config['discord_testing_guild_id'],
        conversations=conversations,
        lock=lock
    )
    discord_bot.config = config
    discord_bot.db_pool = db_pool
    return discord_bot

async def database_init():
    db_pool = await asyncpg.create_pool(
        database='lucy',
        user='postgres',
        command_timeout=30
    )
    return db_pool

async def linkedin_init(db_pool):
    linkedin_bot = LinkedIn(config, db_pool)
    linkedin_bot.config = config
    linkedin_bot.db_pool = db_pool
    return linkedin_bot

async def twitch_init(db_pool, lock):
    token = await ensure_token()
    if not token:
        logger.warning('No valid token available on startup. Please authorize the application.')
        return
    twitch_bot = Vyrtuous(
        config=config,
        db_pool=db_pool,
        conversations=conversations,
        token=token,
        lock=lock
    )
    twitch_bot.config = config
    twitch_bot.db_pool = db_pool
    return twitch_bot

def trigger_discord_start():
    logger.info("Triggering Discord bot start...")
    start_discord_event.set()  # Signal to start Discord bot

def trigger_linkedin_start():
    logger.info("Triggering LinkedIn bot start...")
    start_linkedin_event.set()  # Signal to start LinkedIn bot

def trigger_twitch_start():
    logger.info("Triggering Twitch bot start...")
    start_twitch_event.set()  # Signal to start Twitch bot

async def main():
    increment_version()
    setup_logging(config, PATH_LOG)

    lock = asyncio.Lock()
    db_pool = await database_init()
    app_task = asyncio.create_task(app.run_task(host="0.0.0.0", port=5000))

    print('(You have 30 seconds to authenticate')
    await asyncio.sleep(30)

    discord_bot = await discord_init(db_pool, lock)
    linkedin_bot = await linkedin_init(db_pool)
    twitch_bot = await twitch_init(db_pool, lock)


    discord_task = asyncio.create_task(discord_bot.start(config['discord_token']))
    linkedin_task = asyncio.create_task(linkedin_bot.main())
    twitch_task = asyncio.create_task(twitch_bot.start())
    await asyncio.gather(app_task, discord_task, twitch_task)

def run():
    import asyncio

    asyncio.run(main())

if __name__ == '__main__':
    import threading

    discord_trigger_thread = threading.Thread(target=trigger_discord_start)
    linkedin_trigger_thread = threading.Thread(target=trigger_linkedin_start)
    twitch_trigger_thread = threading.Thread(target=trigger_twitch_start)

    discord_trigger_thread.start()
    linkedin_trigger_thread.start()
    twitch_trigger_thread.start()

    run()

