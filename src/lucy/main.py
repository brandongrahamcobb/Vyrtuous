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
#start_linkedin_event = asyncio.Event()
start_twitch = asyncio.Event()

async def discord_init():
    discord_bot = Lucy(
        command_prefix=config['discord_command_prefix'],
#        db_pool=db_pool,
        initial_extensions=config['discord_cogs'],
        intents=discord.Intents.all(),
        testing_guild_id=config['discord_testing_guild_id'],
    )
    discord_bot.config = config
    await discord_bot.start(config['discord_token'])

async def twitch_init():
    token = await startup()
    if not token:
        logger.warning('No valid token available on startup. Please authorize the application.')
        return
    twitch_bot = Vyrtuous(bot, token)
    await twitch_bot.start()

def trigger_discord_start():
    print("Triggering Discord bot start...")
    start_twitch_event.set()  # Signal to start Discord bot

def trigger_twitch_start():
    print("Triggering Twitch bot start...")
    start_twitch_event.set()  # Signal to start Twitch bot

#def trigger_linkedin_start():
 #   print("Triggering LinkedIn bot start...")
  #  start_linkedin_event.set()  # Signal to start LinkedIn bot

#async def linkedin_init():
 #   linkedin_bot = LinkedInBot(config)
  #  await linkedin_bot.start()

async def main():
    setup_logging(config, PATH_LOG)
#    async with asyncpg.create_pool(
 #       database='lucy',
  #      user='postgres',
   #     command_timeout=30
    #) as db_pool:
    app_task = asyncio.create_task(app.run_task(host="0.0.0.0", port=5000))
    discord_task = asyncio.create_task(discord_init())
 #   linkedin_task = asyncio.create_task(linkedin_init()),
    twitch_task = asyncio.create_task(twitch_init())
    await asyncio.gather(app_task, discord_task, twitch_task)

def run():
    import asyncio
    asyncio.run(main())

if __name__ == '__main__':
    import threading
    twitch_trigger_thread = threading.Thread(target=trigger_twitch_start)
    linkedin_trigger_thread = threading.Thread(target=trigger_linkedin_start)

    twitch_trigger_thread.start()
    linkedin_trigger_thread.start()

    run()

