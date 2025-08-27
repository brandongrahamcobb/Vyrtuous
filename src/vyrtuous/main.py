''' main.py  The purpose of this program is to be the primary executable Vyrtuous.

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
import asyncio
import os
import sys
from pathlib import Path

import asyncpg
from vyrtuous.inc.helpers import *

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.config import Config
from vyrtuous.sec.discord_oauth import discord_app, DiscordOAuth, setup_discord_routes
from vyrtuous.utils.setup_logging import logger, setup_logging

PACKAGE_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGE_ROOT))

async def database_init():
    return await asyncpg.create_pool(host=os.getenv('POSTGRES_HOST'), database=os.getenv('POSTGRES_DATABASE'), user=os.getenv('POSTGRES_USER'), password=os.getenv('POSTGRES_PASSWORD'), command_timeout=30)

async def start_bot(bot, name):
    try:
        if isinstance(bot, DiscordBot):
            await bot.start(bot.api_key)
    except Exception as e:
        logger.error(f"Error running {name}: {e}")

async def main():
    config = Config().get_config()

    setup_logging(config, PATH_LOG)

    db_pool = await database_init()

    lock = asyncio.Lock()

    discord_oauth = DiscordOAuth(config)
    setup_discord_routes(discord_app, discord_oauth)

    discord_quart = asyncio.create_task(discord_app.run_task(host="0.0.0.0", port=3000))
    print("Please authenticate Discord by visiting the following URL:")
    print(discord_oauth.get_authorization_url())

    discord_bot = DiscordBot(
        config=config,
        db_pool=db_pool,
        lock=lock,
        oauth_token=discord_oauth.access_token,
    )


    tasks = [
        asyncio.create_task(start_bot(discord_bot, "DiscordBot")),
        discord_quart,
    ]

    await asyncio.gather(*tasks)

def run():
    asyncio.run(main())

if __name__ == '__main__':
    try:
        run()
    except KeyboardInterrupt:
        logger.info("Shutting down bots and server...")
