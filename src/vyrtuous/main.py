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
import os
import asyncio

import asyncpg
import debugpy

from vyrtuous.inc.helpers import PATH_LOG
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.config import Config
from vyrtuous.utils.setup_logging import logger, setup_logging

async def database_init():
    return await asyncpg.create_pool(
            host=os.getenv('POSTGRES_HOST'),
            database=os.getenv('POSTGRES_DB'),
            user=os.getenv('POSTGRES_USER'),
            password=os.getenv('POSTGRES_PASSWORD'),
            command_timeout=30)

async def main():
#    debugpy.listen(("0.0.0.0", 5678))
#    debugpy.wait_for_client() 
    config = Config().get_config()
    setup_logging(config, PATH_LOG)
    db_pool = await database_init()

    discord_bot = DiscordBot(
        config=config,
        db_pool=db_pool,
    )
    
    await discord_bot.start(config['discord_api_key'])

if __name__ == '__main__':
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        logger.info("Shutting down bots and server...")
