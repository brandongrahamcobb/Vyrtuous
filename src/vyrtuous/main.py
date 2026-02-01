"""main.py The purpose of this program is to be the primary executable for Vyrtuous.

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

import asyncio
import this

import debugpy

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.config import Config
from vyrtuous.db.database import Database
from vyrtuous.inc.helpers import PATH_LOG
from vyrtuous.utils.logger import logger, setup_logging


async def main():

    config = Config().get_config()
    if config["release_mode"] == "False":
        debugpy.listen(("127.0.0.1", 5678))
        debugpy.wait_for_client()

    setup_logging(config, PATH_LOG)
    db_pool = await Database(config=config).database_init()
    discord_bot = DiscordBot(config=config, db_pool=db_pool)
    await discord_bot.start(config["vyrtuous_api_key"])


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        logger.info("Shutting down bots and server...")
