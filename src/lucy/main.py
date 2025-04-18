''' main.py  The purpose of this program is to be the primary executable for Vyrtuous.
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
from lucy.bots.discord_bot import DiscordBot
from lucy.bots.linkedin_bot import LinkedInBot
from lucy.bots.twitch_bot import TwitchBot
from lucy.config import Config
from lucy.utils.handlers.ai_manager import Completions
from lucy.utils.inc.helpers import *
from lucy.utils.inc.increment_version import increment_version
from lucy.utils.inc.setup_logging import setup_logging
from lucy.utils.sec.discord_oauth import discord_app, DiscordOAuth, setup_discord_routes
from lucy.utils.sec.linkedin_oauth import linkedin_app, LinkedInOAuth, setup_linkedin_routes
from lucy.utils.sec.patreon_oauth import patreon_app, PatreonOAuth, setup_patreon_routes
from lucy.utils.sec.twitch_oauth import twitch_app, TwitchOAuth, setup_twitch_routes
from pathlib import Path

import asyncio
import asyncpg
import logging
import sys

PACKAGE_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGE_ROOT))

async def database_init():
    return await asyncpg.create_pool(database='lucy', user='postgres', command_timeout=30)

async def start_bot(bot, name):
    try:
        if isinstance(bot, DiscordBot):
            await bot.start(bot.api_key)
        elif isinstance(bot, LinkedInBot):
            await bot.run()
        elif isinstance(bot, TwitchBot):
            await bot.start()
    except Exception as e:
        logger.error(f"Error running {name}: {e}")

async def main():
    config = Config().get_config()
    increment_version()
    setup_logging(config, PATH_LOG)

    db_pool = await database_init()

    completions = Completions()
    lock = asyncio.Lock()

    discord_oauth = DiscordOAuth(config)
    setup_discord_routes(discord_app, discord_oauth)
    linkedin_oauth = LinkedInOAuth(config)
    setup_linkedin_routes(linkedin_app, linkedin_oauth)
    patreon_oauth = PatreonOAuth(config)
    setup_patreon_routes(patreon_app, patreon_oauth)
    twitch_oauth = TwitchOAuth(config)
    setup_twitch_routes(twitch_app, twitch_oauth)

    discord_quart = asyncio.create_task(discord_app.run_task(host="0.0.0.0", port=5000))
    linkedin_quart = asyncio.create_task(linkedin_app.run_task(host="0.0.0.0", port=5001))
    linkedin_quart = asyncio.create_task(patreon_app.run_task(host="0.0.0.0", port=5003))
    twitch_quart = asyncio.create_task(twitch_app.run_task(host="0.0.0.0", port=5002))

    print("Please authenticate Discord by visiting the following URL:")
    print(discord_oauth.get_authorization_url())
    print("Please authenticate LinkedIn by visiting the following URL:")
    print(linkedin_oauth.get_authorization_url())
    print("Please authenticate Patreon by visiting the following URL:")
    print(patreon_oauth.get_authorization_url())
    print("Please authenticate Twitch by visiting the following URL:")
    print(twitch_oauth.get_authorization_url())
    await asyncio.sleep(8)

    discord_bot = DiscordBot(
        config=config,
        db_pool=db_pool,
        completions=completions,
        lock=lock,
        oauth_token=discord_oauth.access_token,
    )
    linkedin_bot = LinkedInBot(
        config=config,
        db_pool=db_pool,
        completions=completions,
        lock=lock,
        oauth_token=linkedin_oauth.access_token,
    )
    twitch_bot = TwitchBot(
        config=config,
        db_pool=db_pool,
        completions=completions,
        lock=lock,
        oauth_token=twitch_oauth.access_token,
    )

    tasks = [
        asyncio.create_task(start_bot(discord_bot, "DiscordBot")),
        asyncio.create_task(start_bot(linkedin_bot, "LinkedInBot")),
        asyncio.create_task(start_bot(twitch_bot, "TwitchBot")),
        discord_quart,
        linkedin_quart,
        twitch_quart,
    ]

    await asyncio.gather(*tasks)

def run():
    asyncio.run(main())

if __name__ == '__main__':
    try:
        run()
    except KeyboardInterrupt:
        logger.info("Shutting down bots and server...")
