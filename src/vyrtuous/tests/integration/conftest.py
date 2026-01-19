"""conftest.py The purpose of this program is to support glass box, integration and unit testing for Vyrtuous.

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
import asyncpg
import pytest
import pytest_asyncio

from vyrtuous.config import Config
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.tests.integration.mock_discord_bot import MockBot
from vyrtuous.tests.integration.mock_database import dsn

@pytest.fixture
def cf():
    cf = Config().get_config()
    yield cf

@pytest_asyncio.fixture
async def db():
    db_pool = await asyncpg.create_pool(dsn=dsn)
    yield db_pool
    await db_pool.close()

@pytest_asyncio.fixture
async def bot(cf, db):
    bot = MockBot(config=cf, db_pool=db)
    DiscordBot._instance = bot
    await bot.setup_hook()
    return bot