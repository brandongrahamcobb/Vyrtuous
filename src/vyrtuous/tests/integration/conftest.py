
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