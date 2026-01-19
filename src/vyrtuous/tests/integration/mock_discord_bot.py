from discord.ext import commands
import asyncpg
import discord
import pytest
import pytest_asyncio

from vyrtuous.tests.integration.mock_database import dsn
from vyrtuous.inc.helpers import DISCORD_COGS

@pytest_asyncio.fixture
async def bot():
    bot = MockBot.get_instance()
    await bot.setup_hook()
    print("test")
    return bot

@pytest_asyncio.fixture
async def db():
    db_pool = await asyncpg.create_pool(dsn=dsn)
    yield db_pool
    db_pool.close()

class MockBot(commands.Bot):

    def __init__(self):
        intents = discord.Intents.all()
        super().__init__(command_prefix="!", intents=intents)

    @classmethod
    def get_instance(cls):
        return cls()

    async def setup_hook(self):
        for cog in DISCORD_COGS:
            await self.load_extension(cog)