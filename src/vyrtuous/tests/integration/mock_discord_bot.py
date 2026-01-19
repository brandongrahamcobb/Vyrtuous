from discord.ext import commands
import asyncpg
import discord

from vyrtuous.config import Config
from vyrtuous.inc.helpers import DISCORD_COGS


class MockBot(commands.Bot):

    def __init__(self, config: Config, db_pool: asyncpg.Pool):
        intents = discord.Intents.all()
        self.config = config
        self.db_pool = db_pool
        super().__init__(command_prefix="!", help_command=None, intents=intents)

    @classmethod
    def get_instance(cls):
        return cls()

    async def setup_hook(self):
        for cog in DISCORD_COGS:
            await self.load_extension(cog)