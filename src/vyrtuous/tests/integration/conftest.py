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

import asyncio
import asyncpg
import pytest
import pytest_asyncio
from discord.ext.commands import Context, view as cmd_view
from discord import Interaction

from vyrtuous.config import Config
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.tests.integration.mock_discord_bot import MockBot
from vyrtuous.tests.integration.mock_database import dsn


NOT_PRIVILEGED_AUTHOR_SNOWFLAKE_ONE = 10000000000000002


@pytest.fixture
def cf(monkeypatch):
    monkeypatch.setenv("DISCORD_OWNER_ID", str(NOT_PRIVILEGED_AUTHOR_SNOWFLAKE_ONE))
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


def context(bot, message, prefix):
    view = cmd_view.StringView(message.content)
    view.skip_string(prefix)
    command_name = view.get_word()
    ctx = Context(bot=bot, message=message, prefix=prefix, view=view)
    ctx.invoked_with = command_name
    ctx.command = bot.get_command(command_name)
    view.skip_ws()
    return ctx


def interaction(context):
    interaction = Interaction(context=context)
    return interaction
