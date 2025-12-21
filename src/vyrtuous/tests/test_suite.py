
''' test_suite.py The purpose of this program is to provide the shared test variables for tests using Discord objects.
    Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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
from types import SimpleNamespace
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bot.discord_client import DiscordClient
from vyrtuous.config import Config
from vyrtuous.inc.helpers import *
from vyrtuous.tests.make_mock_objects import *
import asyncpg
import discord
import os
import pytest
import pytest_asyncio

privileged_author_obj = make_mock_member(
    bot=True,
    id=PRIVILEGED_AUTHOR_ID,
    name=PRIVILEGED_AUTHOR_NAME
)

not_priveleged_author_obj = make_mock_member(
    bot=False,
    id=NOT_PRIVILEGED_AUTHOR_ID,
    name=NOT_PRIVILEGED_AUTHOR_NAME
)

guild_obj = make_mock_guild(
    channel_defs={
        VOICE_CHANNEL_ONE_ID: (
            VOICE_CHANNEL_ONE_NAME,
            discord.ChannelType.voice
        ),
        VOICE_CHANNEL_TWO_ID: (
            VOICE_CHANNEL_ONE_NAME,
            discord.ChannelType.voice
        ),
        TEXT_CHANNEL_ID: (
            TEXT_CHANNEL_NAME,
            discord.ChannelType.text
        )
    },
    id=GUILD_ID,
    members={
        PRIVILEGED_AUTHOR_ID: privileged_author_obj,
        NOT_PRIVILEGED_AUTHOR_ID: not_priveleged_author_obj
    },
    name=GUILD_NAME,
    owner_id=PRIVILEGED_AUTHOR_ID,
    roles={
        ROLE_ID: SimpleNamespace(
            id=ROLE_ID,
            mention="<@&{ROLE_ID}>"
        )
    }
)

voice_channel_one_obj = make_mock_channel(
    channel_type=discord.ChannelType.voice,
    guild=guild_obj,
    id=VOICE_CHANNEL_ONE_ID,
    name=VOICE_CHANNEL_ONE_NAME
)

voice_channel_two_obj = make_mock_channel(
    channel_type=discord.ChannelType.voice,
    guild=guild_obj,
    id=VOICE_CHANNEL_TWO_ID,
    name=VOICE_CHANNEL_TWO_NAME
)

text_channel_obj = make_mock_channel(
    channel_type=discord.ChannelType.text,
    guild=guild_obj,
    id=TEXT_CHANNEL_ID,
    name=TEXT_CHANNEL_NAME
)

@pytest_asyncio.fixture(scope="function")
async def bot():
    database: Optional[str] = os.getenv('POSTGRES_DB')
    host: Optional[str] = os.getenv('POSTGRES_HOST')
    password: Optional[str] = os.getenv('POSTGRES_PASSWORD')
    user: Optional[str] = os.getenv('POSTGRES_USER')
    dsn = f"postgres://{user}:{password}@{host}:{5432}/{database}"
    if not all([database, host, password, user]):
        pytest.skip("Database environment variables not set")
    db_pool = await asyncpg.create_pool(dsn=dsn)
    config = Config().get_config()
    bot = DiscordBot(config=config, db_pool=db_pool)
    for cog in ('vyrtuous.cogs.admin_commands', 'vyrtuous.cogs.event_listeners', 'vyrtuous.cogs.aliases'):
        await bot.load_extension(cog)
    type(bot).guilds = property(lambda self: [guild_obj])
    bot._state = make_mock_state()
    yield bot
    await db_pool.close()

@pytest_asyncio.fixture(scope="function")
async def client():
    database: Optional[str] = os.getenv('POSTGRES_DB')
    host: Optional[str] = os.getenv('POSTGRES_HOST')
    password: Optional[str] = os.getenv('POSTGRES_PASSWORD')
    user: Optional[str] = os.getenv('POSTGRES_USER')
    dsn = f"postgres://{user}:{password}@{host}:{5432}/{database}"
    if not all([database, host, password, user]):
        pytest.skip("Database environment variables not set")
    db_pool = await asyncpg.create_pool(dsn=dsn)
    config = Config().get_config()
    client = DiscordClient(config=config, db_pool=db_pool)
    type(bot).guilds = property(lambda self: [guild_obj])
    yield client
    await db_pool.close()

@pytest.fixture(scope="function")
def voice_channel_one():
    return voice_channel_one_obj

@pytest.fixture(scope="function")
def voice_channel_two():
    return voice_channel_two_obj

@pytest.fixture(scope="function")
def text_channel():
    return text_channel_obj

@pytest.fixture(scope="function")
def privileged_author():
    return privileged_author_obj

@pytest.fixture(scope="function")
def not_privileged_author():
    return not_priveleged_author_obj


@pytest.fixture(scope="function")
def guild():
    return guild_obj

@pytest.fixture(scope="function")
def config(bot):
    config = bot.config
    yield config

@pytest.fixture(scope="function")
def prefix(config):
    prefix = config['discord_command_prefix']
    yield prefix

async def capturing_send(self, ctx, content=None, embed=None, allowed_mentions=None, **kwargs):
    voice_channel_one.messages.append(content)
    return make_mock_message(
        allowed_mentions=allowed_mentions,
        author=privileged_author_obj,
        channel=voice_channel_one,
        content=content,
        embeds=[embed],
        guild=ctx.guild,
        id=MESSAGE_ID
    )

async def edit(self, **kwargs):
    for k, v in kwargs.items():
        setattr(self, k, v)
    return self