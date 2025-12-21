
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
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bot.discord_client import DiscordClient
from vyrtuous.config import Config
from vyrtuous.inc.helpers import *
import asyncio
import asyncpg
import discord
import os
import pytest
import pytest_asyncio

config = Config().get_config()


def make_member(id, name, bot=True, voice_channel=False):
    return type(
        'MockMember',
        (),
        {
            'bot': bot,
            'id': id,
            'name': name,
            'mention': f'<@{id}>',
            'voice': SimpleNamespace(channel=voice_channel, mute=False)
        }
    )()
class MockMessage:
    def __init__(self, *, content, channel, guild, id, author):
        self.content = content
        self.channel = channel
        self.guild = guild
        self.id = id
        self.author = author
        self.bot = False
        self.attachments = []
        self.embeds = []
        self._state = None

expires_at = datetime.now(timezone.utc) + timedelta(hours=1)
self_member_obj = make_member(id=config['discord_testing_self_member_snowflake'], name="Person 1")
self_member_obj.bot = True
dummy_member_obj = make_member(id=config['discord_testing_dummy_member_snowflake'], name="Person 2")

def make_guild(id, display_name, members, channel_defs, owner_id, roles):
    guild = type(
        'MockGuild',
        (),
        {
            'id': id,
            'display_name': display_name,
            '_channels': {},
            'get_channel': lambda self, channel_id: self._channels.get(channel_id),
            '_members': members,
            'get_member': lambda self, member_id: self._members.get(member_id),
            'owner_id': owner_id,
            "get_role": lambda self, role_id: roles.get(role_id),
        }
    )()
    channels = {}
    for cid, (name, channel_type) in channel_defs.items():
        channels[cid] = make_mock_channel(cid, name, guild, channel_type)
    guild._channels = channels
    for member in members.values():
        client_channel = list(channels.values())[0]
        member.voice.channel = client_channel
        client_channel.members.append(member)
    return guild

def make_mock_channel(id, name, guild, channel_type=None):
    channel = type(
        'MockChannel',
        (),
        {
            'guild': guild,
            'id': id,
            'members': [],
            'mention': f'<@{id}>',
            'messages': [],
            'name': name,
            'type': channel_type
        }
    )()
    return channel

members = {
    self_member_obj.id: self_member_obj,
    dummy_member_obj.id: dummy_member_obj
}

channel_defs = {
    config['discord_testing_first_channel_snowflake']: ("Test Voice Channel", discord.ChannelType.voice),
    config['discord_testing_second_channel_snowflake']: ("Test Stage Channel", discord.ChannelType.stage_voice),
    config['discord_testing_text_channel_snowflake']: ("Test Text Channel", discord.ChannelType.text)
}

mock_role = SimpleNamespace(id=987654321, mention="<@&987654321>")
roles = {mock_role.id: mock_role}

guild_obj = make_guild(
    id=config['discord_testing_guild_snowflake'],
    display_name="Test Guild",
    channel_defs=channel_defs,
    members=members,
    owner_id=self_member_obj.id,
    roles=roles
)

bot_channel_obj = make_mock_channel(
    id=config['discord_testing_first_channel_snowflake'],
    name="Voice Channel 1",
    guild=guild_obj,
    channel_type=discord.ChannelType.voice
)
client_channel_obj = make_mock_channel(
    id=config['discord_testing_second_channel_snowflake'],
    name="Voice Channel 2",
    guild=guild_obj,
    channel_type=discord.ChannelType.voice
)
text_channel_obj = make_mock_channel(
    id=config['discord_testing_text_channel_snowflake'],
    name="text-channel",
    guild=guild_obj,
    channel_type=discord.ChannelType.text
)
new_expires_at = datetime.now(timezone.utc) + timedelta(hours=2)

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
def bot_channel():
    return bot_channel_obj

@pytest.fixture(scope="function")
def client_channel():
    return client_channel_obj

@pytest.fixture(scope="function")
def text_channel():
    return text_channel_obj

@pytest.fixture(scope="function")
def self_member():
    return self_member_obj

@pytest.fixture(scope="function")
def dummy_member():
    return dummy_member_obj


@pytest.fixture(scope="function")
def guild():
    return guild_obj

@pytest.fixture(scope="function")
def _state():
    return None

@pytest.fixture(scope="function")
def config(bot):
    config = bot.config
    yield config

@pytest.fixture(scope="function")
def command():
    command: str = 'cap'
    yield command

@pytest.fixture(scope="function")
def prefix(config):
    prefix = config['discord_command_prefix']
    yield prefix

async def admin_cleanup(guild_id: Optional[int], self_member_id: Optional[int]):
    bot = DiscordBot.get_instance()
    async with bot.db_pool.acquire() as conn:
        await conn.execute('''
            UPDATE users SET administrator_guild_ids=$2 WHERE discord_snowflake=$1
        ''', int(self_member_id), [int(guild_id)])

async def admin_initiation(guild_id: Optional[int], self_member_id: Optional[int]):
    bot = DiscordBot.get_instance()
    async with bot.db_pool.acquire() as conn:
        await conn.execute('''
            INSERT INTO users (discord_snowflake, developer_guild_ids, updated_at, created_at)
            VALUES ($1, $2, NOW(), NOW())
            ON CONFLICT (discord_snowflake) 
            DO UPDATE SET developer_guild_ids = $2, updated_at = NOW()
        ''', int(self_member_id), [int(guild_id)])

