
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
from vyrtuous.inc.helpers import *
from vyrtuous.bot.discord_bot import DiscordBot
import asyncpg
import os
import pytest
import pytest_asyncio

MockGuild = lambda id, display_name, channels, members: type('MockGuild', (), {'id': id, 'display_name': display_name,'_channel': channels or {}, 'get_channel': lambda self, channel_id: self._channels.get(channel_id), '_members': members or {}, 'get_member': lambda self, member_id: self._members.get(member_id)})()
MockChannel = lambda id, name, guild: type('MockChannel', (), {'id': id, 'name': name, 'guild': guild})()
MockMember = lambda id, name: type('MockMember', (), {'id': id, 'name': name})()

expires_at = datetime.now(timezone.utc) + timedelta(hours=1)
member_one = MockMember(id=111111111, name="Person 1")
member_two = MockMember(id=222222222, name="Person 2")
member_three = MockMember(id=333333333, name="Person 3")
member_four = MockMember(id=444444444, name="Person 4")

def make_guild(id, display_name, members, channel_defs):
    guild = type(
        'MockGuild', (), {
            'id': id,
            'display_name': display_name,
            '_members': members or {},
            '_channels': {},  # placeholder, will fill below
            'get_member': lambda self, member_id: self._members.get(member_id),
            'get_channel': lambda self, channel_id: self._channels.get(channel_id)
        }
    )()
    channels = {
        cid: MockChannel(cid, name, guild)
        for cid, name in channel_defs.items()
    }
    guild._channels = channels
    return guild
    
guild = make_guild(
    id=123456789,
    display_name="Test Guild",
    channel_defs={
        987654321: "Test Stage Channel",
        978653421: "Test Stage Channel 2"
    },
    members={
        111111111: member_one,
        222222222: member_two
    }
)

channel_one = guild.get_channel(987654321)
channel_two = guild.get_channel(978653421)
new_expires_at = datetime.now(timezone.utc) + timedelta(hours=2)

@pytest_asyncio.fixture(scope="session")
async def bot_instance():
    database: Optional[str] = os.getenv('POSTGRES_DB')
    host: Optional[str] = os.getenv('POSTGRES_HOST')
    password: Optional[str] = os.getenv('POSTGRES_PASSWORD')
    user: Optional[str] = os.getenv('POSTGRES_USER')
    dsn = f"postgres://{user}:{password}@{host}:{5432}/{database}"
    if not all([database, host, password, user]):
        pytest.skip("Database environment variables not set")
    db_pool = await asyncpg.create_pool(dsn=dsn)
    config = {
        'discord_command_prefix': '!',
        'discord_testing_guild_snowflake': 123456789012345678
    }
    bot = DiscordBot(config=config, db_pool=db_pool)
    type(bot).guilds = property(lambda self: [guild])
    yield bot

@pytest.mark.asyncio
async def test_get_channel(bot_instance):
    assert guild.get_member(111111111) is member_one
    assert guild.get_member(999999999) is None

@pytest.mark.asyncio
async def test_get_channel(bot_instance):
    assert guild.get_channel(987654321) is channel_one
    assert channel_one.guild is guild