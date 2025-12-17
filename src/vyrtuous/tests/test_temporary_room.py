
from datetime import datetime, timedelta, timezone
from unittest.mock import MagicMock
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.utils.temporary_room import TemporaryRoom
import asyncio
import asyncpg
import inspect
import os
import pytest
import pytest_asyncio
import vyrtuous

MockGuild = lambda id, display_name, channels, members: type('MockGuild', (), {'id': id, 'display_name': display_name,'_channel': channels or {}, 'get_channel': lambda self, channel_id: self._channels.get(channel_id), '_members': members or {}, 'get_member': lambda self, member_id: self._members.get(member_id)})()
MockChannel = lambda id, name, guild: type('MockChannel', (), {'id': id, 'name': name, 'guild': guild})()
MockMember = lambda id, name: type('MockMember', (), {'id': id, 'name': name})()

expires_at = datetime.now(timezone.utc) + timedelta(hours=1)
member_one = MockMember(id=111111111, name="Person 1")
member_two = MockMember(id=222222222, name="Person 2")

@pytest.mark.asyncio
async def test_get_channel(bot_instance):
    assert guild.get_member(111111111) is member_one
    assert guild.get_member(999999999) is None

@pytest.mark.asyncio
async def test_get_channel(bot_instance):
    assert guild.get_channel(987654321) is temporary_channel_one
    assert temporary_channel_one.guild is guild

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

temporary_channel_one = guild.get_channel(987654321)
temporary_channel_two = guild.get_channel(978653421)
member_three = MockMember(id=333333333, name="Person 3")
member_four = MockMember(id=444444444, name="Person 4")
new_expires_at = datetime.now(timezone.utc) + timedelta(hours=2)

@pytest_asyncio.fixture(scope="session")
async def bot_instance():
    database: Optional[str] = os.getenv('POSTGRES_DB')
    host: Optional[str] = os.getenv('POSTGRES_HOST')
    password: Optional[str] = os.getenv('POSTGRES_PASSWORD')
    user: Optional[str] = os.getenv('POSTGRES_USER')
    dsn = f"postgres://{user}:{password}@{host}:{5432}/{database}"
    db_pool = await asyncpg.create_pool(dsn=dsn)
    config = {
        'discord_command_prefix': '!',
        'discord_testing_guild_snowflake': 123456789012345678
    }
    bot = DiscordBot(config=config, db_pool=db_pool)
    type(bot).guilds = property(lambda self: [guild])
    yield bot

def test_insert_into_temporary_rooms(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id, guild.id, member_two.id, temporary_channel_two.name, temporary_channel_two.id)
        rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild(guild)
        assert len(rooms) == 2
        channel_ids = {room.channel_id for room in rooms}
        assert channel_ids == {temporary_channel_one.id, temporary_channel_two.id}
        is_temp_rooms = {room.is_temp_room for room in rooms}
        assert is_temp_rooms == {True, True}
        room_owners = {room.room_owner.id for room in rooms}
        assert room_owners == {member_one.id, member_two.id}
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_temporary_room_by_guild_and_member(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id)
        rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild_and_member(guild, member_one)
        assert rooms is not None
        assert len(rooms) == 1
        for room in rooms:
            assert room is not None
            assert room.channel == temporary_channel_one
            assert room.channel_id == temporary_channel_one.id
            assert room.guild == guild
            assert room.is_temp_room == True
            assert room.room_owner == member_one
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_temporary_room_by_guild_and_room_name(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id)
        room = await TemporaryRoom.fetch_temporary_room_by_guild_and_room_name(guild, temporary_channel_one.name)
        assert room is not None
        assert room.channel == temporary_channel_one
        assert room.channel_id == temporary_channel_one.id
        assert room.guild == guild
        assert room.is_temp_room == True
        assert room.room_owner == member_one
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_all_guilds_with_temporary_rooms(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id, guild.id, member_two.id, temporary_channel_two.name, temporary_channel_two.id)
        guilds = await TemporaryRoom.fetch_all_guilds_with_temporary_rooms()
        assert guild in guilds
        for g in guilds:
            if g.id == 123456789:
                rooms = guilds[g]
                assert len(rooms) == 2
                channel_ids = {room.channel_id for room in rooms}
                assert channel_ids == {temporary_channel_one.id, temporary_channel_two.id}
                channels = {room.channel for room in rooms}
                assert channels == {temporary_channel_one, temporary_channel_two}
                is_temp_rooms = {room.is_temp_room for room in rooms}
                assert is_temp_rooms == {True, True}
                room_owners = {room.room_owner for room in rooms}
                assert room_owners == {member_one, member_two}
                room_snowflakes = {room.channel_id for room in rooms}
                assert room_snowflakes == {temporary_channel_one.id, temporary_channel_two.id}
                for room in rooms:
                    if room.channel_id == temporary_channel_one.id:
                        assert room.room_name == temporary_channel_one.name
                    elif room.channel_id == temporary_channel_two.id:
                        assert room.room_name == temporary_channel_two.name
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())
    
def test_delete_temporary_room_by_channel(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id, guild.id, member_two.id, temporary_channel_two.name, temporary_channel_two.id)
        await TemporaryRoom.delete_temporary_room_by_channel(temporary_channel_one)
        async with bot.db_pool.acquire() as conn:
            room_one = await conn.fetchrow(
                'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2',
                guild.id, temporary_channel_one.id
            )
        async with bot.db_pool.acquire() as conn:
            room_two = await conn.fetchrow(
                'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2',
                guild.id, temporary_channel_two.id
            )
        assert room_one is None
        assert room_two['room_name'] == temporary_channel_two.name
        assert room_two['owner_snowflake'] == member_two.id
        assert room_two['room_snowflake'] == temporary_channel_two.id
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_update_temporary_room_owner_snowflake(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id)
        room = TemporaryRoom(guild, temporary_channel_one.id, member_one)
        room.update_temporary_room_owner_snowflake(member_two.id)
        async with bot.db_pool.acquire() as conn:
            room_one = await conn.fetchrow(
                'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2',
                guild.id, temporary_channel_one.id
            )
        assert room_one['room_snowflake'] == temporary_channel_one.id
        assert room_one['room_name'] == temporary_channel_one.name
        assert room_one['owner_snowflake'] == member_one.id
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_update_temporary_room_name_and_room_snowflake(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id)
        room = TemporaryRoom(guild, temporary_channel_one.id, member_one)
        room.update_temporary_room_name_and_room_snowflake(temporary_channel_two.name, temporary_channel_two.id)
        async with bot.db_pool.acquire() as conn:
            room_one = await conn.fetchrow(
                'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2',
                guild.id, temporary_channel_one.id
            )
        assert room_one['room_snowflake'] == temporary_channel_one.id
        assert room_one['room_name'] == temporary_channel_one.name
        assert room_one['owner_snowflake'] == member_one.id
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_temporary_rooms_by_guild(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4), ($1, $5, $6, $7)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id, member_two.id, temporary_channel_two.name, temporary_channel_two.id)
        rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild(guild)
        assert rooms is not None
        assert len(rooms) == 2
        channel_ids = {room.channel_id for room in rooms}
        assert channel_ids == {temporary_channel_one.id, temporary_channel_two.id}
        channels = {room.channel for room in rooms}
        assert channels == {temporary_channel_one, temporary_channel_two}
        is_temp_rooms = {room.is_temp_room for room in rooms}
        assert is_temp_rooms == {True, True}
        room_owners = {room.room_owner.id for room in rooms}
        assert room_owners == {member_one.id, member_two.id}
        room_snowflakes = {room.channel_id for room in rooms}
        assert room_snowflakes == {temporary_channel_one.id, temporary_channel_two.id}
        for room in rooms:
            if room.channel_id == temporary_channel_one.id:
                assert room.room_name == temporary_channel_one.name
            elif room.channel_id == temporary_channel_two.id:
                assert room.room_name == temporary_channel_two.name
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_temporary_room_by_channel(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id, guild.id, member_two.id, temporary_channel_two.name, temporary_channel_two.id)
        room_one = await TemporaryRoom.fetch_temporary_room_by_channel(temporary_channel_one)
        assert room_one is not None
        assert room_one.channel_id == temporary_channel_one.id
        assert room_one.guild == guild
        assert room_one.is_temp_room == True
        assert room_one.room_name == temporary_channel_one.name
        assert room_one.room_owner == member_one
        room_two = await TemporaryRoom.fetch_temporary_room_by_channel(temporary_channel_two)
        assert room_two is not None
        assert room_two.channel_id == temporary_channel_two.id
        assert room_two.guild == guild
        assert room_two.is_temp_room == True
        assert room_two.room_name == temporary_channel_two.name
        assert room_two.room_owner == member_two
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())