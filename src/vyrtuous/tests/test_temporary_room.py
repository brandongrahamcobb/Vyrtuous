from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.utils.temporary_room import TemporaryRoom
import asyncio
import asyncpg
from datetime import datetime, timedelta, timezone
import inspect
import os
import pytest
import pytest_asyncio
import vyrtuous

MockGuild = lambda id, display_name: type('MockGuild', (), {'id': id, 'display_name': display_name})()
MockChannel = lambda id, name, guild: type('MockChannel', (), {'id': id, 'name': name, 'guild': guild})()
MockMember = lambda id, name: type('MockMember', (), {'id': id, 'name': name})()

expires_at = datetime.now(timezone.utc) + timedelta(hours=1)
guild = MockGuild(id=123456789, display_name="Test Guild")
temorary_channel_one = MockChannel(id=987654321, name="Test Stage Channel", guild=guild)
temorary_channel_two = MockChannel(id=978653421, name="Test Stage Channel 2", guild=guild)
member_one = MockMember(id=111111111, name="Person 1")
member_two = MockMember(id=222222222, name="Person 2")
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
        'discord_testing_guild_id': 123456789012345678
    }
    bot = DiscordBot(config=config, db_pool=db_pool)
    yield bot
    
def test_insert_into_temporary_rooms(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_id = 123456789')
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$3, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id, guild.id, member_two.id, temporary_channel_two.name, temporary_channel_two.id)
        rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild(guild)
        assert len(rooms) == 2
        room_owners = {room.room_owner.id for room in rooms}
        assert room_owners == {member_one.id, member_two.id}
        channel_ids = {room.channel_id for room in rooms}
        assert channel_ids == {temporary_channel_one.id, temporary_channel_two.id}
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def fetch_all_guilds_with_temporary_rooms(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_id = 123456789')
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$3, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id, guild.id, member_two.id, temporary_channel_two.name, temporary_channel_two.id)
        guilds = await TemporaryRoom.fetch_all_guilds_with_temporary_rooms()
        assert guild.id in guilds
        for guild in guilds:
            if guild.id == 123456789:
                rooms = guilds[guild]
                assert len(rooms) == 2
                channel_ids = {room.channel_id for room in rooms}
                assert channel_ids == {temporary_channel_one.id, temporary_channel_two.id}
                is_temp_rooms = {room.is_temp_room for room in rooms}
                assert is_temp_rooms == {True, True}
                room_owners = {room.room_owner.id for room in rooms}
                assert room_owners == {member_one.id, member_two.id}
                room_snowflakes = {room.room_snowflake for room in rooms}
                assert room_snowflakes == {temporary_channel_one.id, temporary_channel_two.id}
                for room in rooms:
                    if room.channel_id == temporary_channel_one.id:
                        assert room.room_name == temporary_channel_one.name
                    elif room.channel_id == temporary_channel_two.id:
                        assert room.room_name == temporary_channel_two.name
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_delete_temporary_room_by_channel(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$3, room_snowflake=$4
            ''', guild.id, member_one.id, temporary_channel_one.name, temporary_channel_one.id)

        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$3, room_snowflake=$4
            ''', guild.id, member_two.id, temporary_channel_two.name, temporary_channel_two.id)
        stage = await Stage.fetch_stage_by_guild_and_stage_name(guild, stage_name)
        assert stage.channel_id == stage_channel_one.id
        assert stage.channel_name == stage_name
        assert stage.expires_at == expires_at
        assert stage.guild_id == guild.id
        assert stage.initiator_id == member_one.id
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_stage_temporary_coordinator_ids_by_guild_and_stage_name(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO stage_coordinators (channel_id, discord_snowflake, guild_id, room_name)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
            ''', stage_channel_one.id, member_one.id, guild.id, stage_name, stage_channel_one.id, member_two.id, guild.id, stage_name)
        coordinator_ids = await Stage.fetch_stage_temporary_coordinator_ids_by_guild_and_stage_name(guild, stage_name)
        assert coordinator_ids == {111111111, 222222222}
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

@pytest.mark.asyncio
async def test_send_stage_ask_to_speak_message(bot_instance):
    pytest.skip("Requires Discord API mocking")

def test_update_stage_by_channel_and_temporary_coordinator_ids(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (channel_id, expires_at, guild_id, initiator_id, room_name)
                VALUES ($1, $2, $3, $4, $5)
            ''', stage_channel_one.id, expires_at, guild.id, member_one.id, stage_name)
            await conn.execute('''
                INSERT INTO stage_coordinators (channel_id, discord_snowflake, guild_id, room_name)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
            ''', stage_channel_one.id, member_three.id, guild.id, stage_name, stage_channel_one.id, member_four.id, guild.id, stage_name)
        stage = await Stage.fetch_stage_by_channel(stage_channel_one)
        await stage.update_stage_by_channel_and_temporary_coordinator_ids(stage_channel_one, temporary_coordinator_ids)
        updated_stage = await Stage.fetch_stage_by_channel(stage_channel_one)
        coordinator_ids = await Stage.fetch_stage_temporary_coordinator_ids_by_channel(stage_channel_one)
        assert coordinator_ids == {member_three.id, member_four.id}
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_update_stage_by_channel_initiator_and_temporary_coordinator_ids(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (channel_id, expires_at, guild_id, initiator_id, room_name)
                VALUES ($1, $2, $3, $4, $5)
            ''', stage_channel_one.id, expires_at, guild.id, member_one.id, stage_name)
            await conn.execute('''
                INSERT INTO stage_coordinators (channel_id, discord_snowflake, guild_id, room_name)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
            ''', stage_channel_one.id, member_three.id, guild.id, stage_name, stage_channel_one.id, member_four.id, guild.id, stage_name)
        stage = await Stage.fetch_stage_by_channel(stage_channel_one)
        await stage.update_stage_by_channel_initiator_and_temporary_coordinator_ids(stage_channel_one, member_one.id, temporary_coordinator_ids)
        updated_stage = await Stage.fetch_stage_by_channel(stage_channel_one)
        coordinator_ids = await Stage.fetch_stage_temporary_coordinator_ids_by_channel(stage_channel_one)
        assert coordinator_ids == {member_three.id, member_four.id}
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())