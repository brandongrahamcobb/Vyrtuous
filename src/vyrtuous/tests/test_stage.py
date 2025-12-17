from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.utils.stage import Stage
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
stage_channel_one = MockChannel(id=987654321, name="Test Stage Channel", guild=guild)
stage_channel_two = MockChannel(id=978653421, name="Test Stage Channel", guild=guild)
member_one = MockMember(id=111111111, name="Person 1")
member_two = MockMember(id=222222222, name="Person 2")
member_three = MockMember(id=333333333, name="Person 3")
member_four = MockMember(id=444444444, name="Person 4")
new_expires_at = datetime.now(timezone.utc) + timedelta(hours=2)
stage_name = "Test Stage Channel"
temporary_coordinator_ids = [member_one.id, member_two.id]

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
    
def test_fetch_stage_temporary_coordinator_ids_by_channel(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO stage_coordinators (channel_id, discord_snowflake, guild_id, room_name)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
            ''', stage_channel_one.id, member_one.id, guild.id, stage_name, stage_channel_two.id, 222222222, guild.id, stage_name)
        coordinator_ids = await Stage.fetch_stage_temporary_coordinator_ids_by_guild_and_stage_name(guild, stage_name)
        assert len(coordinator_ids) == 2
        assert coordinator_ids == {111111111, 222222222}
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_stage_by_channel(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (channel_id, expires_at, guild_id, initiator_id, room_name)
                VALUES ($1, $2, $3, $4, $5)
            ''', stage_channel_one.id, expires_at, guild.id, member_one.id, stage_name)
        stage = await Stage.fetch_stage_by_channel(stage_channel_one)
        assert stage.channel_id == stage_channel_one.id
        assert stage.channel_name == stage_name
        assert stage.expires_at == expires_at
        assert stage.guild_id == guild.id
        assert stage.initiator_id == member_one.id
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_stage_by_guild_and_stage_name(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (channel_id, expires_at, guild_id, initiator_id, room_name)
                VALUES ($1, $2, $3, $4, $5)
            ''', stage_channel_one.id, expires_at, guild.id, member_one.id, stage_name)
        stage = await Stage.fetch_stage_by_guild_and_stage_name(guild, stage_name)
        assert stage.channel_id == stage_channel_one.id
        assert stage.channel_name == stage_name
        assert stage.expires_at == expires_at
        assert stage.guild_id == guild.id
        assert stage.initiator_id == member_one.id
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_stage_temporary_coordinator_ids_by_guild_and_stage_name(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO stage_coordinators (channel_id, discord_snowflake, guild_id, room_name)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
            ''', stage_channel_one.id, member_one.id, guild.id, stage_name, stage_channel_one.id, member_two.id, guild.id, stage_name)
        coordinator_ids = await Stage.fetch_stage_temporary_coordinator_ids_by_guild_and_stage_name(guild, stage_name)
        assert len(coordinator_ids) == 2
        assert coordinator_ids == {111111111, 222222222}
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

@pytest.mark.asyncio
async def test_send_stage_ask_to_speak_message(bot_instance):
    pytest.skip("Requires Discord API mocking")

def test_update_stage_by_channel_and_temporary_coordinator_ids(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (channel_id, expires_at, guild_id, initiator_id, room_name)
                VALUES ($1, $2, $3, $4, $5)
            ''', stage_channel_one.id, expires_at, guild.id, member_one.id, stage_name)
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO stage_coordinators (channel_id, discord_snowflake, guild_id, room_name)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
            ''', stage_channel_one.id, member_three.id, guild.id, stage_name, stage_channel_one.id, member_four.id, guild.id, stage_name)
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT expires_at, initiator_id FROM active_stages
                WHERE channel_id = $1 AND guild_id = $2 AND room_name = $3
            ''', stage_channel_one.id, guild.id, stage_name)
        stage = Stage(row['expires_at'], stage_channel_one.id, stage_name, guild.id, row['initiator_id'])
        await stage.update_stage_by_channel_and_temporary_coordinator_ids(stage_channel_one, temporary_coordinator_ids)
        async with bot.db_pool.acquire() as conn:
            second_row = await conn.fetch('''
                SELECT discord_snowflake FROM stage_coordinators
                WHERE channel_id = $1 AND guild_id = $2
            ''', stage_channel_one.id, guild.id)
        if second_row:
            temporary_stage_coordinator_ids = {c['discord_snowflake'] for c in second_row}
        assert len(temporary_stage_coordinator_ids) == 2
        assert temporary_stage_coordinator_ids == {member_three.id, member_four.id}
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_update_stage_by_channel_initiator_and_temporary_coordinator_ids(bot_instance):
    async def inner():
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (channel_id, expires_at, guild_id, initiator_id, room_name)
                VALUES ($1, $2, $3, $4, $5)
            ''', stage_channel_one.id, expires_at, guild.id, member_one.id, stage_name)
            await conn.execute('''
                INSERT INTO stage_coordinators (channel_id, discord_snowflake, guild_id, room_name)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
            ''', stage_channel_one.id, member_three.id, guild.id, stage_name, stage_channel_one.id, member_four.id, guild.id, stage_name)
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT expires_at, initiator_id FROM active_stages
                WHERE channel_id = $1 AND guild_id = $2 AND room_name = $3
            ''', stage_channel_one.id, guild.id, stage_name)
        stage = Stage(row['expires_at'], stage_channel_one.id, stage_name, guild.id, row['initiator_id'])
        await stage.update_stage_by_channel_initiator_and_temporary_coordinator_ids(stage_channel_one, member_one.id, temporary_coordinator_ids)
        async with bot.db_pool.acquire() as conn:
            second_row = await conn.fetch('''
                SELECT discord_snowflake FROM stage_coordinators
                WHERE channel_id = $1 AND guild_id = $2
            ''', stage_channel_one.id, guild.id)
        if second_row:
            temporary_stage_coordinator_ids = {c['discord_snowflake'] for c in second_row}
        assert temporary_stage_coordinator_ids == {member_three.id, member_four.id}
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())