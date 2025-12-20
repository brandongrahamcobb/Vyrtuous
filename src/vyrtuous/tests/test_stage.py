
''' test_stage.py The purpose of this program is to provide the tests for the Stage module.
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
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.tests.test_suite import *
from vyrtuous.utils.stage import Stage
import asyncio
import asyncpg
from datetime import datetime, timedelta, timezone
import inspect
import os
import pytest
import pytest_asyncio
import vyrtuous

stage_name = "Test Stage Channel"
temporary_coordinator_ids = [self_member_obj.id, dummy_member_obj.id]
    
def test_fetch_stage_temporary_coordinator_ids_by_guild_id_and_channel_name(bot_instance):
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
            ''', channel_one.id, member_one.id, guild.id, stage_name, channel_two.id, 222222222, guild.id, stage_name)
        coordinator_ids = await Stage.fetch_stage_temporary_coordinator_ids_by_guild_id_and_channel_name(guild_id=guild.id, channel_name=stage_name)
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
            ''', channel_one.id, expires_at, guild.id, member_one.id, stage_name)
        stage = await Stage.fetch_stage_by_guild_id_channel_id_and_channel_name(guild_id=guild.id, channel_id=channel_one.id, channel_name=channel_one.name)
        assert stage.channel_id == channel_one.id
        assert stage.channel_name == stage_name
        assert stage.expires_at == expires_at
        assert stage.guild_id == guild.id
        assert stage.initiator_id == member_one.id
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_stage_by_guild_id_and_channel_name(bot_instance):
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
            ''', channel_one.id, expires_at, guild.id, member_one.id, stage_name)
        stage = await Stage.fetch_stage_by_guild_id_and_channel_name(guild=guild.id, channel_name=stage_name)
        assert stage.channel_id == channel_one.id
        assert stage.channel_name == stage_name
        assert stage.expires_at == expires_at
        assert stage.guild_id == guild.id
        assert stage.initiator_id == member_one.id
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())

def test_fetch_stage_temporary_coordinator_ids_channel_id(bot_instance):
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
            ''', channel_one.id, member_one.id, guild.id, stage_name, channel_one.id, member_two.id, guild.id, stage_name)
        coordinator_ids = await Stage.fetch_stage_temporary_coordinator_ids_by_channel_id(channel_id=channel_one.id)
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
            ''', channel_one.id, expires_at, guild.id, member_one.id, stage_name)
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO stage_coordinators (channel_id, discord_snowflake, guild_id, room_name)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
            ''', channel_one.id, member_three.id, guild.id, stage_name, channel_one.id, member_four.id, guild.id, stage_name)
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT expires_at, initiator_id FROM active_stages
                WHERE channel_id = $1 AND guild_id = $2 AND room_name = $3
            ''', channel_one.id, guild.id, stage_name)
        stage = Stage(row['expires_at'], channel_one.id, stage_name, guild.id, row['initiator_id'])
        await stage.update_stage_by_channel_id_name(channel_id=channel_one.id, channel_name=channel_one.name)
        async with bot.db_pool.acquire() as conn:
            second_row = await conn.fetch('''
                SELECT discord_snowflake FROM stage_coordinators
                WHERE channel_id = $1 AND guild_id = $2
            ''', channel_one.id, guild.id)
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
            ''', channel_one.id, expires_at, guild.id, member_one.id, stage_name)
            await conn.execute('''
                INSERT INTO stage_coordinators (channel_id, discord_snowflake, guild_id, room_name)
                VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
            ''', channel_one.id, member_three.id, guild.id, stage_name, channel_one.id, member_four.id, guild.id, stage_name)
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT expires_at, initiator_id FROM active_stages
                WHERE channel_id = $1 AND guild_id = $2 AND room_name = $3
            ''', channel_one.id, guild.id, stage_name)
        stage = Stage(row['expires_at'], channel_one.id, stage_name, guild.id, row['initiator_id'])
        await stage.update_stage_by_channel_id_name_initiator_id_and_temporary_coordinator_ids(channel_id=channel_one.id, channel_name=channel_one.name, stage_initiator_id=member_one.id, temporary_stage_coordinator_ids=temporary_coordinator_ids)
        async with bot.db_pool.acquire() as conn:
            second_row = await conn.fetch('''
                SELECT discord_snowflake FROM stage_coordinators
                WHERE channel_id = $1 AND guild_id = $2
            ''', channel_one.id, guild.id)
        if second_row:
            temporary_stage_coordinator_ids = {c['discord_snowflake'] for c in second_row}
        assert temporary_stage_coordinator_ids == {member_three.id, member_four.id}
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM stage_coordinators WHERE guild_id = 123456789')
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_stages WHERE guild_id = 123456789')
    asyncio.get_event_loop().run_until_complete(inner())