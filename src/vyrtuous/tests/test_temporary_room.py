
# ''' test_temporary_room.py The purpose of this program is to provide the tess for the TemporaryRoom module.
#     Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
# '''
# from vyrtuous.bot.discord_bot import DiscordBot
# from vyrtuous.tests.black_box.test_suite import *
# from vyrtuous.utils.temporary_room import TemporaryRoom
# import asyncio
# import asyncpg
# import pytest

# def test_insert_into_temporary_rooms(bot_instance):
#     async def inner():
#         bot = DiscordBot.get_instance()
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('''
#                 INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
#                 VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
#                 ON CONFLICT (guild_snowflake, room_name)
#                 DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
#             ''', guild.id, member_one.id, channel_one.name, channel_one.id, guild.id, member_two.id, channel_two.name, channel_two.id)
#         rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild(guild)
#         assert len(rooms) == 2
#         channel_ids = {room.channel_id for room in rooms}
#         assert channel_ids == {channel_one.id, channel_two.id}
#         is_temp_rooms = {room.is_temp_room for room in rooms}
#         assert is_temp_rooms == {True, True}
#         room_owners = {room.room_owner.id for room in rooms}
#         assert room_owners == {member_one.id, member_two.id}
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#     asyncio.get_event_loop().run_until_complete(inner())

# def test_fetch_temporary_room_by_guild_and_member(bot_instance):
#     async def inner():
#         bot = DiscordBot.get_instance()
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('''
#                 INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
#                 VALUES ($1, $2, $3, $4)
#                 ON CONFLICT (guild_snowflake, room_name)
#                 DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
#             ''', guild.id, member_one.id, channel_one.name, channel_one.id)
#         rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild_and_member(guild, member_one)
#         assert rooms is not None
#         assert len(rooms) == 1
#         for room in rooms:
#             assert room is not None
#             assert room.channel == channel_one
#             assert room.channel_id == channel_one.id
#             assert room.guild == guild
#             assert room.is_temp_room == True
#             assert room.room_owner == member_one
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#     asyncio.get_event_loop().run_until_complete(inner())

# def test_fetch_temporary_room_by_guild_and_room_name(bot_instance):
#     async def inner():
#         bot = DiscordBot.get_instance()
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('''
#                 INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
#                 VALUES ($1, $2, $3, $4)
#                 ON CONFLICT (guild_snowflake, room_name)
#                 DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
#             ''', guild.id, member_one.id, channel_one.name, channel_one.id)
#         room = await TemporaryRoom.fetch_temporary_room_by_guild_and_room_name(guild, channel_one.name)
#         assert room is not None
#         assert room.channel == channel_one
#         assert room.channel_id == channel_one.id
#         assert room.guild == guild
#         assert room.is_temp_room == True
#         assert room.room_owner == member_one
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#     asyncio.get_event_loop().run_until_complete(inner())

# def test_fetch_all_guilds_with_temporary_rooms(bot_instance):
#     async def inner():
#         bot = DiscordBot.get_instance()
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('''
#                 INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
#                 VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
#                 ON CONFLICT (guild_snowflake, room_name)
#                 DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
#             ''', guild.id, member_one.id, channel_one.name, channel_one.id, guild.id, member_two.id, channel_two.name, channel_two.id)
#         guilds = await TemporaryRoom.fetch_all_guilds_with_temporary_rooms()
#         assert guild in guilds
#         for g in guilds:
#             if g.id == 123456789:
#                 rooms = guilds[g]
#                 assert len(rooms) == 2
#                 channel_ids = {room.channel_id for room in rooms}
#                 assert channel_ids == {channel_one.id, channel_two.id}
#                 channels = {room.channel for room in rooms}
#                 assert channels == {channel_one, channel_two}
#                 is_temp_rooms = {room.is_temp_room for room in rooms}
#                 assert is_temp_rooms == {True, True}
#                 room_owners = {room.room_owner for room in rooms}
#                 assert room_owners == {member_one, member_two}
#                 room_snowflakes = {room.channel_id for room in rooms}
#                 assert room_snowflakes == {channel_one.id, channel_two.id}
#                 for room in rooms:
#                     if room.channel_id == channel_one.id:
#                         assert room.room_name == channel_one.name
#                     elif room.channel_id == channel_two.id:
#                         assert room.room_name == channel_two.name
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#     asyncio.get_event_loop().run_until_complete(inner())
    
# def test_delete_temporary_room_by_channel(bot_instance):
#     async def inner():
#         bot = DiscordBot.get_instance()
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('''
#                 INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
#                 VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
#                 ON CONFLICT (guild_snowflake, room_name)
#                 DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
#             ''', guild.id, member_one.id, channel_one.name, channel_one.id, guild.id, member_two.id, channel_two.name, channel_two.id)
#         await TemporaryRoom.delete_temporary_room_by_channel(channel_one)
#         async with bot.db_pool.acquire() as conn:
#             room_one = await conn.fetchrow(
#                 'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2',
#                 guild.id, channel_one.id
#             )
#         async with bot.db_pool.acquire() as conn:
#             room_two = await conn.fetchrow(
#                 'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2',
#                 guild.id, channel_two.id
#             )
#         assert room_one is None
#         assert room_two['room_name'] == channel_two.name
#         assert room_two['owner_snowflake'] == member_two.id
#         assert room_two['room_snowflake'] == channel_two.id
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#     asyncio.get_event_loop().run_until_complete(inner())

# def test_update_temporary_room_owner_snowflake(bot_instance):
#     async def inner():
#         bot = DiscordBot.get_instance()
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('''
#                 INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
#                 VALUES ($1, $2, $3, $4)
#                 ON CONFLICT (guild_snowflake, room_name)
#                 DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
#             ''', guild.id, member_one.id, channel_one.name, channel_one.id)
#         room = TemporaryRoom(guild, channel_one.id, member_one)
#         await room.update_temporary_room_owner_snowflake(member_two)
#         async with bot.db_pool.acquire() as conn:
#             room_one = await conn.fetchrow(
#                 'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2',
#                 guild.id, channel_one.id
#             )
#         assert room_one['room_snowflake'] == channel_one.id
#         assert room_one['room_name'] == channel_one.name
#         assert room_one['owner_snowflake'] == member_two.id
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#     asyncio.get_event_loop().run_until_complete(inner())

# def test_update_temporary_room_name_and_room_snowflake(bot_instance):
#     async def inner():
#         bot = DiscordBot.get_instance()
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('''
#                 INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
#                 VALUES ($1, $2, $3, $4)
#                 ON CONFLICT (guild_snowflake, room_name)
#                 DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
#             ''', guild.id, member_one.id, channel_one.name, channel_one.id)
#         room = TemporaryRoom(guild, channel_one.id, member_one)
#         await room.update_temporary_room_name_and_room_snowflake(channel_two, channel_two.name)
#         async with bot.db_pool.acquire() as conn:
#             room_one = await conn.fetchrow(
#                 'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2',
#                 guild.id, channel_one.id
#             )
#         assert room_one['room_snowflake'] == channel_one.id
#         assert room_one['room_name'] == channel_one.name
#         assert room_one['owner_snowflake'] == member_one.id
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#     asyncio.get_event_loop().run_until_complete(inner())

# def test_fetch_temporary_rooms_by_guild(bot_instance):
#     async def inner():
#         bot = DiscordBot.get_instance()
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('''
#                 INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
#                 VALUES ($1, $2, $3, $4), ($1, $5, $6, $7)
#                 ON CONFLICT (guild_snowflake, room_name)
#                 DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
#             ''', guild.id, member_one.id, channel_one.name, channel_one.id, member_two.id, channel_two.name, channel_two.id)
#         rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild(guild)
#         assert rooms is not None
#         assert len(rooms) == 2
#         channel_ids = {room.channel_id for room in rooms}
#         assert channel_ids == {channel_one.id, channel_two.id}
#         channels = {room.channel for room in rooms}
#         assert channels == {channel_one, channel_two}
#         is_temp_rooms = {room.is_temp_room for room in rooms}
#         assert is_temp_rooms == {True, True}
#         room_owners = {room.room_owner.id for room in rooms}
#         assert room_owners == {member_one.id, member_two.id}
#         room_snowflakes = {room.channel_id for room in rooms}
#         assert room_snowflakes == {channel_one.id, channel_two.id}
#         for room in rooms:
#             if room.channel_id == channel_one.id:
#                 assert room.room_name == channel_one.name
#             elif room.channel_id == channel_two.id:
#                 assert room.room_name == channel_two.name
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#     asyncio.get_event_loop().run_until_complete(inner())

# def test_fetch_temporary_room_by_channel(bot_instance):
#     async def inner():
#         bot = DiscordBot.get_instance()
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('''
#                 INSERT INTO temporary_rooms (guild_snowflake, owner_snowflake, room_name, room_snowflake)
#                 VALUES ($1, $2, $3, $4), ($5, $6, $7, $8)
#                 ON CONFLICT (guild_snowflake, room_name)
#                 DO UPDATE SET owner_snowflake=$2, room_snowflake=$4
#             ''', guild.id, member_one.id, channel_one.name, channel_one.id, guild.id, member_two.id, channel_two.name, channel_two.id)
#         room_one = await TemporaryRoom.fetch_temporary_room_by_channel(channel_one)
#         assert room_one is not None
#         assert room_one.channel_id == channel_one.id
#         assert room_one.guild == guild
#         assert room_one.is_temp_room == True
#         assert room_one.room_name == channel_one.name
#         assert room_one.room_owner == member_one
#         room_two = await TemporaryRoom.fetch_temporary_room_by_channel(channel_two)
#         assert room_two is not None
#         assert room_two.channel_id == channel_two.id
#         assert room_two.guild == guild
#         assert room_two.is_temp_room == True
#         assert room_two.room_name == channel_two.name
#         assert room_two.room_owner == member_two
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake = 123456789')
#     asyncio.get_event_loop().run_until_complete(inner())