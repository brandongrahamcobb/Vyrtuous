
# ''' test_cap.py The purpose of this program is to provide the tests for the Cap module.
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
# from datetime import timedelta
# from vyrtuous.tests.black_box.test_suite import *
# from vyrtuous.utils.cap import Cap
# import asyncio
# import pytest

# duration_seconds_one = timedelta(hours=7).total_seconds()
# moderation_type_one = 'mute'
# duration_seconds_two = timedelta(hours=7).total_seconds()
# moderation_type_two = 'ban'

# def test_fetch_caps_for_channel(bot_instance):
#     async def inner():
#         bot = DiscordBot.get_instance()
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('DELETE FROM active_caps WHERE guild_id = 123456789')
#         async with bot.db_pool.acquire() as conn:
#             await conn.execute('''
#                 INSERT INTO active_caps (channel_id, duration_seconds, guild_id, moderation_type, room_name)
#                 VALUES ($1, $2, $3, $4, $5), ($1, $6, $3, $7, $8)
#                 ON CONFLICT (channel_id, guild_id, moderation_type, room_name)
#                 DO UPDATE SET duration_seconds=$2, moderation_type=$4
#             ''', channel_one.id, duration_seconds_one, guild.id, moderation_type_one, channel_one.name, duration_seconds_two, moderation_type_two, channel_one.name)
#         caps = await Cap.fetch_caps_for_channel(guild.id, channel_one.id)
#         assert len(caps) == 2
#         assert caps[0][0] == duration_seconds_one
#         assert caps[0][1] == moderation_type_one
#         assert caps[1][0] == duration_seconds_two
#         assert caps[1][1] == moderation_type_two
#     asyncio.get_event_loop().run_until_complete(inner())