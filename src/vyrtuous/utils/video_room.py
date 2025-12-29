''' video_rooms.py A utility module for managing video rooms in the Vyrtuous Discord bot.

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
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot

class VideoRoom:
        
    def __init__(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        self.bot = DiscordBot.get_instance()
        self.channel_mention = f"<#{channel_snowflake}>"
        self.channel_snowflake = channel_snowflake
        self.guild_snowflake = guild_snowflake
        self.is_video_room: Optional[bool] = True

    async def create(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO video_rooms (channel_snowflake, guild_snowflake)
                VALUES ($1, $2)
                ON CONFLICT (channel_snowflake, guild_snowflake)
                DO NOTHING
            ''', self.channel_snowflake, self.guild_snowflake)
            
    @classmethod
    async def fetch_all(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, guild_snowflake FROM video_rooms
            ''')
        video_rooms = []
        if rows:
            for row in rows:
                video_rooms.append(VideoRoom(channel_snowflake=row['channel_snowflake'], guild_snowflake=row['guild_snowflake']))
        return video_rooms

    @classmethod
    async def fetch_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT channel_snowflake FROM video_rooms WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
        video_room = []
        if row:
            video_room.append(VideoRoom(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake))
        return video_room
            
    @classmethod
    async def delete_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM video_rooms WHERE channel_snowflake = $1 AND guild_snowflake = $2
            ''', channel_snowflake, guild_snowflake)
    
    @classmethod
    async def fetch_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake FROM video_rooms WHERE guild_snowflake=$1
            ''', guild_snowflake)
        video_rooms = []
        if rows:
            for row in rows:
                video_rooms.append(VideoRoom(channel_snowflake=row['channel_snowflake'], guild_snowflake=guild_snowflake))
        return video_rooms
