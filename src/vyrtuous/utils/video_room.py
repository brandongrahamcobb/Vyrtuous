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
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.setup_logging import logger
import asyncio

class VideoRoom:

    video_rooms = []
        
    def __init__(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        self.bot = DiscordBot.get_instance()
        self.channel_mention = f"<#{channel_snowflake}>"
        self.channel_snowflake = channel_snowflake
        self.cooldowns = {}
        self.emojis = Emojis()
        self.guild_snowflake = guild_snowflake
        self.is_video_room: Optional[bool] = True
        self.video_tasks = {}

    @classmethod
    async def enforce_video(cls, member, channel, delay):
        await asyncio.sleep(delay)
        if not member.voice:
            return
        if member.voice.channel != channel:
            return
        if member.voice.self_video:
            return
        try:
            await member.move_to(None)
        except Exception as e:
            logger.info("Unable to enforce video by kicking the user.")
        try:
            await member.send(f"{cls.emojis.get_random_emoji()} You were kicked from {channel.mention} because your video feed stopped. {channel.mention} is a video-only channel.")
        except Exception as e:
            logger.info("Unable to send a message to enforce video.")

    @classmethod
    def cancel_task(cls, key):
        task = cls.video_tasks.pop(key, None)
        if task:
            task.cancel()

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
