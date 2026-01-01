''' temporary_rooms.py A utility module for managing temporary rooms in the Vyrtuous Discord bot.

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

class TemporaryRoom:
        
    def __init__(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int], room_name: Optional[str]):
        self.bot = DiscordBot.get_instance()
        self.channel_mention = f"<#{channel_snowflake}>"
        self.channel_snowflake = channel_snowflake
        self.guild_snowflake = guild_snowflake
        self.is_temp_room: Optional[bool] = True
        self.member_mention = f"<@{member_snowflake}>"
        self.member_snowflake = member_snowflake
        self.room_name = room_name

    async def create(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (channel_snowflake, guild_snowflake, member_snowflake, room_name)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT DO NOTHING
            ''', self.channel_snowflake, self.guild_snowflake, self.member_snowflake, self.room_name)

    @classmethod
    async def fetch_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT channel_snowflake, member_snowflake, room_name
                FROM temporary_rooms
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
        if not row:
            return None
        return TemporaryRoom(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=row['member_snowflake'], room_name=row['room_name'])
            
    @classmethod
    async def fetch_by_guild_and_member(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, guild_snowflake, member_snowflake, room_name, updated_at
                FROM temporary_rooms
                WHERE guild_snowflake=$1 AND member_snowflake=$2
            ''', guild_snowflake, member_snowflake)
        temporary_rooms = []
        if rows:
            for row in rows:
                temporary_rooms.append(TemporaryRoom(hannel_snowflake=row['channel_snowflake'], guild_snowflake=guild_snowflake, member_snowflake=member_snowflake, room_name=row['room_name']))
        return temporary_rooms
            
    @classmethod
    async def fetch_by_guild_and_room_name(cls, guild_snowflake: Optional[int], room_name: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT channel_snowflake, created_at, guild_snowflake, member_snowflake, room_name, updated_at
                FROM temporary_rooms
                WHERE guild_snowflake=$1 AND room_name=$2
            ''', guild_snowflake, room_name)
        if not row:
            return None
        return TemporaryRoom(channel_snowflake=row['channel_snowflake'], guild_snowflake=guild_snowflake, member_snowflake=row['member_snowflake'], room_name=row['room_name'])
    
    @classmethod
    async def update_owner(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE temporary_rooms
                SET member_snowflake=$3
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 
            ''', channel_snowflake, guild_snowflake, member_snowflake)
            
    @classmethod
    async def update_by_source_and_target(cls, guild_snowflake: Optional[int], source_channel_snowflake: Optional[int], room_name: Optional[str], target_channel_snowflake: Optional[int], ):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE temporary_rooms
                SET room_name=$3, channel_snowflake=$4
                WHERE guild_snowflake=$1 AND channel_snowflake=$2
            ''', guild_snowflake, target_channel_snowflake, room_name, source_channel_snowflake)
            
    @classmethod
    async def delete_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM temporary_rooms
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
    
    @classmethod
    async def fetch_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, guild_snowflake, member_snowflake, room_name, updated_at
                FROM temporary_rooms WHERE guild_snowflake=$1 ORDER BY room_name
            ''', guild_snowflake)
        temporary_rooms = []
        if rows:
            for row in rows:
                temporary_rooms.append(TemporaryRoom(channel_snowflake=row['channel_snowflake'], guild_snowflake=guild_snowflake, member_snowflake=row['member_snowflake'], room_name=row['room_name']))
        return temporary_rooms
            
    @classmethod
    async def fetch_all(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, guild_snowflake, member_snowflake, room_name, updated_at
                FROM temporary_rooms
            ''')
        temporary_rooms = []
        if rows:
            for row in rows:
                temporary_rooms.append(TemporaryRoom(channel_snowflake=row['channel_snowflake'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], room_name=row['room_name']))
        return temporary_rooms
