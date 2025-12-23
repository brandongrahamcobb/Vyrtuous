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
import discord
import asyncpg

class TemporaryRoom:
        
    def __init__(self, guild: discord.Guild, channel_id: Optional[str], room_name: Optional[str], room_owner: discord.Member):
        self.bot = DiscordBot.get_instance()
        self.channel_id: Optional[int] = channel_id
        self.channel_mention = f"<@{channel_id}>"
        self.guild = guild
        self.is_temp_room: Optional[bool] = True
        self.room_name = room_name
        self.room_owner = room_owner

    async def insert_into_temporary_rooms(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, room_name, owner_snowflake, room_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$3, room_snowflake=$4
            ''', self.guild.id, self.room_name, self.room_owner.id, self.channel_id)

    @classmethod
    async def fetch_temporary_room_by_channel(cls, channel: discord.abc.GuildChannel):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            room = await conn.fetchrow(
                'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2',
                channel.guild.id, channel.id
            )
            if not room:
                return None
            member = channel.guild.get_member(room['owner_snowflake'])
            return TemporaryRoom(guild=channel.guild, channel_id=channel.id, room_name=channel.name, room_owner=member)
            
    @classmethod
    async def fetch_temporary_rooms_by_guild_and_member(cls, guild: discord.Guild, member: discord.Member):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND owner_snowflake=$2',
                guild.id, member.id
            )
            if not rows:
                return None
            temporary_rooms = []
            for row in rows:
                temporary_rooms.append(TemporaryRoom(guild=guild, channel_id=row['room_snowflake'], room_name=row['room_name'], room_owner=member))
            return temporary_rooms
            
    @classmethod
    async def fetch_temporary_room_by_guild_and_room_name(cls, guild: discord.Guild, room_name: str):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            room = await conn.fetchrow(
                'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_name=$2',
                guild.id, room_name
            )
            if not room:
                return None
            member = guild.get_member(room['owner_snowflake'])
            return TemporaryRoom(guild=guild, channel_id=room['room_snowflake'], room_name=room['room_name'], room_owner=member)
            
    async def update_temporary_room_owner_snowflake(self, member: discord.Member):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                'UPDATE temporary_rooms SET owner_snowflake=$1 WHERE guild_snowflake=$2 AND room_snowflake=$3 AND room_name=$4',
                member.id, self.guild.id, self.channel_id, self.room_name
            )
            
    async def update_temporary_room_name_and_room_snowflake(self, channel: discord.abc.GuildChannel, room_name: Optional[str]):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                'UPDATE temporary_rooms SET room_name=$3, room_snowflake=$4 WHERE guild_snowflake=$1 AND room_name=$2',
                self.guild.id, channel.name, room_name, channel.id
            )
            
    @classmethod
    async def delete_temporary_room_by_channel(cls, channel: discord.abc.GuildChannel):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                'DELETE FROM temporary_rooms WHERE guild_snowflake = $1 AND room_snowflake = $2',
                channel.guild.id, channel.id
            )
    
    @classmethod
    async def fetch_temporary_rooms_by_guild(cls, guild: discord.Guild):
        try:
            bot = DiscordBot.get_instance()
            async with bot.db_pool.acquire() as conn:
                rows = await conn.fetch(
                    'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 ORDER BY room_name',
                    guild.id
                )
                if not rows:
                    return None
                temporary_rooms = []
                for row in rows:
                    member = guild.get_member(row['owner_snowflake'])
                    channel = guild.get_channel(row['room_snowflake'])
                    temporary_rooms.append(TemporaryRoom(guild=guild, channel_id=row['room_snowflake'], room_name=row['room_name'], room_owner=member))
                return temporary_rooms
            raise Exception('No temporary rooms found for {guild.name}.')
        except Exception:
            raise
            
    @classmethod
    async def fetch_all_guilds_with_temporary_rooms(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            guilds = {}
            for guild in bot.guilds:
                if not guild:
                    return None
                rows = await conn.fetch(
                    'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1',
                    guild.id
                )
                temporary_rooms = []
                for row in rows:
                    member = guild.get_member(row['owner_snowflake'])
                    temporary_rooms.append(TemporaryRoom(guild=guild, channel_id=row['room_snowflake'], room_name=row['room_name'], room_owner=member))
                guilds[guild] = temporary_rooms
            if not guilds:
                return None
            return guilds
