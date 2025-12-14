
from typing import Optional

from vyrtuous.bot.discord_bot import DiscordBot
import discord
import asyncpg

class TemporaryRoom:
    
        
    def __init__(self, guild: discord.Guild, channel: discord.abc.GuildChannel, room_owner: discord.Member):
        self.channel = channel
        self.bot = DiscordBot.get_instance()
        self.guild = guild
        self.is_temp_room: Optional[bool] = True
        self.room_name: Optional[str] = channel.name
        self.room_owner = room_owner

    def load_channel(self, channel_obj: discord.abc.GuildChannel):
        if channel_obj and isinstance(channel_obj, discord.abc.GuildChannel):
            self.channel = channel_obj
        else:
            raise ValueError(f"Invalid channel.")
        
    def load_guild(self, guild_obj: discord.Guild):
        if guild_obj and isinstance(guild_obj, discord.Guild):
            self.guild = guild_obj
        else:
            raise ValueError(f"Invalid guild.")
        
    def load_room_name(self, room_name_str: Optional[str]):
        if room_name_str and room_name_str.strip():
            self.room_name = room_name_str
        else:
            raise ValueError(f"Invalid room name.")
        
    def load_room_owner(self, room_owner_obj: discord.Member):
        if room_owner_obj:
            self.room_owner = room_owner_obj
        else:
            raise ValueError(f"Invalid room owner.")

    async def insert_into_temporary_rooms(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, room_name, owner_snowflake, room_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_snowflake, room_name)
                DO UPDATE SET owner_snowflake=$3, room_snowflake=$4
            ''', self.guild.id, self.channel.name, self.room_owner.id, self.channel.id)

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
            return TemporaryRoom(guild=channel.guild, channel=channel, room_owner=member)
            
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
                channel = guild.get_channel(row['room_snowflake'])
                temporary_rooms.append(TemporaryRoom(guild=guild, channel=channel, room_owner=member))
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
            channel = guild.get_channel(room['room_snowflake']) # TODO: NoneType after the channel is deleted and can't be used as a reference anymore.
            member = guild.get_member(room['owner_snowflake'])
            return TemporaryRoom(guild=guild, channel=channel, room_owner=member)
            
    async def update_temporary_room_owner_snowflake(self, member: discord.Member):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                'UPDATE temporary_rooms SET owner_snowflake=$1 WHERE guild_snowflake=$2 AND room_snowflake=$3 AND room_name=$4',
                member.id, self.guild.id, self.channel.id, self.channel.name
            )
            
    async def update_temporary_room_name_and_room_snowflake(self, channel: discord.abc.GuildChannel):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                'UPDATE temporary_rooms SET room_name=$3, room_snowflake=$4 WHERE guild_snowflake=$1 AND room_name=$2',
                self.guild.id, self.channel.name, channel.name, channel.id
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
                temporary_rooms.append(TemporaryRoom(guild=guild, channel=channel, room_owner=member))
            return temporary_rooms
            
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
                    channel = guild.get_channel(row['room_snowflake'])
                    temporary_rooms.append(TemporaryRoom(guild=guild, channel=channel, room_owner=member))
                guilds[guild] = temporary_rooms
            if not guilds:
                return None
            return guilds
