
from typing import Optional

from vyrtuous.bot.discord_bot import DiscordBot
import discord
import asyncpg

class TemporaryRoom:
    
        
    def __init__(self, guild: discord.Guild, channel: discord.abc.GuildChannel, room_owner: discord.User):
        self.channel = channel
        self.bot = DiscordBot.get_instance()
        self.guild = guild
        self.is_temp_room: Optional[bool] = True
        self.room_name: Optional[str] = channel.name
        self.room_owner = room_owner

    def get_temporary_room_name_and_channel_snowflake_by_command(self, alias_category: Optional[str], command_name: Optional[str]) -> dict[str, str]:
        temporary_room_aliases_by_category = self.get_temporary_room_aliases_by_category(alias_category)
        temporary_room_alias_channel_snowflake_and_room_name_pair = temporary_room_aliases_by_category.get(command_name)
        return temporary_room_alias_channel_snowflake_and_room_name_pair
        
    def get_temporary_room_aliases_by_category(self, alias_category: Optional[str]) -> dict[str, dict[str, str]]:
        temporary_room_aliases_by_category = bot.command_aliases.get(
            self.guild.id,
            self.temporary_room_alias.command_aliases.default_factory()
        ).get('channel_aliases', {}).get(alias_category, {})
        return temporary_room_aliases_by_category
    
    def get_all_temporary_room_aliases(self) -> dict[str, dict[str, dict[str, str]]]:
        aliases = self.temporary_room_alias.command_aliases.get(self.guild.id, {}).get('channel_aliases', {})
        return aliases
        
    def get_temporary_room_channel_snowflake(self, temporary_room_name_and_channel_snowflake_pair: dict[str, str]) -> Optional[str]:
        temporary_room_channel_snowflake = temporary_room_name_and_channel_snowflake_pair.get('channel_id')
        return temporary_room_channel_snowflake
        
    def get_temporary_room_name(self, temporary_room_name_and_channel_snowflake_pair: dict[str, str]) -> Optional[str]:
        temporary_room_name = temporary_room_name_and_channel_snowflake_pair.get('room_name')
        return temporary_room_name

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
        
    def load_room_owner(self, room_owner_obj: discord.User):
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
    async def fetch_temporary_room_by_channel(self, guild: discord.Guild, channel: discord.abc.GuildChannel):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            room = await conn.fetch(
                'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2',
                guild.id, channel.id
            )
            member = guild.get_member(room['owner_snowflake'])
            return TemporaryRoom(guild=guild, channel=channel, room_owner=member)
            
    @classmethod
    async def fetch_temporary_rooms_by_member(self, guild: discord.Guild, member: discord.User):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND owner_snowflake=$2',
                guild.id, member.id
            )
            temporary_rooms = []
            for row in rows:
                channel = guild.get_channel(row['room_snowflake'])
                temporary_rooms.append(TemporaryRoom(guild=guild, channel=channel, room_owner=member))
            return temporary_rooms
            
    async def update_temporary_room(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                'UPDATE temporary_rooms SET owner_snowflake=$1 WHERE guild_snowflake=$2 AND room_snowflake=$3 AND room_name=$4',
                self.room_owner.id, self.guild.id, self.channel.id, self.channel.name
            )
            
    @classmethod
    async def delete_temporary_room(self, guild: discord.Guild, channel: discord.abc.GuildChannel):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                'DELETE FROM temporary_rooms WHERE guild_snowflake = $1 AND room_snowflake = $2',
                guild.id, channel.id
            )
    
    @classmethod
    async def fetch_temporary_rooms_by_guild(self, guild: discord.Guild):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT owner_snowflake, room_name, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 ORDER BY room_name',
                guild.id
            )
            temporary_rooms = []
            for row in rows:
                member = guild.get_member(row['owner_snowflake'])
                channel = guild.get_channel(row['room_snowflake'])
                temporary_rooms.append(TemporaryRoom(guild=guild, channel=channel, room_owner=member))
            return temporary_rooms
            
    @classmethod
    async def fetch_all_guilds_with_temporary_rooms(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            guilds = {}
            for guild in bot.guilds:
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
            return guilds
