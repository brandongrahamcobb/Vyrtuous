from typing import defaultdict, Optional
from vyrtuous.bot.discord_bot import DiscordBot
import discord


class Alias:
    
    def __init__(self, guild_id: Optional[int], channel_id: Optional[int], alias_type: Optional[str], alias_name: Optional[str], role_id: Optional[int]):
        self.bot = DiscordBot.get_instance()
        self.alias_type = alias_type
        self.alias_name = alias_name
        self.channel_id = channel_id
        self.guild_id = guild_id
        self.alias_cog = self.bot.get_cog("Aliases")
        self.handlers = {
            'ban': self.alias_cog.handle_ban_alias,
            'cow': self.alias_cog.handle_cow_alias,
            'uncow': self.alias_cog.handle_uncow_alias,
            'unban': self.alias_cog.handle_unban_alias,
            'flag': self.alias_cog.handle_flag_alias,
            'unflag': self.alias_cog.handle_unflag_alias,
            'mute': self.alias_cog.handle_voice_mute_alias,
            'unmute': self.alias_cog.handle_unmute_alias,
            'tmute': self.alias_cog.handle_text_mute_alias,
            'untmute': self.alias_cog.handle_untextmute_alias,
            'role': self.alias_cog.handle_role_alias,
            'unrole': self.alias_cog.handle_unrole_alias
        }
        self.handler = self.handlers[alias_type]
        self.role_id = role_id
        
    @classmethod
    def format_aliases_by_channel(self, aliases, channel: discord.abc.GuildChannel) -> list[str]:
        if not aliases:
            return []
        lines = []
        lines.append(f'**Aliases in {channel.mention}**')
        for alias in aliases:
            if alias.role_id:
                lines.append(f'`{alias.alias_name}` → <@&{alias.role_id}>')
            else:
                lines.append(f'`{alias.alias_name}`')
        return lines
        
    @classmethod
    def format_all_aliases(self, aliases) -> list[str]:
        if not aliases:
            return []
        grouped = defaultdict(list)
        lines = []
        for alias in aliases:
            grouped[alias.channel_id].append(alias)
        for channel, channel_aliases in grouped.items():
            lines.append(f'**Aliases in {channel.mention}**')
            for alias in channel_aliases:
                if alias.role_id:
                    lines.append(f'`{alias.alias_name}` → <@&{alias.role_id}>')
                else:
                    lines.append(f'`{alias.alias_name}`')
        return lines

    def load_alias_name(self, alias_name: str):
        if not alias_name or not isinstance(alias_name, str):
            raise ValueError("Invalid alias name.")
        if self.alias_type is None:
            raise ValueError("Alias type must be set before alias name.")
        self.alias_name = alias_name

    async def insert_into_command_aliases(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id, role_id)
                VALUES ($1, $2, $3, $4, $5)
            ''', self.guild_id, self.alias_type, self.alias_name, self.channel_id, self.role_id)

    @classmethod
    async def fetch_command_aliases_by_channel(self, channel: discord.abc.GuildChannel):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT guild_id, alias_type, alias_name, channel_id, role_id FROM command_aliases WHERE guild_id=$1 AND channel_id=$2',
                channel.guild.id, channel.id
            )
            if not rows:
                return None
            aliases = []
            for row in rows:
                aliases.append(Alias(guild_id=channel.guild.id, channel_id=channel.id, alias_type=row['alias_type'], alias_name=row['alias_name'], role_id=row['role_id']))
            return aliases
    
    @classmethod
    async def fetch_command_aliases_by_channel_id(self, guild_id: Optional[int], channel_id: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT guild_id, alias_type, alias_name, channel_id, role_id FROM command_aliases WHERE guild_id=$1 AND channel_id=$2',
                guild_id, channel_id
            )
            if not rows:
                return None
            aliases = []
            for row in rows:
                aliases.append(Alias(guild_id=guild_id, channel_id=channel_id, alias_type=row['alias_type'], alias_name=row['alias_name'], role_id=row['role_id']))
            return aliases
    
    @classmethod
    async def fetch_command_alias_by_guild_and_alias_name(self, guild: discord.Guild, alias_name: str):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT guild_id, alias_type, alias_name, channel_id, role_id FROM command_aliases WHERE guild_id=$1 AND alias_name=$2',
                guild.id, alias_name
            )
            if not row:
                return None
            channel = guild.get_channel(row['channel_id'])
            return Alias(guild_id=guild.id, channel_id=channel.id, alias_type=row['alias_type'], alias_name=row['alias_name'], role_id=row['role_id'])
            
    @classmethod
    async def fetch_command_aliases_by_role(self, guild: discord.Guild, role: discord.Role):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT guild_id, alias_type, alias_name, channel_id, role_id FROM command_aliases WHERE guild_id=$1 AND role_id=$2',
                guild.id, role.id
            )
            if not rows:
                return None
            aliases = []
            for row in rows:
                channel = guild.get_channel(row['channel_id'])
                aliases.append(Alias(guild_id=guild.id, channel_id=channel.id, alias_type=row['alias_type'], alias_name=row['alias_name'], role_id=row['role_id']))
            return aliases
            
    @classmethod
    async def delete_command_alias_by_guild_and_alias_name(self, guild: discord.Guild, alias_name: str):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                'DELETE FROM command_aliases WHERE guild_id=$1 AND alias_name=$2',
                guild.id, alias_name
            )

    @classmethod
    async def delete_all_command_aliases_by_channel(self, channel: discord.abc.GuildChannel):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                'DELETE FROM command_aliases WHERE guild_id=$1 AND channel_id=$2',
                channel.guild.id, channel.id
            )
    
    @classmethod
    async def fetch_command_aliases_by_guild(self, guild: discord.Guild):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT guild_id, alias_type, alias_name, channel_id, role_id FROM command_aliases WHERE guild_id=$1 ORDER BY alias_name',
                guild.id
            )
            if not rows:
                return None
            aliases = []
            for row in rows:
                channel = guild.get_channel(row['channel_id'])
                aliases.append(Alias(guild_id=guild.id, channel_id=channel.id, alias_type=row['alias_type'], alias_name=row['alias_name'], role_id=row['role_id']))
            return aliases
    
    @classmethod
    async def fetch_all_guilds_with_command_aliases(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            guilds = {}
            for guild in bot.guilds:
                if not guild:
                    return None
                rows = await conn.fetch(
                    'SELECT guild_id, alias_type, alias_name, channel_id, role_id FROM command_aliases WHERE guild_id=$1',
                    guild.id
                )
                aliases = []
                for row in rows:
                    channel = guild.get_channel(row['channel_id'])
                    aliases.append(Alias(guild_id=guild.id, channel_id=channel.id, alias_type=row['alias_type'], alias_name=row['alias_name'], role_id=row['role_id']))
                guilds[guild] = aliases
            if not guilds:
                return None
            return guilds
            
    async def update_command_aliases_with_channel(self, channel: discord.abc.GuildChannel):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                'UPDATE command_aliases SET channel_id=$2 WHERE guild_id=$1 AND alias_type=$3 AND alias_name=$4',
                self.guild_id, channel.id, self.alias_type, self.alias_name
            )
