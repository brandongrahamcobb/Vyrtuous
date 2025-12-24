''' alias.py The purpose of this program is to provide alias utilities.
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
from collections import defaultdict
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot

import discord

class Alias:
    
    def __init__(self, alias_name: Optional[str], alias_type: Optional[str], channel_snowflake: Optional[int], guild_snowflake: Optional[int], role_snowflake: Optional[int]):
        self.alias_type = alias_type
        self.alias_name = alias_name
        self.bot = DiscordBot.get_instance()
        self.channel_snowflake = channel_snowflake
        self.guild_snowflake = guild_snowflake
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
        self.role_snowflake = role_snowflake
        
    @classmethod
    def format_aliases(self, aliases) -> list[str]:
        if not aliases:
            return []
        grouped = defaultdict(list)
        lines = []
        for alias in aliases:
            grouped[(alias.channel_snowflake, alias.alias_type)].append(alias)
        for (alias.channel_snowflake, alias.alias_type), channel_aliases in grouped.items():
            lines.append(f'**{alias.alias_type.capitalize()}**')
            for alias in channel_aliases:
                if alias.role_snowflake:
                    lines.append(f'`{alias.alias_name}` → <@&{alias.role_snowflake}>')
                else:
                    lines.append(f'`{alias.alias_name}`')
        return lines

    async def create(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO command_aliases (alias_name, alias_type, channel_snowflake, created_at, guild_snowflake, role_snowflake)
                VALUES ($1, $2, $3, NOW(), $4, $5)
            ''', self.alias_name, self.alias_type, self.channel_snowflake, self.guild_snowflake, self.role_snowflake)

    @classmethod
    async def fetch_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT alias_name, alias_type, channel_snowflake, created_at, guild_snowflake, role_snowflake, updated_at
                FROM command_aliases WHERE channel_snowflake=$1 AND guild_snowflake=$2
                ORDER BY alias_name
            ''', channel_snowflake, guild_snowflake)
        if rows:
            aliases = []
            for row in rows:
                aliases.append(Alias(alias_type=row['alias_type'], alias_name=row['alias_name'], channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, role_snowflake=row['role_snowflake']))
            return aliases
    
    @classmethod
    async def fetch_by_guild_and_name(cls, alias_name: Optional[str], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT alias_name, alias_type, channel_snowflake, created_at, guild_snowflake, role_snowflake, updated_at
                FROM command_aliases WHERE alias_name=$1 AND guild_snowflake=$2
                ORDER BY alias_name
            ''', alias_name, guild_snowflake)
        if row:
            return Alias(alias_name=alias_name, alias_type=row['alias_type'], channel_snowflake=row['channel_snowflake'], guild_snowflake=guild_snowflake, role_snowflake=row['role_snowflake'])
            
    @classmethod
    async def fetch_by_guild_and_role(cls, guild_snowflake: Optional[int], role_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT alias_name, alias_type, channel_snowflake, created_at, guild_snowflake, role_snowflake, updated_at
                FROM command_aliases WHERE guild_snowflake=$1 AND role_snowflake=$2
                ORDER BY alias_name
            ''', guild_snowflake, role_snowflake)
        if rows:
            aliases = []
            for row in rows:
                aliases.append(Alias(alias_name=row['alias_name'], alias_type=row['alias_type'], channel_snowflake=row['channel_snowflake'], guild_snowflake=guild_snowflake, role_snowflake=role_snowflake))
            return aliases
            
    @classmethod
    async def delete_by_guild_and_name(cls, alias_name: Optional[str], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM command_aliases
                WHERE alias_name=$1 AND guild_snowflake=$2
            ''', alias_name, guild_snowflake)

    @classmethod
    async def delete_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM command_aliases WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
    
    @classmethod
    async def fetch_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT alias_name, alias_type, channel_snowflake, created_at, guild_snowflake, role_snowflake, updated_at
                FROM command_aliases WHERE guild_snowflake=$1
                ORDER BY alias_name
            ''', guild_snowflake)
        if rows:
            aliases = []
            for row in rows:
                aliases.append(Alias(alias_name=row['alias_name'], alias_type=row['alias_type'], channel_snowflake=row['channel_snowflake'], guild_snowflake=guild_snowflake, role_snowflake=row['role_snowflake']))
            return aliases
    
    @classmethod
    async def fetch_all(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            for guild in bot.guilds:
                rows = await conn.fetch('''
                    SELECT alias_name, alias_type, channel_snowflake, created_at, guild_snowflake, role_snowflake, updated_at
                    FROM command_aliases WHERE guild_snowflake=$1
                ''', guild.id)
        if rows:
            aliases = []
            for row in rows:
                aliases.append(Alias(alias_name=row['alias_name'], alias_type=row['alias_type'], channel_snowflake=row['channel_snowflake'], guild_snowflake=row['guild_snowflake'], role_snowflake=row['role_snowflake']))
            return aliases
            
    @classmethod
    async def update_by_source_and_target(cls, source_channel_snowflake: Optional[int], target_channel_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE command_aliases SET channel_snowflake=$2
                WHERE channel_snowflake=$1
            ''', source_channel_snowflake, target_channel_snowflake)

    @property
    def alias_type(self):
        return self._alias_type
    
    @alias_type.setter
    def alias_type(self, alias_type: Optional[str]):
        if alias_type not in ('cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute', 'role', 'unrole'):
            raise ValueError("Invalid alias_type.")
        self._alias_type = alias_type
