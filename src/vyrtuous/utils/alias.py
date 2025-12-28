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
from datetime import datetime, timedelta, timezone
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.check_service import *
from vyrtuous.service.member_service import MemberService
from vyrtuous.utils.reason import Reason
from vyrtuous.utils.ban import Ban
from vyrtuous.utils.flag import Flag
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.voice_mute import VoiceMute

import discord

class Alias:
    
    MODERATION_TABLES = {
        "Ban": "active_bans",
        "Flag": "active_flags",
        "TextMute": "active_text_mutes",
        "VoiceMute": "active_voice_mutes",
    }

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
        self.channel_service = ChannelService()
        self.member_service = MemberService()
        logger.info("test")
        
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
        logger.info("test")
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
    async def update_reason(cls, channel_snowflake: int, guild_snowflake: int, member_snowflake: int, moderation_type, updated_reason: str):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute(f'''
                UPDATE {cls.get_table_name_by_moderation_type(moderation_type)}
                SET reason=$4, updated_at=NOW()
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND member_snowflake=$3
            ''', channel_snowflake, guild_snowflake, member_snowflake, updated_reason)

    @classmethod
    async def update_duration(cls, channel_snowflake: int, expires_at, guild_snowflake: int, member_snowflake: int, moderation_type):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute(f'''
                UPDATE {cls.get_table_name_by_moderation_type(moderation_type)}
                SET expires_at=$2, updated_at=NOW()
                WHERE channel_snowflake=$1 AND guild_snowflake=$3 AND member_snowflake=$4
            ''', channel_snowflake, expires_at, guild_snowflake, member_snowflake)
            
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
        
    @classmethod
    async def check_channel_member_and_permission(cls, alias, member_string, message):
        try:
            channel_obj = await cls.channel_service.resolve_channel(message, alias.channel_snowflake)
            member_obj = await cls.member_service.resolve_member(message, member_string)
            await has_equal_or_higher_role(message, channel_snowflake=channel_obj.id, guild_snowflake=alias.guild_snowflake, ember_snowflake=member_obj.id, sender_snowflake=message.author.id)
            if member_obj.id == message.guild.me.id:
                raise Exception("\U000026A0\U0000FE0F You cannot {alias.alias_type} {message.guild.me.mention}.")
        except:
            raise

    @classmethod
    def get_table_name_by_moderation_type(cls, moderation_type):
        try:
            return cls.MODERATION_TABLES[type(moderation_type).__name__]
        except KeyError:
            raise ValueError(f"Unknown moderation type: {moderation_type}")

    @classmethod 
    async def get_existing_moderation(cls, alias, channel_snowflake, guild_snowflake, member_snowflake):
        match alias.alias_type:
            case "ban":
                return await Ban.fetch_by_channel_guild_and_member(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
            case "flag":
                return await Flag.fetch_by_channel_guild_and_member(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
            case "tmute":
                return await TextMute.fetch_by_channel_guild_and_member(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
            case "voice_mute":
                return await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake, target="user")
