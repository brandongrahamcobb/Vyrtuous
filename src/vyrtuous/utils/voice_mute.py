''' cap.py The purpose of this program is to provide the Cap utility class.
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
from datetime import datetime
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot

class VoiceMute:

    def __init__(self, channel_snowflake: Optional[int], expires_at: Optional[datetime], guild_snowflake: Optional[int], member_snowflake: Optional[int], target: Optional[str]):
        self.channel_snowflake = channel_snowflake
        self.expires_at = expires_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.target = target
 
    @classmethod
    async def delete_by_channel_and_guild(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
              await conn.execute('''
                DELETE FROM active_voice_mutes
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)

    @classmethod
    async def delete_by_guild_and_member(self, guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
              await conn.execute('''
                DELETE FROM active_voice_mutes
                WHERE guild_snowflake=$1 AND member_snowflake=$2 
            ''', guild_snowflake, member_snowflake)

    @classmethod
    async def delete_by_channel_guild_member_and_target(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int], target: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
              await conn.execute('''
                DELETE FROM active_voice_mutes
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND member_snowflake=$3 AND target=$4
            ''', channel_snowflake, guild_snowflake, member_snowflake, target)

    @classmethod
    async def delete_by_channel_guild_and_target(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int], target: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
              await conn.execute('''
                DELETE FROM active_voice_mutes
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND target=$3
            ''', channel_snowflake, guild_snowflake, target)

    async def create(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_voice_mutes (channel_snowflake, created_at, expires_at, guild_snowflake, member_snowflake, target)
                VALUES ($1, NOW(), $2, $3, $4, $5)
                ON CONFLICT (channel_snowflake, guild_snowflake, member_snowflake)
            ''', self.channel_snowflake, self.expires_at, self.guild_snowflake, self.member_snowflake, self.target)

    @classmethod
    async def update_by_source_and_target(cls, source_channel_snowflake: Optional[int], target_channel_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE active_voice_mutes SET channel_snowflake=$2 WHERE channel_snowflake = $1
            ''', source_channel_snowflake, target_channel_snowflake)

    @classmethod
    async def fetch_by_channel_guild_member_and_target(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int], target: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT channel_snowflake, created_at, guild_snowflake, member_snowflake, reason, target, updated_at
                FROM active_voice_mutes
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND member_snowflake=$3 AND target=$4
            ''', channel_snowflake, guild_snowflake, member_snowflake, target)
        if row:
            return VoiceMute(channel_snowflake=row['channel_snowflake'], expires_at=row['expires_at'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason'], target=row['target'])

    @classmethod
    async def fetch_by_channel_guild_and_target(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], target: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, expires_at, guild_snowflake, member_snowflake, reason, target, updated_at
                FROM active_voice_mutes
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND target=$3
            ''', channel_snowflake, guild_snowflake, target)
        voice_mutes = []
        for row in rows:
            voice_mutes.append(VoiceMute(channel_snowflake=row['channel_snowflake'], expires_at=row['expires_at'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason'], target=row['target']))
        return voice_mutes

    @classmethod
    async def fetch_by_guild_and_target(cls, guild_snowflake: Optional[int], target: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, expires_at, guild_snowflake, member_snowflake, reason, target, updated_at
                FROM active_voice_mutes
                WHERE guild_snowflake=$1 AND target=$2
            ''', guild_snowflake, target)
        voice_mutes = []
        for row in rows:
            voice_mutes.append(VoiceMute(channel_snowflake=row['channel_snowflake'], expires_at=row['expires_at'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason'], target=row['target']))
        return voice_mutes

    @classmethod
    async def fetch_by_expired(cls, now: Optional[datetime]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, expires_at, guild_snowflake, member_snowflake, reason, target
                FROM active_voice_mutes
                WHERE expires_at IS NOT NULL AND expires_at <= $1
            ''', now)
        expired_voice_mutes = []
        for row in rows:
            expired_voice_mutes.append(VoiceMute(channel_snowflake=row['channel_snowflake'], expires_at=row['expires_at'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason'], target=row['target']))
        return expired_voice_mutes

    @classmethod
    async def fetch_by_guild_member_and_target(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int], target: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT channel_snowflake, guild_snowflake, member_snowflake, reason, target
                FROM active_voice_mutes
                WHERE guild_snowflake=$1 AND member_snowflake=$2 AND target=$3
            ''', guild_snowflake, member_snowflake, target)
        if row:
            return VoiceMute(channel_snowflake=row['channel_snowflake'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason'], target=row['target'])

    @property
    def target(self):
        return self._target

    @target.setter
    def target(self, target):
        if target not in ["room", "user"]:
            raise ValueError("Invalid target.")
        self._target = target