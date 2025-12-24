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
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot

class Cap:

    def __init__(self, channel_snowflake: Optional[int], duration: Optional[int], guild_snowflake: Optional[int], moderation_type: Optional[str]):
        self.channel_snowflake = channel_snowflake
        self.duration_seconds = duration
        self.guild_snowflake = guild_snowflake
        self.moderation_type = moderation_type
 
    @classmethod
    async def fetch_by_channel_and_guild(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int]) -> list[tuple[str, str]]:
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT duration_seconds, moderation_type FROM active_caps WHERE guild_snowflake=$1 AND channel_snowflake=$2',
                guild_snowflake, channel_snowflake
            )
            return [(r['duration_seconds'], r['moderation_type']) for r in rows]

    @classmethod
    async def delete_by_channel_guild_and_moderation_type(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int], moderation_type: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
              await conn.execute('''
                DELETE FROM active_caps
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND moderation_type=$3
            ''', channel_snowflake, guild_snowflake, moderation_type)

    @classmethod
    async def fetch_by_channel_guild_and_moderation_type(cls, guild_snowflake: Optional[int], channel_snowflake: Optional[int], moderation_type: Optional[str]) -> list[tuple[str, str]]:
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT duration_seconds, moderation_type
                FROM active_caps
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND moderation_type=$3
            ''', channel_snowflake, guild_snowflake, moderation_type)
        if row:
            return (row['duration_seconds'], row['moderation_type'])

    async def grant(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_caps (channel_snowflake, duration_seconds, guild_snowflake, moderation_type)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (channel_snowflake, guild_snowflake, moderation_type)
                DO UPDATE SET duration_seconds = EXCLUDED.duration_seconds
            ''', self.channel_snowflake, self.duration_seconds, self.guild_snowflake, self.moderation_type)

    @classmethod
    async def update_by_source_and_target(cls, source_channel_snowflake: Optional[int], target_channel_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE active_caps SET channel_snowflake=$2 WHERE channel_snowflake = $1
            ''', source_channel_snowflake, target_channel_snowflake)

    @classmethod
    async def fetch_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, duration_seconds, guild_snowflake, updated_at
                FROM active_caps
                WHERE guild_snowflake=$1
            ''', guild_snowflake)
        caps = []
        for row in rows:
            caps.append(Cap(channel_snowflake=row['channel_snowflake'], duration=row['duration_seconds'], guild_snowflake=row['guild_snowflake']))
        return caps
    
    @classmethod
    async def update_by_channel_and_duration(cls, channel_snowflake: Optional[int], duration: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE active_caps SET duration_seconds=$2 WHERE channel_snowflake=$1
            ''', channel_snowflake, duration)

    @property
    def moderation_type(self):
        return self._moderation_type
    
    @moderation_type.setter
    def moderation_type(self, moderation_type: Optional[str]):
        if moderation_type not in ('ban', 'mute', 'tmute'):
            raise ValueError("Invalid moderation_type.")
        self._moderation_type = moderation_type