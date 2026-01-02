''' coordinator.py The purpose of this program is to inherit from the user class as a coordinator.

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

class Coordinator:

    PLURAL = "Coordinators"
    SINGULAR = "Coordinator"

    def __init__(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        self.bot = DiscordBot.get_instance()
        self.channel_snowflake = channel_snowflake
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None

    @classmethod
    async def update_by_source_and_target(cls, source_channel_snowflake: Optional[int], target_channel_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE coordinators
                SET channel_snowflake=$2
                WHERE channel_snowflake=$1
            ''', source_channel_snowflake, target_channel_snowflake)

    @classmethod
    async def delete_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM coordinators
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)

    @classmethod
    async def delete_by_channel_guild_and_member(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM coordinators
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND member_snowflake=$3
            ''', channel_snowflake, guild_snowflake, member_snowflake)

    async def grant(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO coordinators (channel_snowflake, created_at, guild_snowflake, member_snowflake)
                VALUES ($1, NOW(), $2, $3)
                ON CONFLICT DO NOTHING
            ''', self.channel_snowflake, self.guild_snowflake, self.member_snowflake)

    @classmethod
    async def fetch_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT created_at, channel_snowflake, guild_snowflake, member_snowflake, updated_at
                FROM coordinators
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
        coordinators = []
        if rows:
            for row in rows:
                coordinators.append(Coordinator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=row["member_snowflake"]))
        return coordinators

    @classmethod
    async def fetch_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT created_at, channel_snowflake, guild_snowflake, member_snowflake, updated_at
                FROM coordinators
                WHERE guild_snowflake=$1
            ''', guild_snowflake)
            coordinators = []
            for row in rows:
                coordinators.append(
                    Coordinator(channel_snowflake=row['channel_snowflake'], guild_snowflake=guild_snowflake, member_snowflake=row['member_snowflake'])
                )
            return coordinators

    @classmethod
    async def fetch_by_channel_guild_and_member(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT created_at, channel_snowflake, guild_snowflake, member_snowflake, updated_at
                FROM coordinators
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND member_snowflake=$3
            ''', channel_snowflake, guild_snowflake, member_snowflake)
        coordinator = None
        if row:
            coordinator = Coordinator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
        return coordinator
    
    @classmethod
    async def fetch_by_guild_and_member(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT created_at, channel_snowflake, guild_snowflake, member_snowflake, updated_at
                FROM coordinators
                WHERE guild_snowflake=$1 AND member_snowflake=$2
            ''', guild_snowflake, member_snowflake)
        coordinators = None
        if rows:
            for row in rows:
                coordinators.append(Coordinator(channel_snowflake=row['channel_snowflake'], guild_snowflake=guild_snowflake, member_snowflake=member_snowflake))
        return coordinators
    
    async def revoke(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM coordinators
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND member_snowflake=$3
            ''', self.channel_snowflake, self.guild_snowflake, self.member_snowflake)

    @classmethod
    async def fetch_all(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT created_at, channel_snowflake, guild_snowflake, member_snowflake, updated_at
                FROM coordinators
            ''')
        coordinators = []
        if rows:
            for row in rows:
                coordinators.append(Coordinator(channel_snowflake=row['channel_snowflake'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake']))
        return coordinators