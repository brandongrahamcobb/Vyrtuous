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

class Coordinator:

    PLURAL = "Coordinators"
    SINGULAR = "Coordinator"

    def __init__(self, channel_id: Optional[str], guild_id: Optional[int], member_id: Optional[str]):
        self.bot = DiscordBot.get_instance()
        self.channel_id: Optional[int] = channel_id
        self.guild_id = guild_id
        self.member_id: Optional[int] = member_id

    @classmethod
    async def update_source_channel_id_to_target_channel_id(cls, source_channel_id: Optional[int], target_channel_id: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE coordinators SET channel_snowflake=$2 WHERE channel_snowflake=$1
            ''', source_channel_id, target_channel_id)

    @classmethod
    async def delete_channel_id_for_member_id(cls, channel_id: Optional[int], member_id: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM coordinators WHERE channel_snowflake=$1 AND member_snowflake=$2
            ''', channel_id, member_id)

    @classmethod
    async def delete_channel_ids_for_channel_id(cls, channel_id: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM coordinators WHERE channel_snowflake=$1
            ''', channel_id)

    async def set_channel_id_for_member(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO coordinators (channel_snowflake, guild_snowflake, member_snowflake)
                VALUES ($1,$2,$3)
                ON CONFLICT DO NOTHING
            ''', self.channel_id, self.guild_id, self.member_id)

    @classmethod
    async def fetch_channel_ids_for_guild_id_and_member_id(cls, guild_id: Optional[int], member_id: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake
                FROM coordinators
                WHERE guild_snowflake=$1 AND member_snowflake=$2
            ''', guild_id, member_id)
            return [row['channel_snowflake'] for row in rows]

    @classmethod
    async def fetch_discord_snowflakes_for_channel_id(cls, guild_id: Optional[int], channel_id: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT member_snowflake
                FROM coordinators
                WHERE guild_snowflake=$1 AND channel_snowflake=$2
            ''', guild_id, channel_id)
            return [row['member_snowflake'] for row in rows]

    @classmethod
    async def fetch_all_members_in_guild(cls, guild_id: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, member_snowflake
                FROM coordinators
                WHERE guild_snowflake=$1
            ''', guild_id)
            coordinators = []
            for row in rows:
                moderators.append(
                    Coordinator(channel_id=row['channel_snowflake'],guild_id=guild_id,member_id=row['member_snowflake'])
                )
            return coordinators

    @classmethod
    async def revoke(cls, channel_id: Optional[int], member_id: Optional[int]):
        await cls.delete_channel_id_for_member_id(channel_id=channel_id, member_id=member_id)

    @classmethod
    async def grant(cls, channel_id: Optional[int], guild_id: Optional[int], member_id: Optional[int]):
        coordinator = Coordinator(channel_id=channel_id, guild_id=guild_id, member_id=member_id)
        await coordinator.set_channel_id_for_member()