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

class Moderator:

    PLURAL = "Moderators"
    SINGULAR = "Moderator"

    def __init__(self, channel_snowflake: Optional[str], guild_snowflake: Optional[int], member_snowflake: Optional[str]):
        self.bot = DiscordBot.get_instance()
        self.channel_snowflake: Optional[int] = channel_snowflake
        self.guild_snowflake = guild_snowflake
        self.member_snowflake: Optional[int] = member_snowflake

    @classmethod
    async def update_by_source_and_target(cls, source_channel_snowflake: Optional[int], target_channel_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE moderators SET channel_snowflake=$2 WHERE channel_snowflake=$1
            ''', source_channel_snowflake, target_channel_snowflake)

    @classmethod
    async def delete_by_channel_and_member(cls, channel_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM moderators WHERE channel_snowflake=$1 AND member_snowflake=$2
            ''', channel_snowflake, member_snowflake)

    async def create(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderators (channel_snowflake, created_at, guild_snowflake, member_snowflake)
                VALUES ($1, NOW(), $2, $3)
                ON CONFLICT DO NOTHING
            ''', self.channel_snowflake, self.guild_snowflake, self.member_snowflake)

    @classmethod
    async def fetch_channels_by_guild_and_member(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake
                FROM moderators
                WHERE guild_snowflake=$1 AND member_snowflake=$2
            ''', guild_snowflake, member_snowflake)
            return [row['channel_snowflake'] for row in rows]

    @classmethod
    async def delete_channel(cls, channel_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM moderators WHERE channel_snowflake=$1
            ''', channel_snowflake)

    @classmethod
    async def fetch_members_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT member_snowflake
                FROM moderators
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
        moderators = []
        if rows:
            for row in rows:
                moderators.append(Moderator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=row["member_snowflake"]))

    @classmethod
    async def fetch_members_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, member_snowflake
                FROM moderators
                WHERE guild_snowflake=$1
            ''', guild_snowflake)
            moderators = []
            for row in rows:
                moderators.append(
                    Moderator(channel_snowflake=row['channel_snowflake'], guild_snowflake=guild_snowflake, member_snowflake=row['member_snowflake'])
                )
            return moderators

    @classmethod
    async def revoke(cls, channel_snowflake: Optional[int], member_snowflake: Optional[int]):
        await cls.delete_by_channel_and_member(channel_snowflake=channel_snowflake, member_snowflake=member_snowflake)

    @classmethod
    async def grant(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        moderator = Moderator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
        await moderator.create()