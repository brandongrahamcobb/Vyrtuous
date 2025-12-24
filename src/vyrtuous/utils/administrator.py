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

class Administrator:

    PLURAL = "Administrators"
    SINGULAR = "Administrator"

    def __init__(self, guild_snowflakes: list[int|None], member_snowflake: Optional[int], role_snowflakes: list[int|None]):
        self.bot = DiscordBot.get_instance()
        self.guild_snowflakes = guild_snowflakes
        self.member_snowflake: Optional[int] = member_snowflake
        self.role_snowflakes = role_snowflakes

    @classmethod
    async def revoke(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int], role_snowflake: Optional[int]):
        await cls.delete_by_guild_member_and_role(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake, role_snowflake=role_snowflake)

    @classmethod
    async def grant(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int], role_snowflake: Optional[int]):
        administrator = Administrator(guild_snowflakes=[guild_snowflake], member_snowflake=member_snowflake, role_snowflakes=[role_snowflake])
        await administrator.create()

    @classmethod
    async def delete_by_guild_member_and_role(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int], role_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM administrators
                WHERE guild_snowflake = $1 AND member_snowflake = $2 and role_snowflake = $3
            ''', guild_snowflake, member_snowflake, role_snowflake)
    
    async def create(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            for guild_id in self.guild_snowflakes:
                for role_id in self.role_snowflakes:
                    await conn.execute('''
                        INSERT INTO administrators (created_at, guild_snowflake, member_snowflake, role_snowflake)
                        VALUES (NOW(), $1, $2, $3)
                        ON CONFLICT DO NOTHING
                    ''', guild_id, self.member_snowflake, role_id)

    @classmethod
    async def fetch_members(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT
                    array_agg(guild_snowflake ORDER BY guild_snowflake) AS guild_snowflakes,
                    member_snowflake,
                    array_agg(role_snowflake ORDER BY role_snowflake) AS role_snowflakes
                FROM administrators
                GROUP BY member_snowflake
            ''')
        administrators = []
        for row in rows:
            administrators.append(
                Administrator(guild_snowflakes=row["guild_snowflakes"], member_snowflake=row["member_snowflake"], role_snowflakes=row["role_snowflakes"])
            )
        return administrators

    @classmethod
    async def fetch_members_by_role(cls, role_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT
                    guild_snowflake,
                    member_snowflake
                FROM administrators
                WHERE role_snowflake = $1
            ''', role_snowflake)
        administrators = []
        for row in rows:
            administrators.append(Administrator(guild_snowflakes=[row["guild_snowflake"]], member_snowflake=row["member_snowflake"], role_snowflakes=[role_snowflake]))
        return administrators

    @classmethod
    async def fetch_members_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT
                    member_snowflake,
                    role_snowflake
                FROM administrators
                WHERE guild_snowflake = $1
            ''', guild_snowflake)
        administrators = []
        for row in rows:
            administrators.append(Administrator(guild_snowflakes=[guild_snowflake], member_snowflake=row["member_snowflake"], role_snowflakes=[row["role_snowflake"]]))
        return administrators

    @classmethod
    async def fetch_member(cls, member_snowflake):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT
                    array_agg(guild_snowflake ORDER BY guild_snowflake) AS guild_snowflakes,
                    member_snowflake,
                    array_agg(role_snowflake ORDER BY role_snowflake) AS role_snowflakes
                FROM administrators
                WHERE member_snowflake = $1
                GROUP BY member_snowflake
            ''', member_snowflake)
        if not row:
            return None
        return Administrator(guild_snowflakes=row["guild_snowflakes"], member_snowflake=row["member_snowflake"], role_snowflakes=row["role_snowflakes"])

    @classmethod
    async def update_guild_and_role_for_member(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int], role_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE administrators SET guild_snowflake=$1, role_snowflake=$3 WHERE discord_snowflake=$2
            ''', guild_snowflake, member_snowflake, role_snowflake)
            