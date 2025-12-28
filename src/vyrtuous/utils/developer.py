
''' test_admin_helpers.py The purpose of this program is to provide the shared test variables for tests using Discord objects.
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
from vyrtuous.inc.helpers import *

class Developer:

    PLURAL = "Developers"
    SINGULAR = "Developer"

    def __init__(self, guild_snowflake: Optional[int], member_snowflake: Optional[str]):
        self.bot = DiscordBot.get_instance()
        self.guild_snowflake = guild_snowflake
        self.member_snowflake: Optional[int] = member_snowflake
        self.member_mention: Optional[str] = f"<@{member_snowflake}>"

    async def grant(self):
        async with self.bot.db_pool.acquire() as conn:   
            await conn.execute('''
                INSERT INTO developers (created_at, guild_snowflake, member_snowflake)
                VALUES (NOW(), $1, $2)
                ON CONFLICT DO NOTHING
            ''', self.guild_snowflake, self.member_snowflake)

    async def revoke(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM developers
                WHERE guild_snowflake = $1 AND member_snowflake = $2
            ''', self.guild_snowflake, self.member_snowflake)
    
    @classmethod
    async def fetch_all(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT guild_snowflake, member_snowflake
                FROM developers
                ORDER BY member_snowflake
            ''')
        developers = []
        for row in rows:
            developers.append(Developer(guild_snowflake=row["guild_snowflake"], member_snowflake=row["member_snowflake"]))
        return developers

    @classmethod
    async def fetch_guilds_by_member(cls, member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT guild_snowflake FROM developers WHERE member_snowflake=$1
            ''', member_snowflake)
        guild_snowflakes = []
        for row in rows:
            guild_snowflakes.append(row['guild_snowflake'])
        return guild_snowflakes

    @classmethod
    async def fetch_members_by_guild(cls, guild_snowflake):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT guild_snowflake, member_snowflake
                FROM developers
                WHERE guild_snowflake = $1
            ''', guild_snowflake)
        developers = []
        for row in rows:
            developers.append(Developer(guild_snowflake=row["guild_snowflake"], member_snowflake=row["member_snowflake"]))
        return developers
    
    @classmethod
    async def delete_by_guild_and_member(cls, guild_snowflake: int, member_snowflake: int):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM developers
                WHERE guild_snowflake = $1 AND member_snowflake = $2
            ''', guild_snowflake, member_snowflake)
