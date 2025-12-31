''' administrator.py The purpose of this program is to inherit from the user class as an administrator.

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

class Administrator:

    PLURAL = "Administrators"
    SINGULAR = "Administrator"

    def __init__(self, guild_snowflake: list[int|None], member_snowflake: Optional[int], role_snowflake: list[int|None]):
        self.bot = DiscordBot.get_instance()
        self.guild_snowflake = guild_snowflake
        self.member_snowflake: Optional[int] = member_snowflake
        self.member_mention: Optional[str] = f"<@{member_snowflake}>"
        self.role_snowflake = role_snowflake
    
    async def grant(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO administrators (created_at, guild_snowflake, member_snowflake, role_snowflake)
                VALUES (NOW(), $1, $2, $3)
                ON CONFLICT DO NOTHING
            ''', self.guild_snowflake, self.member_snowflake, self.role_snowflake)

    async def revoke(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM administrators
                WHERE guild_snowflake = $1 AND member_snowflake = $2 and role_snowflake = $3
            ''', self.guild_snowflake, self.member_snowflake, self.role_snowflake)

    @classmethod
    async def fetch_all(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT guild_snowflake, member_snowflake, role_snowflake
                FROM administrators
            ''')
        administrators = []
        for row in rows:
            administrators.append(
                Administrator(guild_snowflake=row["guild_snowflake"], member_snowflake=row["member_snowflake"], role_snowflake=row["role_snowflake"])
            )
        return administrators

    @classmethod
    async def fetch_member(cls, member_snowflake):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT guild_snowflake, member_snowflake, role_snowflake
                FROM administrators
                WHERE member_snowflake = $1
            ''', member_snowflake)
        if not row:
            return None
        return Administrator(guild_snowflake=row["guild_snowflake"], member_snowflake=row["member_snowflake"], role_snowflake=row["role_snowflake"])

    @classmethod
    async def update_guild_and_role_for_member(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int], role_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE administrators SET guild_snowflake=$1, role_snowflake=$3 WHERE member_snowflake=$2
            ''', guild_snowflake, member_snowflake, role_snowflake)

    @classmethod
    async def fetch_by_guild_and_member(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT guild_snowflake, member_snowflake, role_snowflake FROM administrators WHERE guild_snowflake=$1 AND member_snowflake=$2
            ''', guild_snowflake, member_snowflake)
        administrators = []
        if rows:
            for row in rows:
                administrators.append(Administrator(guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], role_snowflake=row['row_snowflake']))
        return administrators
    
    @classmethod
    async def fetch_by_guild_and_role(cls, guild_snowflake: Optional[int], role_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT guild_snowflake, member_snowflake, role_snowflake FROM administrators WHERE guild_snowflake=$1 AND role_snowflake=$2
            ''', guild_snowflake, role_snowflake)
        administrators = []
        if rows:
            for row in rows:
                administrators.append(Administrator(guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], role_snowflake=row['row_snowflake']))
        return administrators
    
    @classmethod
    async def fetch_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT guild_snowflake, member_snowflake, role_snowflake FROM administrators WHERE guild_snowflake=$1
            ''', guild_snowflake)
        administrators = []
        if rows:
            for row in rows:
                administrators.append(Administrator(guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], role_snowflake=row['row_snowflake']))
        return administrators
            