''' developer.py The purpose of this program is to inherit from the user class as a developer.

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
from datetime import datetime, timezone
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
import uuid

class DeveloperLog:

    def __init__(self, channel_snowflake: Optional[int], developer_snowflakes: list[int], guild_snowflake: Optional[int], message_snowflake: Optional[int], created_at: Optional[datetime] = None, id: Optional[str] = None, notes: Optional[str] = None, resolved: bool = False, updated_at: Optional[datetime] = None):
        self.bot = DiscordBot.get_instance()
        self.channel_snowflake = channel_snowflake
        self.created_at: datetime = created_at or datetime.now(timezone.utc)
        self.developer_snowflakes = developer_snowflakes
        self.id = id
        self.guild_snowflake = guild_snowflake
        self.message_snowflake = message_snowflake
        self.notes = notes
        self.resolved = resolved
        self.updated_at: datetime = updated_at or datetime.now(timezone.utc)

    async def create(self):
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT id
                FROM developer_logs
            ''')
        real_uuid = str(uuid.uuid4())
        while any(real_uuid == row['id'] for row in rows):
            real_uuid = str(uuid.uuid4())
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO developer_logs (channel_snowflake, created_at, developer_snowflakes, guild_snowflake, id, message_snowflake, notes, resolved, updated_at)
                VALUES ($1, NOW(), $2, $3, $4, $5, $6, $7, NOW())
                ON CONFLICT DO NOTHING
            ''', self.channel_snowflake, self.developer_snowflakes, self.guild_snowflake, real_uuid, self.message_snowflake, self.notes, self.resolved)
            self.created_at = datetime.now(timezone.utc)
            self.updated_at = datetime.now(timezone.utc)
        return DeveloperLog(channel_snowflake=self.channel_snowflake, created_at=self.created_at, developer_snowflakes=self.developer_snowflakes, guild_snowflake=self.guild_snowflake, id=real_uuid, message_snowflake=self.message_snowflake, notes=self.notes, resolved=self.resolved, updated_at=self.updated_at)

    async def resolve(self):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE developer_logs
                SET resolved=TRUE, updated_at=NOW()
                WHERE id=$1;
            ''', self.id)
            self.resolved = True
            self.updated_at = datetime.now(timezone.utc)
    
    async def append(self, notes: str):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE developer_logs
                SET notes = COALESCE(notes, '') || $1,
                    updated_at=NOW()
                WHERE id=$2
            ''', f"\n{notes}", self.id)
            self.notes = self.notes + f'{notes}'
            self.updated_at = datetime.now(timezone.utc)
    
    async def overwrite(self, notes: str):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE developer_logs
                SET notes=$1, updated_at=NOW()
                WHERE id=$2
            ''', notes, self.id)
            self.notes = notes
            self.updated_at = datetime.now(timezone.utc)

    @classmethod
    async def fetch_all(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, developer_snowflakes,  guild_snowflake, id, message_snowflake, notes, resolved, updated_at
                FROM developer_logs
                ORDER BY created_at
            ''')
        developer_logs = []
        if rows:
            for row in rows:
                developer_logs.append(DeveloperLog(channel_snowflake=row['channel_snowflake'], created_at=row['created_at'], developer_snowflakes=row['developer_snowflakes'], guild_snowflake=row["guild_snowflake"], id=row['id'], message_snowflake=row["message_snowflake"], notes=row['notes'], resolved=row["resolved"], updated_at=row['updated_at']))
        return developer_logs

    @classmethod
    async def fetch_unresolved_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, developer_snowflakes,  guild_snowflake, id, message_snowflake, notes, resolved, updated_at
                FROM developer_logs
                WHERE guild_snowflake=$1 AND resolved=FALSE
            ''', guild_snowflake)
        developer_logs = []
        if rows:
            for row in rows:
                developer_logs.append(DeveloperLog(channel_snowflake=row['channel_snowflake'], created_at=row['created_at'], developer_snowflakes=row['developer_snowflakes'], guild_snowflake=guild_snowflake, id=row['id'], message_snowflake=row["message_snowflake"], notes=row['notes'], resolved=row["resolved"], updated_at=row['updated_at']))
        return developer_logs
    
    @classmethod
    async def fetch_unresolved_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, developer_snowflakes,  guild_snowflake, id, message_snowflake, notes, resolved, updated_at
                FROM developer_logs
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND resolved=FALSE
            ''', channel_snowflake, guild_snowflake)
        developer_logs = []
        if rows:
            for row in rows:
                developer_logs.append(DeveloperLog(channel_snowflake=channel_snowflake, created_at=row['created_at'], developer_snowflakes=row['developer_snowflakes'], guild_snowflake=guild_snowflake, id=row['id'], message_snowflake=row["message_snowflake"], notes=row['notes'], resolved=row["resolved"], updated_at=row['updated_at']))
        return developer_logs
    
    @classmethod
    async def fetch_all_unresolved(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, developer_snowflakes,  guild_snowflake, id, message_snowflake, notes, resolved, updated_at
                FROM developer_logs
                WHERE resolved=FALSE
            ''')
        developer_logs = []
        if rows:
            for row in rows:
                developer_logs.append(DeveloperLog(channel_snowflake=row['channel_snowflake'], created_at=row['created_at'], developer_snowflakes=row['developer_snowflakes'], guild_snowflake=row["guild_snowflake"], id=row['id'], message_snowflake=row["message_snowflake"], notes=row['notes'], resolved=False, updated_at=row['updated_at']))
        return developer_logs
    
    @classmethod
    async def delete_by_resolved(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM developer_logs
                WHERE resolved=TRUE
            ''')
    
    @classmethod
    async def fetch_unresolved_by_reference(cls, id: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT channel_snowflake, created_at, developer_snowflakes,  guild_snowflake, id, message_snowflake, notes, resolved, updated_at
                FROM developer_logs
                WHERE id=$1 AND resolved=FALSE
            ''', id)
        issue = None
        if row:
            issue = DeveloperLog(channel_snowflake=row['channel_snowflake'], created_at=row['created_at'], developer_snowflakes=row['developer_snowflakes'], guild_snowflake=row["guild_snowflake"], id=id, message_snowflake=row["message_snowflake"], notes=row['notes'], resolved=row['resolved'], updated_at=row['updated_at'])
        return issue

    @classmethod
    async def fetch_resolved_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, developer_snowflakes,  guild_snowflake, id, message_snowflake, notes, resolved, updated_at
                FROM developer_logs
                WHERE guild_snowflake=$1 AND resolved=TRUE
            ''', guild_snowflake)
        developer_logs = []
        if rows:
            for row in rows:
                developer_logs.append(DeveloperLog(channel_snowflake=row['channel_snowflake'], created_at=row['created_at'], developer_snowflakes=row['developer_snowflakes'], guild_snowflake=guild_snowflake, id=row['id'], message_snowflake=row["message_snowflake"], notes=row['notes'], resolved=row["resolved"], updated_at=row['updated_at']))
        return developer_logs
    
    @classmethod
    async def fetch_resolved_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, developer_snowflakes,  guild_snowflake, id, message_snowflake, notes, resolved, updated_at
                FROM developer_logs
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND resolved=TRUE
            ''', channel_snowflake, guild_snowflake)
        developer_logs = []
        if rows:
            for row in rows:
                developer_logs.append(DeveloperLog(channel_snowflake=channel_snowflake, created_at=row['created_at'], developer_snowflakes=row['developer_snowflakes'], guild_snowflake=guild_snowflake, id=row['id'], message_snowflake=row["message_snowflake"], notes=row['notes'], resolved=row["resolved"], updated_at=row['updated_at']))
        return developer_logs
    
    @classmethod
    async def fetch_all_resolved(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, developer_snowflakes,  guild_snowflake, id, message_snowflake, notes, resolved, updated_at
                FROM developer_logs
                WHERE resolved=TRUE
            ''')
        developer_logs = []
        if rows:
            for row in rows:
                developer_logs.append(DeveloperLog(channel_snowflake=row['channel_snowflake'], created_at=row['created_at'], developer_snowflakes=row['developer_snowflakes'], guild_snowflake=row["guild_snowflake"], id=row['id'], message_snowflake=row["message_snowflake"], notes=row['notes'], resolved=True, updated_at=row['updated_at']))
        return developer_logs
    
    @classmethod
    async def fetch_resolved_by_reference(cls, id: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT channel_snowflake, created_at, developer_snowflakes,  guild_snowflake, id, message_snowflake, notes, resolved, updated_at
                FROM developer_logs
                WHERE id=$1 AND resolved=TRUE
            ''', id)
        issue = None
        if row:
            issue = DeveloperLog(channel_snowflake=row['channel_snowflake'], created_at=row['created_at'], developer_snowflakes=row['developer_snowflakes'], guild_snowflake=row["guild_snowflake"], id=id, message_snowflake=row["message_snowflake"], notes=row['notes'], resolved=row['resolved'], updated_at=row['updated_at'])
        return issue
    
    async def assign(self, member_snowflake: Optional[int]):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE developer_logs
                SET developer_snowflakes = array_append(COALESCE(developer_snowflakes, '{}'), $1),
                    updated_at = NOW()
                WHERE id = $2 AND NOT $1 = ANY(COALESCE(developer_snowflakes, '{}'))
            ''', member_snowflake, self.id)
            self.updated_at = datetime.now(timezone.utc)

    async def unassign(self, member_snowflake: Optional[int]):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE developer_logs
                SET developer_snowflakes = array_remove(COALESCE(developer_snowflakes, '{}'), $1),
                    updated_at = NOW()
                WHERE id = $2
            ''', member_snowflake, self.id)
            self.updated_at = datetime.now(timezone.utc)