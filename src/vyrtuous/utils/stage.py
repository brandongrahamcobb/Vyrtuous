''' stage.py The purpose of this program is to provide the Stage utility class.

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
from vyrtuous.bot.discord_bot import DiscordBot
from typing import Optional

import discord
import time

class Stage:

    PLURAL = "Stages"
    SINGULAR = "Stage"

    def __init__(self, channel_snowflake: Optional[int], expires_in: Optional[datetime], guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        self.channel_snowflake = channel_snowflake
        self.expires_in = expires_in
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake

    @classmethod
    async def fetch_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], ):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, updated_at
                FROM active_stages
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
            if row:
                return Stage(channel_snowflake=channel_snowflake, expires_in=row['expires_in'], guild_snowflake=guild_snowflake, member_snowflake=row['member_snowflake'])

    async def send_stage_ask_to_speak_message(self, join_log: dict[int, discord.Member], member: discord.Member):
        bot = DiscordBot.get_instance()
        now = time.time()
        join_log[member.id] = [t for t in join_log[member.id] if now - t < 300]
        if len(join_log[member.id]) < 1:
            join_log[member.id].append(now)
            embed = discord.Embed(
                title=f"\U0001F399 {self.channel_snowflake} — Stage Mode",
                description=f"Ends <t:{int(self.expires_in.timestamp())}:R>",
                color=discord.Color.green()
            )
            embed.add_field(name="\u200b", value="**Ask to speak!**", inline=False)
            await bot.get_channel(self.channel_snowflake).send(embed=embed)

    @classmethod
    async def update_duration(cls, channel_snowflake: int, expires_in, guild_snowflake: int):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute(f'''
                UPDATE active_stages
                SET expires_in=$2, updated_at=NOW()
                WHERE channel_snowflake=$1 AND guild_snowflake=$3
            ''', channel_snowflake, expires_in, guild_snowflake)

    @classmethod
    async def update_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE active_stages
                SET channel_snowflake=$1
                WHERE guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)

    @classmethod
    async def fetch_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT
                    s.channel_snowflake,
                    s.member_snowflake,
                    s.expires_in,
                    COALESCE(COUNT(v.member_snowflake), 0) AS active_mutes
                FROM active_stages s
                LEFT JOIN active_voice_mutes v
                    ON s.guild_snowflake = v.guild_snowflake
                   AND s.channel_snowflake = v.channel_snowflake
                   AND v.target = 'room'
                   AND (v.expires_in IS NULL OR v.expires_in > NOW())
                WHERE s.guild_snowflake=$1
                GROUP BY
                    s.channel_snowflake,
                    s.member_snowflake,
                    s.expires_in
                ORDER BY s.channel_snowflake
            ''', guild_snowflake)
        stages = []
        if rows:
            for row in rows:
                stages.append(Stage(channel_snowflake=row['channel_snowflake'], expires_in=row['expires_in'], guild_snowflake=guild_snowflake, member_snowflake=row['member_snowflake']))

    @classmethod
    async def update_by_source_and_target(cls, source_channel_snowflake: Optional[int], target_channel_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE active_stages
                SET channel_snowflake=$2
                WHERE channel_snowflake=$1
            ''', source_channel_snowflake, target_channel_snowflake)

    @classmethod
    async def fetch_by_guild_and_channel(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT
                    s.channel_snowflake,
                    s.member_snowflake,
                    s.expires_in,
                    COALESCE(COUNT(v.member_snowflake), 0) AS active_mutes
                FROM active_stages s
                LEFT JOIN active_voice_mutes v
                    ON s.guild_snowflake = v.guild_snowflake
                   AND s.channel_snowflake = v.channel_snowflake
                   AND v.target = 'room'
                   AND (v.expires_in IS NULL OR v.expires_in > NOW())
                WHERE s.channel_snowflake=$1
                  AND s.guild_snowflake=$2
                GROUP BY
                    s.channel_snowflake,
                    s.member_snowflake,
                    s.expires_in
            ''', channel_snowflake, guild_snowflake)
        if row:
            return Stage(channel_snowflake=channel_snowflake, expires_in=row['expires_in'], guild_snowflake=guild_snowflake, member_snowflake=row['member_snowflake'])

    async def create(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (channel_snowflake, expires_in, guild_snowflake, member_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT DO NOTHING
            ''', self.channel_snowflake, self.expires_in, self.guild_snowflake, self.member_snowflake)

    @classmethod
    async def fetch_by_expired(cls, now: Optional[datetime]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, updated_at
                FROM active_stages
                WHERE expires_in IS NOT NULL AND expires_in <= $1
            ''', now)
        expired_stages = []
        for row in rows:
            expired_stages.append(Stage(channel_snowflake=row['channel_snowflake'], expires_in=row['expires_in'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake']))
        return expired_stages
    
    @classmethod
    async def fetch_all(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, updated_at
                FROM active_stages
            ''')
        stages = []
        for row in rows:
            stages.append(Stage(channel_snowflake=row['channel_snowflake'], expires_in=row['expires_in'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake']))
        return stages

    @classmethod
    async def delete_by_channel_and_guild(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM active_stages
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)