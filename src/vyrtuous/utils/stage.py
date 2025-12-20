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

    def __init__(self, stage_expires_at: datetime, stage_channel_id: Optional[int], stage_channel_name: Optional[str], stage_guild_id: Optional[str], stage_initiator_id: Optional[int]):
        self.channel_id: Optional[int] = stage_channel_id
        self.channel_name: Optional[str] = stage_channel_name
        self.expires_at: Optional[datetime] = stage_expires_at
        self.guild_id: Optional[int] = stage_guild_id
        self.initiator_id: Optional[int] = stage_initiator_id

    @classmethod
    async def fetch_stage_by_guild_id_channel_id_and_channel_name(cls, guild_id: Optional[int], channel_id: Optional[int], channel_name: Optional[str]):
        try:
            bot = DiscordBot.get_instance()
            async with bot.db_pool.acquire() as conn:
                row = await conn.fetchrow('''
                    SELECT expires_at, initiator_id FROM active_stages
                    WHERE channel_id = $1 AND guild_id = $2 AND room_name = $3
                ''', channel_id, guild_id, channel_name)
                if row:
                    return Stage(row['expires_at'], channel_id,  channel_name, guild_id, row['initiator_id'])
                raise Exception(f'No active stage found for <#{channel_id}>.')
        except Exception:
            raise

    @classmethod
    async def fetch_stage_by_guild_id_and_channel_name(cls, guild_id: Optional[int], channel_name: Optional[str]):
        try:
            bot = DiscordBot.get_instance()
            async with bot.db_pool.acquire() as conn:
                row = await conn.fetchrow('''
                    SELECT channel_id, expires_at, initiator_id FROM active_stages
                    WHERE guild_id = $1 AND room_name = $2
                ''', guild_id, channel_name)
                if row:
                    return Stage(row['expires_at'], row['channel_id'], channel_name, guild_id, row['initiator_id'])
                raise Exception(f'No active stage for `{channel_name}`.')
        except Exception:
            raise
    
    @classmethod
    async def fetch_stage_temporary_coordinator_ids_by_guild_id_and_channel_name(cls, guild_id: Optional[int], channel_name: Optional[str]):
        try:
            bot = DiscordBot.get_instance()
            async with bot.db_pool.acquire() as conn:
                row = await conn.fetch('''
                    SELECT discord_snowflake FROM stage_coordinators
                    WHERE guild_id = $1 AND room_name = $2
                ''', guild_id, channel_name)
                if row:
                    temporary_stage_coordinator_ids = {c['discord_snowflake'] for c in row}
                    return temporary_stage_coordinator_ids
                raise Exception('No temporary stage coordinators for this stage.')
        except Exception:
            raise

    @classmethod
    async def fetch_stage_temporary_coordinator_ids_by_channel_id(cls, channel_id: Optional[int]):
        try:
            bot = DiscordBot.get_instance()
            async with bot.db_pool.acquire() as conn:
                row = await conn.fetch('''
                    SELECT discord_snowflake FROM stage_coordinators
                    WHERE channel_id = $1
                ''', channel_id)
                if row:
                    temporary_stage_coordinator_ids = {c['discord_snowflake'] for c in row}
                    return temporary_stage_coordinator_ids
                raise Exception('No temporary stage coordinators for this stage.')
        except Exception:
            raise

    async def send_stage_ask_to_speak_message(self, join_log: dict[int, discord.Member], member: discord.Member):
        bot = DiscordBot.get_instance()
        now = time.time()
        join_log[member.id] = [t for t in join_log[member.id] if now - t < 300]
        if len(join_log[member.id]) < 1:
            join_log[member.id].append(now)
            embed = discord.Embed(
                title=f"\U0001F399 {self.stage_channel_id} — Stage Mode",
                description=f"Ends <t:{int(self.expires_at.timestamp())}:R>",
                color=discord.Color.green()
            )
            embed.add_field(name="\u200b", value="**Ask to speak!**", inline=False)
            await bot.get_channel(self.channel_id).send(embed=embed)

    async def update_stage_by_channel_id_name(self, channel_id: Optional[int], channel_name: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE active_stages SET channel_id=$1
                WHERE guild_id=$2 AND room_name=$3
            ''', channel_id, self.guild_id, channel_name)
            await conn.execute('''
                UPDATE stage_coordinators SET channel_id=$1
                WHERE guild_id=$2 AND room_name=$3
            ''', channel_id.id, self.guild_id, channel_name)

    async def update_stage_by_channel_id_name_initiator_id_and_temporary_coordinator_ids(self, channel_id: Optional[int], channel_name: Optional[str], stage_initiator_id: Optional[int], temporary_stage_coordinator_ids: set[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE active_stages SET channel_id=$1, initiator_id=$3
                WHERE guild_id=$2 AND room_name=$4
            ''', channel_id, self.guild_id, stage_initiator_id, channel_name)
            await conn.execute('''
                UPDATE stage_coordinators SET channel_id=$1
                WHERE discord_snowflake=ANY($2) AND guild_id=$3 AND room_name=$4
            ''', channel_id, list(temporary_stage_coordinator_ids), self.guild_id, channel_name)