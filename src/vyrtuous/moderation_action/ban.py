''' ban.py The purpose of this program is to inherit from the moderation event class to provide the ban moderation.

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
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.moderation_action.moderation_action import ModerationAction
from vyrtuous.utils.history import History

class Ban(ModerationAction):

    PLURAL = "Bans"
    SINGULAR = "Ban"

    def __init__(self, channel_snowflake: Optional[int], expires_in: Optional[datetime], guild_snowflake: Optional[int], member_snowflake: Optional[int], reason: Optional[str]):
        self.channel_snowflake = channel_snowflake
        self.expires_in = expires_in
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.reason = reason
 
    @classmethod
    async def delete_by_channel_and_guild(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
              await conn.execute('''
                DELETE FROM active_bans
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
    
    @classmethod
    async def delete_by_channel_guild_and_member(self, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
              await conn.execute('''
                DELETE FROM active_bans
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND member_snowflake=$3
            ''', channel_snowflake, guild_snowflake, member_snowflake)

    @classmethod
    async def delete_by_guild_and_member(self, guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
              await conn.execute('''
                DELETE FROM active_bans
                WHERE guild_snowflake=$1 AND member_snowflake=$2
            ''', guild_snowflake, member_snowflake)

    @classmethod
    async def delete_by_guild(self, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
              await conn.execute('''
                DELETE FROM active_bans
                WHERE guild_snowflake=$1
            ''', guild_snowflake)

    async def create(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_bans (channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason)
                VALUES ($1, NOW(), $2, $3, $4, $5)
                ON CONFLICT DO NOTHING
            ''', self.channel_snowflake, self.expires_in, self.guild_snowflake, self.member_snowflake, self.reason)

    @classmethod
    async def update_by_source_and_target(cls, source_channel_snowflake: Optional[int], target_channel_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE active_bans
                SET channel_snowflake=$2
                WHERE channel_snowflake=$1
            ''', source_channel_snowflake, target_channel_snowflake)
    
    @classmethod
    async def expired_ban(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason, updated_at
                FROM active_bans
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND member_snowflake=$3
            ''', channel_snowflake, guild_snowflake, member_snowflake)
        if row:
            return Ban(channel_snowflake=row['channel_snowflake'], expires_in=row['expires_in'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason'])

    @classmethod
    async def fetch_by_channel_guild_and_member(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason, updated_at
                FROM active_bans
                WHERE channel_snowflake=$1 AND guild_snowflake=$2 AND member_snowflake=$3
            ''', channel_snowflake, guild_snowflake, member_snowflake)
        if row:
            return Ban(channel_snowflake=row['channel_snowflake'], expires_in=row['expires_in'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason'])

    @classmethod
    async def fetch_by_guild_and_member(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason, updated_at
                FROM active_bans
                WHERE guild_snowflake=$1 AND member_snowflake=$2
            ''', guild_snowflake, member_snowflake)
        bans = []
        if rows:
            for row in rows:
                bans.append(Ban(channel_snowflake=row['channel_snowflake'], expires_in=row['expires_in'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason']))
        return bans

    @classmethod
    async def fetch_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason, updated_at
                FROM active_bans
                WHERE guild_snowflake=$1
            ''', guild_snowflake)
        bans = []
        if rows:
            for row in rows:
                bans.append(Ban(channel_snowflake=row['channel_snowflake'], expires_in=row['expires_in'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason']))
        return bans

    @classmethod
    async def fetch_all(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason, updated_at
                FROM active_bans
            ''')
        bans = []
        if rows:
            for row in rows:
                bans.append(Ban(channel_snowflake=row['channel_snowflake'], expires_in=row['expires_in'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason']))
        return bans
    
    @classmethod
    async def fetch_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason, updated_at
                FROM active_bans
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
        bans = []
        if rows:
            for row in rows:
                bans.append(Ban(channel_snowflake=channel_snowflake, expires_in=row['expires_in'], guild_snowflake=guild_snowflake, member_snowflake=row['member_snowflake'], reason=row['reason']))
        return bans

    @classmethod
    async def fetch_by_expired(cls, now: Optional[datetime]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason
                FROM active_bans
                WHERE expires_in IS NOT NULL AND expires_in <= $1
            ''', now)
        expired_bans = []
        if rows:
            for row in rows:
                expired_bans.append(Ban(channel_snowflake=row['channel_snowflake'], expires_in=row['expires_in'], guild_snowflake=row['guild_snowflake'], member_snowflake=row['member_snowflake'], reason=row['reason']))
        return expired_bans

    @classmethod
    async def clear_by_channel_guild_highest_role_and_modification(cls, ctx_interaction_or_message, channel_snowflake: Optional[int], guild_snowflake: Optional[int], highest_role: Optional[str], is_modification: bool):
        bans = await cls.fetch_by_channel_and_guild(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake)
        await cls.delete_by_channel_and_guild(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake)
        if bans:
            for ban in bans:
                await History.save_entry(
                    ctx_interaction_or_message=ctx_interaction_or_message,
                    action_type='unban',
                    channel_snowflake=channel_snowflake,
                    duration=None,
                    highest_role=highest_role,
                    is_modification=is_modification,
                    member_snowflake=voice_mute.member_snowflake,
                    reason="Clear command"
                )
    
    @classmethod
    async def clear_by_guild_highest_role_member_and_modification(cls, ctx_interaction_or_message, guild_snowflake: Optional[int], highest_role: Optional[str], is_modification: bool, member_snowflake: Optional[int]):
        bans = await cls.fetch_by_guild_and_member(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
        await cls.delete_by_guild_and_member(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
        if bans:
            for ban in bans:
                await History.save_entry(
                    ctx_interaction_or_message=ctx_interaction_or_message,
                    action_type='unban',
                    channel_snowflake=ban.channel_snowflake,
                    duration=None,
                    highest_role=highest_role,
                    is_modification=is_modification,
                    member_snowflake=member_snowflake,
                    reason="Clear command"
                )