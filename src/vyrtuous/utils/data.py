"""!/bin/python3

data.py The purpose of this program is to manage statistics of Vyrtuous.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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
"""

from datetime import datetime

from vyrtuous.bot.discord_bot import DiscordBot


class Data:

    @classmethod
    async def save(
        cls,
        identifier: str,
        channel_members_voice_count: int,
        channel_snowflake: int,
        executor_highest_role: str,
        executor_member_snowflake: int,
        expires_at: datetime,
        guild_members_offline_and_online_member_count: int,
        guild_members_online_count: int,
        guild_members_voice_count: int,
        guild_snowflake: int,
        is_modification: bool,
        target_member_snowflake: int,
        target_highest_role: str,
        reason: str,
    ):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                """
                INSERT INTO moderation_logs (infraction_type, channel_members_voice_count, channel_snowflake, executor_highest_role, executor_member_snowflake, expires_at, guild_members_offline_and_online_member_count, guild_members_online_count, guild_members_voice_count, guild_snowflake, is_modification, target_highest_role, target_member_snowflake, reason)
                VALUES ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14)
            """,
                identifier,
                channel_members_voice_count,
                channel_snowflake,
                executor_highest_role,
                executor_member_snowflake,
                expires_at,
                guild_members_offline_and_online_member_count,
                guild_members_online_count,
                guild_members_voice_count,
                guild_snowflake,
                is_modification,
                target_highest_role,
                target_member_snowflake,
                reason,
            )
