''' server_mute.py The purpose of this program is to inherit from the moderation event class to provide the server mute moderation.

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
from vyrtuous.moderation_action.moderation_action import ModerationAction

class ServerMute(ModerationAction):

    PLURAL = "Server Mutes"
    SINGULAR = "Server Mute"

    def __init__(self, guild_snowflake: Optional[int], member_snowflake: Optional[int], reason: Optional[str]):
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.reason = reason

    async def create(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_server_voice_mutes (created_at, guild_snowflake, member_snowflake, reason)
                VALUES (NOW(), $1, $2, $3)
                ON CONFLICT DO NOTHING
            ''', self.guild_snowflake, self.member_snowflake, self.reason)

    @classmethod
    async def delete_by_guild_and_member(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM active_server_voice_mutes
                WHERE guild_snowflake=$1 AND member_snowflake=$2
            ''', guild_snowflake, member_snowflake)

    @classmethod
    async def fetch_by_member(self, member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT created_at, guild_snowflake, member_snowflake, reason, updated_at
                FROM active_server_voice_mutes WHERE member_snowflake=$1
            ''', member_snowflake)
        if row:
            return ServerMute(guild_snowflake=row['guild_snowflake'], member_snowflake=member_snowflake, reason=row["reason"])