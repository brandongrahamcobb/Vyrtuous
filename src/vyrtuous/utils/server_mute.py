from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
import asyncpg

class ServerMute:

    PLURAL = "Server Mutes"
    SINGULAR = "Server Mute"

    def __init__(self, guild_snowflake: Optional[int], member_snowflake: Optional[int], reason: Optional[str]):
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.reason = reason

    async def grant(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_server_voice_mutes (created_at, guild_snowflake, member_snowflake, reason)
                VALUES (NOW(), $1, $2, $3)
                ON CONFLICT (guild_snowflake, member_snowflake) DO UPDATE
                SET reason = EXCLUDED.reason
            ''', self.guild_snowflake, self.member_snowflake, self.reason)

    @classmethod
    async def delete_by_guild_and_member(cls, guild_snowflake: Optional[int], member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_server_voice_mutes WHERE guild_snowflake = $1 AND member_snowflake = $2', guild_snowflake, member_snowflake)

    @classmethod
    async def fetch_by_member(self, member_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT guild_snowflake, reason FROM active_server_voice_mutes WHERE member_snowflake = $1
            ''', member_snowflake)
        if row:
            return ServerMute(guild_snowflake=row['guild_snowflake'], member_snowflake=member_snowflake, reason=row["reason"])