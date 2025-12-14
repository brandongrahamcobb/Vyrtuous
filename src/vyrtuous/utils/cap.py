
from vyrtuous.bot.discord_bot import DiscordBot

class Cap:
 
    @classmethod
    async def get_caps_for_channel(self, guild_id: int, channel_id: int) -> list[tuple[str, str]]:
        self.bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT moderation_type, duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2',
                guild_id, channel_id
            )
            return [(r['moderation_type'], r['duration']) for r in rows]
