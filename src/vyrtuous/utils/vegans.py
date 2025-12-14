
from vyrtuous.bot.discord_bot import DiscordBot
import discord

class Vegans:
    
    vegans = set()
        
    def __init__(self):
        self.bot = DiscordBot.get_instance()
        self.state: bool = False

    async def unrestrict(self, guild: discord.Guild, member: discord.Member):
        uid = member.id
        async with self.bot.db_pool.acquire() as conn:
            ban_rows = await conn.fetch('SELECT channel_id FROM active_bans WHERE discord_snowflake=$1', uid)
            mute_rows = await conn.fetch('SELECT channel_id FROM active_voice_mutes WHERE discord_snowflake=$1', uid)
            text_rows = await conn.fetch('SELECT channel_id FROM active_text_mutes WHERE discord_snowflake=$1', uid)
        for r in ban_rows:
            try: await guild.unban(discord.Object(id=uid), reason='Toggle OFF')
            except: pass
        for r in mute_rows:
            ch = guild.get_channel(r['channel_id'])
            if ch and member.voice and member.voice.mute: await member.edit(mute=False)
        for r in text_rows:
            ch = guild.get_channel(r['channel_id'])
            text_mute_role = discord.utils.get(guild.roles, name='TextMuted')
            if ch and text_mute_role and text_mute_role in member.roles: await member.remove_roles(text_mute_role)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_bans WHERE discord_snowflake=$1', uid)
            await conn.execute('DELETE FROM active_voice_mutes WHERE discord_snowflake=$1', uid)
            await conn.execute('DELETE FROM active_text_mutes WHERE discord_snowflake=$1', uid)

    @classmethod
    def get_vegans(cls):
        return cls.vegans
