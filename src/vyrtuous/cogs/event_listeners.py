''' event_listeners.py

    Copyright (C) 2024  github.com/brandongrahamcobb

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
import discord
import os
from discord.ext import commands
from vyrtuous.inc.helpers import *

from vyrtuous.service.discord_message_service import DiscordMessageService


class EventListeners(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.db_pool)

    @commands.Cog.listener()
    async def on_voice_state_update(
        self,
        member: discord.Member,
        before: discord.VoiceState,
        after: discord.VoiceState
    ) -> None:
        if member.bot:
            return
        user_id = member.id
        before_channel = before.channel
        after_channel = after.channel
        if before.channel == after.channel or after.channel is None:
            return
        guild = member.guild
        channel = after.channel
        channel_id = channel.id
        if os.getenv("DEVELOPMENT") == "False":
            async with self.db_pool.acquire() as conn:
                if before.mute and not after.mute and before_channel:
                    row = await conn.fetchrow("""
                        SELECT source FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2
                    """, user_id, before_channel.id)
                    if row and row['source'] == 'manual':
                        await conn.execute("""
                            DELETE FROM active_mutes
                            WHERE user_id = $1 AND channel_id = $2
                        """, user_id, before_channel.id)
                        await conn.execute("""
                            UPDATE users
                            SET manual_mute_channels = array_remove(manual_mute_channels, $2),
                                updated_at = NOW()
                            WHERE user_id = $1
                        """, user_id, before_channel.id)
                if not before.mute and after.mute and after_channel:
                    row = await conn.fetchrow("""
                        SELECT source FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2
                    """, user_id, after_channel.id)
                    if not row:
                        await conn.execute("""
                            INSERT INTO active_mutes (user_id, channel_id, source, issuer_id)
                            VALUES ($1, $2, 'manual', $3)
                            ON CONFLICT DO NOTHING
                        """, user_id, after_channel.id, member.guild.owner_id)
                        await conn.execute("""
                            INSERT INTO users (user_id, manual_mute_channels)
                            VALUES ($1, ARRAY[$2]::BIGINT[])
                            ON CONFLICT (user_id) DO UPDATE
                            SET manual_mute_channels = (
                                SELECT ARRAY(
                                    SELECT DISTINCT unnest(u.manual_mute_channels || ARRAY[$2])
                                    FROM users u WHERE u.user_id = EXCLUDED.user_id
                                )
                            ),
                            updated_at = NOW()
                        """, user_id, after_channel.id)
                if after_channel:
                    row = await conn.fetchrow("""
                        SELECT source FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2
                    """, user_id, after_channel.id)
                    if row:
                        if row['source'] in ('manual', 'bot') and not after.mute:
                            await member.edit(mute=True)
                    else:
                        if after.mute:
                            await member.edit(mute=False)
        try:
            async with self.bot.db_pool.acquire() as conn:
                is_flagged = await conn.fetchval(
                    '''
                    SELECT 1 FROM users
                    WHERE user_id = $1
                      AND $2 = ANY(flagged_channel_ids)
                    ''',
                    member.id, channel_id
                )
                if not is_flagged:
                    return
                if isinstance(channel, discord.VoiceChannel):
                    await channel.send(
                    f'⚠️ <@{member.id}> has joined voice channel <#{channel.id}> and is flagged.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
        except Exception as e:
            print(f"🔥 Error in on_voice_state_update: {e}")

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)
                
    @commands.Cog.listener()
    async def on_command_error(self, ctx, error):
        if isinstance(error, commands.MissingRequiredArgument):
            missing = error.param.name
            await ctx.reply(f"❌ Missing required argument: `{missing}`")
            return

    
async def setup(bot: commands.Bot):
    await bot.add_cog(EventListeners(bot))

