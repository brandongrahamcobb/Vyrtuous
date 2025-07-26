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
from discord.ext import commands
from vyrtuous.service.discord_message_service import DiscordMessageService
from vyrtuous.inc.helpers import *
from vyrtuous.utils.setup_logging import logger

import discord
import inspect

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
    
        async with self.db_pool.acquire() as conn:
            # Case 1: User was muted and now is not, in a channel
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
    
            # Case 2: User just became muted in a channel
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
    
            # Case 3: User joined or moved to a new channel
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
    



#   @commands.Cog.listener()
#    async def on_voice_state_update(self, member: discord.Member, before: discord.VoiceState, after: discord.VoiceState) -> None:
#        if member.bot:
#            return
#        user_id = member.id
#        before_channel = before.channel
#        after_channel = after.channel
#        async with self.db_pool.acquire() as conn:
#            if before.mute and not after.mute and before_channel:
#                row = await conn.fetchrow("""
#                    SELECT source FROM active_mutes
#                    WHERE user_id = $1 AND channel_id = $2
#                """, user_id, before_channel.id)
#                if row and row['source'] == 'manual':
#                    await conn.execute("""
#                        DELETE FROM active_mutes
#                        WHERE user_id = $1 AND channel_id = $2
#                    """, user_id, before_channel.id)
#                    await conn.execute("""
#                        UPDATE users
#                        SET manual_mute_channels = array_remove(manual_mute_channels, $2),
#                            updated_at = NOW()
#                        WHERE user_id = $1
#                    """, user_id, before_channel.id)
#            if not before.mute and after.mute and after_channel:
#                row = await conn.fetchrow("""
#                    SELECT source FROM active_mutes
#                    WHERE user_id = $1 AND channel_id = $2
#                """, user_id, after_channel.id)
#                if not row:
#                    await conn.execute("""
#                        INSERT INTO active_mutes (user_id, channel_id, source)
#                        VALUES ($1, $2, 'manual')
#                        ON CONFLICT DO NOTHING
#                    """, user_id, after_channel.id)
#                    await conn.execute("""
#                        INSERT INTO users (user_id, manual_mute_channels)
#                        VALUES ($1, ARRAY[$2]::BIGINT[])
#                        ON CONFLICT (user_id) DO UPDATE
#                        SET manual_mute_channels = (
#                            SELECT ARRAY(
#                                SELECT DISTINCT unnest(u.manual_mute_channels || ARRAY[$2])
#                                FROM users u WHERE u.user_id = EXCLUDED.user_id
#                            )
#                        ),
#                        updated_at = NOW()
#                    """, user_id, after_channel.id)
#            if after_channel:
#                row = await conn.fetchrow("""
#                    SELECT source FROM active_mutes
#                    WHERE user_id = $1 AND channel_id = $2
#                """, user_id, after_channel.id)
#                if row:
#                    if row['source'] in ('manual', 'bot') and not after.mute:
#                        await member.edit(mute=True)
#                else:
#                    if after.mute:
#                        await member.edit(mute=False)

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)
                
    async def send_command_help(self, ctx: commands.Context, cmd: commands.Command) -> None:
        embed = discord.Embed(
            title=f'/{cmd.name}',
            description=cmd.help or 'No description provided.',
            color=discord.Color.blue()
        )
    
        sig = inspect.signature(cmd.callback)
        parameters = list(sig.parameters.items())[2:]  # Skip self and ctx
    
        for name, param in parameters:
            is_optional = param.default != inspect.Parameter.empty
            annotation = (
                param.annotation.__name__
                if hasattr(param.annotation, '__name__')
                else str(param.annotation)
            )
            label = 'Optional' if is_optional else 'Required'
            embed.add_field(
                name=f'`{name}`',
                value=f'Type: `{annotation}`\n{label}',
                inline=False
            )
    
        await ctx.send(embed=embed)
    
async def setup(bot: commands.Bot):
    await bot.add_cog(EventListeners(bot))

