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

#    async def on_voice_state_update(self, member: discord.Member, before: discord.VoiceState, after: discord.VoiceState) -> None:
#        if member.bot:
#            return
#        user_id = member.id
#        before_channel = before.channel
#        after_channel = after.channel
#        async with self.db_pool.acquire() as conn:
#            if before.mute and not after.mute and before_channel:
#                row = await conn.fetchrow("""
#                    SELECT source, issuer_id FROM active_mutes
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
#                    SELECT source, issuer_id FROM active_mutes
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
#                    SELECT source, issuer_id FROM active_mutes
#                    WHERE user_id = $1 AND channel_id = $2
#                """, user_id, after_channel.id)
#                if row:
#                    if row['source'] == 'owner':
#                        if not after.mute:
#                            await member.edit(mute=True)
#                    elif row['source'] in ('manual', 'bot') and not after.mute:
#                        await member.edit(mute=True)
#                else:
#                    if after.mute:
#                        await member.edit(mute=False)
#        
#            if after_channel and before_channel != after_channel:
#                is_flagged = await conn.fetchval("""
#                    SELECT flagged FROM users WHERE user_id = $1 AND flagged = TRUE
#                """, user_id)
#                if is_flagged:
#                    linked_text_channel = discord.utils.get(
#                        after_channel.guild.text_channels,
#                        name=after_channel.name
#                    )
#                    if not linked_text_channel and after_channel.category:
#                        for tc in after_channel.category.text_channels:
#                            if tc.permissions_for(after_channel.guild.me).send_messages:
#                                linked_text_channel = tc
#                                break
#                    if linked_text_channel:
#                        try:
#                            await linked_text_channel.send(
#                                f'🚩 **Flagged user joined**: {member.mention} has joined {after_channel.mention}.'
#                            )
#                        except discord.Forbidden:
#                            pass

    @commands.Cog.listener()
    async def on_voice_state_update(self, member: discord.Member, before: discord.VoiceState, after: discord.VoiceState) -> None:
        if member.bot:
            return
        user_id = member.id
        before_channel = before.channel
        after_channel = after.channel
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
                        INSERT INTO active_mutes (user_id, channel_id, source)
                        VALUES ($1, $2, 'manual')
                        ON CONFLICT DO NOTHING
                    """, user_id, after_channel.id)
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
        
    @commands.Cog.listener()
    async def on_command_error(self, ctx, error) -> None:
#        if isinstance(error, commands.MissingRequiredArgument):
#            await self.send_command_help(ctx, ctx.command)
        if isinstance(error, commands.CheckFailure):
            await self.handler.send_message(ctx, content=str(error))
            
async def setup(bot: commands.Bot):
    await bot.add_cog(EventListeners(bot))
