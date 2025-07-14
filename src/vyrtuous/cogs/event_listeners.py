''' indica.py The purpose of the program is to be an extension for a Discord bot for listeners.
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
import os
import shutil
import traceback

import discord
from discord.ext import commands
from vyrtuous.utils.handlers.message_service import MessageService
from vyrtuous.utils.handlers.mute_service import MuteService
from vyrtuous.utils.handlers.predicator import Predicator
from vyrtuous.utils.inc.helpers import *
from vyrtuous.utils.inc.setup_logging import logger
from vyrtuous.utils.inc.helpers import *


    
class Indica(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.handler = MessageService(self.bot, self.config, self.db_pool)
        self.predicator = Predicator(self.bot)
        self.user_messages = {}

    @commands.after_invoke
    async def after_invoke(ctx) -> None:
        if hasattr(bot, 'db_pool'):
            await bot.db_pool.close()
            
    @commands.Cog.listener()
    async def on_command_error(self, ctx, error) -> None:
        if isinstance(error, commands.CheckFailure):
            await self.handler.send_message(ctx, content=str(error))

    @commands.Cog.listener()
    async def on_voice_state_update(self, member: discord.Member, before: discord.VoiceState, after: discord.VoiceState) -> None:
        if member.bot:
            return
        user_id = member.id
        guild_id = member.guild.id
        before_channel = before.channel
        after_channel = after.channel
        async with self.db_pool.acquire() as conn:
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
            if before_channel and not after_channel and before.mute:
                await conn.execute("""
                    INSERT INTO users (user_id, mute_channel_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (user_id) DO UPDATE
                    SET mute_channel_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(u.mute_channel_ids || EXCLUDED.mute_channel_ids)
                        )
                        FROM users u WHERE u.user_id = EXCLUDED.user_id
                    ),
                    updated_at = NOW()
                """, user_id, before_channel.id)
                await member.edit(mute=False)
            if after_channel:
                row = await conn.fetchrow("""
                    SELECT source FROM active_mutes
                    WHERE user_id = $1 AND channel_id = $2
                """, user_id, after_channel.id)
                if row and row['source'] == 'bot':
                    if not after.mute:
                        await member.edit(mute=True)
                elif row and row['source'] == 'manual':
                    pass
                else:
                    if after.mute:
                        await member.edit(mute=False)
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
                    
    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)
                
#    @commands.Cog.listener()
#    async def on_message(self, message):
#        try:
#            if message.author.id == message.is_system():
#                return
#            ctx = await self.bot.get_context(message)
#            author = ctx.author.name
#            return
#        except Exception as e:
#            logger.error(traceback.format_exc())
#        finally:
#            try:
#                shutil.rmtree(DIR_TEMP)
#                os.makedirs(DIR_TEMP, exist_ok=True)
#            except Exception as cleanup_error:
#                logger.error(f"Error cleaning up temporary files: {cleanup_error}")

    @commands.Cog.listener()
    async def on_ready(self):
        bot_user = self.bot.user
        bot_name = bot_user.name
        bot_id = bot_user.id
        guild_count = len(self.bot.guilds)
        info = (
            f'\n=============================\n'
            f'bot Name: {bot_name}\n'
            f'bot ID: {bot_id}\n'
            f'Connected Guilds: {guild_count}\n'
            f'============================='
        )
        guild_info = '\n'.join(
            [f'- {guild.name} (ID: {guild.id})' for guild in self.bot.guilds]
        )
        stats_message = f'{info}\n\nGuilds:\n{guild_info}'
        print(stats_message)

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
