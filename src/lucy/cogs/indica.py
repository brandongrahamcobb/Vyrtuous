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
from discord.ext import commands
from lucy.utils.helpers import *
from lucy.utils.game import Game
from lucy.utils.message import Message
from lucy.utils.setup_logging import logger

import discord
import os
import shutil
import time
import traceback
import uuid

class Indica(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.conversations = bot.conversations
        self.db_pool = bot.db_pool
        self.game = Game(self.bot)
        self.handler = Message(self.bot, self.config, self.conversations, self.db_pool)
        self.user_messages = {}

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)

    @commands.Cog.listener()
    async def on_member_update(self, before: discord.Member, after: discord.Member):
        if before.nick == after.nick:
            return
        try:
            flagged = await self.game.moderate_name(after.nick or after.name)
            if flagged:
                try:
                    await after.edit(nick=None, reason="Nickname reverted due to moderation violation.")
                except discord.Forbidden:
                    logger.error(f"[moderation] Missing permissions to revert nickname for {after}.")
                except discord.HTTPException as e:
                    logger.error(f"[moderation] HTTPException while reverting nickname for {after}: {e}")
                return
            user_data = await self.game.get_user(after.id)
            if not user_data:
                return
            faction_name = user_data["faction_name"]
            if not faction_name:
                return
            expected_nick = f"[{faction_name}] {after.name}"
            if after.nick != expected_nick:
                try:
                    await after.edit(nick=expected_nick, reason="Enforcing faction nickname format.")
                except discord.Forbidden:
                    logger.error(f"[faction] Cannot change nickname for {after} due to missing permissions.")
                except discord.HTTPException as e:
                    logger.error(f"[faction] HTTPException while changing nickname for {after}: {e}")
        except Exception as e:
            logger.error(traceback.format_exc())
        finally:
            try:
                shutil.rmtree(DIR_TEMP)
                os.makedirs(DIR_TEMP, exist_ok=True)
            except Exception as cleanup_error:
                logger.error(f"Error cleaning up temporary files: {cleanup_error}")

    @commands.Cog.listener()
    async def on_message(self, message):
        try:
            if message.author.bot or message.is_system():
                return
            ctx = await self.bot.get_context(message)
            author = ctx.author.name
            current_time = time.time()
            ctx = await self.bot.get_context(message)
            author = ctx.author.name
            user_id = ctx.author.id
            if user_id not in self.user_messages:
                self.user_messages[user_id] = []
                self.user_messages[user_id].append(current_time)
                self.user_messages[user_id] = [t for t in self.user_messages[user_id] if current_time - t < 5]
            if len(self.user_messages[user_id]) > 5:
                await self.handler.send_message(ctx, content=f"{message.author.mention}, stop spamming!")
                await message.delete()
            self.handler.handle_users(message.author.name)
            await self.game.distribute_xp(ctx.author.id)
            await self.handler.ai_handler(ctx)
        except Exception as e:
            logger.error(traceback.format_exc())
        finally:
            try:
                shutil.rmtree(DIR_TEMP)
                os.makedirs(DIR_TEMP, exist_ok=True)
            except Exception as cleanup_error:
                logger.error(f"Error cleaning up temporary files: {cleanup_error}")

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
