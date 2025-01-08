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
from collections import defaultdict
from discord.ext import commands, tasks
from lucy.utils.backup import perform_backup, setup_backup_directory
from lucy.utils.create_https_completion import Conversations
from lucy.utils.create_https_moderation import create_https_moderation
from lucy.utils.create_moderation import create_moderation
from lucy.utils.nlp_utils import NLPUtils
from lucy.utils.load_contents import load_contents
from lucy.utils.message import Message
from os.path import abspath, dirname, exists, expanduser, join

import asyncio
import datetime
import discord
import json
import os
import random
import subprocess
import traceback
from lucy.utils.helpers import *

def at_home():
    async def predicate(ctx):
        return ctx.guild is not None and ctx.guild.id == 1300517536001036348
    return commands.check(predicate)

class Indica(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.conversations = Conversations()
        self.handler = Message(self.config, self.conversations)

    def cog_unload(self):
        """Stop the backup task when the Cog is unloaded."""
        self.backup_task.cancel()

    @tasks.loop(hours=24)  # Adjust the interval as needed
    async def backup_task(self):
        """Periodic task to back up the database."""
        try:
            # Ensure the backup directory exists
            backup_dir = setup_backup_directory("./backups")

            # Perform the database backup
            backup_file = perform_backup(
                db_user="postgres",
                db_name="lucy",
                db_host="localhost",
                backup_dir=backup_dir
            )

            logger.info(f"Backup completed successfully: {backup_file}")
        except Exception as e:
            logger.error(f"Error during database backup: {e}")

    @backup_task.before_loop
    async def before_backup(self):
        """Wait until the bot is ready before starting the task."""
        await self.bot.wait_until_ready()

    @commands.Cog.listener()
    @at_home()
    async def on_message_edit(self, before, after):

        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)

    @commands.Cog.listener()
    async def on_message(self, message):
        logger.info(f"Received message: {message.content}")
        try:
            if message.author.bot:
                return
            if message.attachments:
                array = await self.handler.process_array(message.content, attachments=message.attachments)
            else:
                array = await self.handler.process_array(message.content)
            # Chat
            if message.guild.id == 1300517536001036348:
                if self.bot.user in message.mentions:
                    if self.config['openai_chat_completion']:
                        async for chat_completion in self.handler.generate_chat_completion(custom_id=message.author.id, array=array, sys_input=OPENAI_CHAT_SYS_INPUT):
                            await message.reply(chat_completion)

            # Moderate Text and Images
            if self.config['openai_chat_moderation']:
                async for moderation_completion in self.handler.generate_moderation_completion(custom_id=message.author.id, array=array):
                    full_response = json.loads(moderation_completion)
                    results = full_response.get('results', [])
                    flagged = results[0].get('flagged', False)
                    carnism_flagged = results[0]['categories'].get('carnism', False)
                    carnism_score = results[0]['category_scores'].get('carnism', 0)
                    total_carnism_score = sum(arg['category_scores'].get('carnism', 0) for arg in results)
                    if carnism_flagged or flagged:
                        if message.guild.id == 1300517536001036348:
                            guild = await self.bot.fetch_guild(self.config['discord_testing_guild_id'])
                            role = guild.get_role(self.config['discord_role_pass'])
                            if not role in message.author.roles and flagged:
                                await message.delete()
                        NLPUtils.append_to_other_jsonl(PATH_TRAINING, carnism_score, message.content, message.author.id)
        except Exception as e:
            logger.error(traceback.format_exc())
            await message.reply(e)

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
