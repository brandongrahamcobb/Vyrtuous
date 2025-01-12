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
from discord.ext import commands, menus, tasks
from discord.utils import get
from lucy.utils.create_https_moderation import create_https_moderation
from lucy.utils.create_moderation import create_moderation
from lucy.utils.load_contents import load_contents
from lucy.utils.message import Message
from lucy.utils.nlp_utils import NLPUtils
from lucy.utils.predicator import Predicator
from os.path import abspath, dirname, exists, expanduser, join

import asyncio
import datetime
import discord
import json
import os
import pytz
import random
import shutil
import subprocess
import traceback
from lucy.utils.helpers import *

class Indica(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.conversations = bot.conversations
        self.db_pool = bot.db_pool
        self.handler = Message(self.config, self.conversations)
        self.predicator = Predicator(self.bot)

    @commands.Cog.listener()
    @commands.check(lambda ctx: self.predicator(release_mode))
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)

    async def handle_moderation(self, message: discord.Message):
        if lambda ctx: self.predicator.at_home():
            if not (lambda ctx: self.predicator.is_vegan(message.author)):
                user_id = message.author.id
                async with self.db_pool.acquire() as connection:
                    async with connection.transaction():
                        row = await connection.fetchrow("SELECT flagged_count FROM moderation_counts WHERE user_id = $1", user_id)
                        if row:
                            flagged_count = row['flagged_count'] + 1
                            await connection.execute("UPDATE moderation_counts SET flagged_count = $1 WHERE user_id = $2", flagged_count, user_id)
                        else:
                            flagged_count = 1
                            await connection.execute("INSERT INTO moderation_counts (user_id, flagged_count) VALUES ($1, $2)", user_id, flagged_count)
                if flagged_count == 1:
                    await message.reply("Warning: Your message has been flagged.")
                elif flagged_count in [2, 3, 4]:
                    await message.delete()
                    if flagged_count == 4:
                        await message.author.send("Warning: Your message has been flagged again.")
                elif flagged_count == 5:
                    await message.delete()
                    await message.author.send("You have been timed out for 30 seconds due to repeated violations.")
                    await message.author.timeout(duration=30)  # Timeout for 30 seconds

    @commands.Cog.listener()
    async def on_message(self, message):
        logger.info(f'Received message: {message.content}')
        try:
            if message.author.bot:
                return
            ctx = await self.bot.get_context(message)
            array = await self.handler.process_array(
                message.content, attachments=message.attachments
            )
            if not array:
                logger.error("Invalid 'messages': The array is empty or improperly formatted.")
                if (lambda ctx: self.predicator.at_home()):
                    await message.reply("Your message must include text or valid attachments.")
                return
            logger.info(f"Final payload for processing: {json.dumps(array, indent=2)}")
            for item in array:
                # Moderation
                if self.config['openai_chat_moderation']:
                    async for moderation_completion in create_moderation(input_array=[item]):
                        try:
                            full_response = json.loads(moderation_completion)
                            results = full_response.get('results', [])
                            if results and results[0].get('flagged', False):
                                await self.handle_moderation(message)
                                return
                            else:
                                async for moderation_completion in self.handler.generate_moderation_completion(custom_id=message.author.id, array=array):
                                    chat_moderation = json.loads(moderation_completion)
                                    carnism_results = chat_moderation.get('results', [])
                                    carnism_flagged = carnism_results[0]['categories'].get('carnism', False)
                                    if carnism_flagged:
                                        carnism_score = carnism_results[0]['category_scores'].get('carnism', 0)
                                        await self.handle_moderation(message)
                                        NLPUtils.append_to_jsonl(PATH_TRAINING, carnism_score, message.content, message.author.id)
                                        return
                        except Exception as e:
                            logger.error(traceback.format_exc())
                            if (lambda ctx: self.predicator.at_home()):
                                await message.reply(f'An error occurred: {e}')
                # Chat completion
                if (lambda ctx: self.predicator.is_vegan(message.author)) and message.guild.id != 730907954345279591:
                    if self.config['openai_chat_completion'] and self.bot.user in message.mentions:
                        async for chat_completion in self.handler.generate_chat_completion(
                            custom_id=message.author.id, array=[item], sys_input=OPENAI_CHAT_SYS_INPUT
                        ):
                            if len(chat_completion) > 2000:
                                with open(f'temp.txt', 'w') as f:
                                    f.write(chat_completion)
                                await message.reply(file=discord.File(f'temp.txt'))
                                os.remove(f'temp.txt')
                            else:
                                await message.reply(chat_completion)
        except Exception as e:
            logger.error(traceback.format_exc())
            if (lambda ctx: self.predicator.at_home()):
                await message.reply(f'An error occurred: {e}')
        finally:
            try:
                shutil.rmtree(DIR_TEMP)
                os.makedirs(DIR_TEMP, exist_ok=True)
                logger.info("Temporary files cleaned up successfully.")
            except Exception as cleanup_error:
                logger.error(f"Error cleaning up temporary files: {cleanup_error}")

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
