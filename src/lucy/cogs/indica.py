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
from lucy.utils.combine import combine
from lucy.utils.create_https_moderation import create_https_moderation
from lucy.utils.create_moderation import create_moderation
from lucy.utils.draw_fingerprint import draw_fingerprint
from lucy.utils.helpers import *
from lucy.utils.game import Game
from lucy.utils.load_contents import load_contents
from lucy.utils.message import Message
from lucy.utils.nlp_utils import NLPUtils
from lucy.utils.predicator import Predicator
from lucy.utils.tag import TagManager
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
import uuid
import yaml

class Indica(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.conversations = bot.conversations
        self.db_pool = bot.db_pool
        self.handler = Message(self.config, self.conversations)
        self.predicator = Predicator(self.bot)
        self.game = Game(self.bot)

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if await self.predicator.is_at_home_func(before.guild.id):
            if before.content != after.content:
                ctx = await self.bot.get_context(after)
                if ctx.command:
                    await self.bot.invoke(ctx)

    async def ai_handler(self, ctx: commands.Context):
        array = await self.handler.process_array(ctx.message.content, attachments=ctx.message.attachments)
        
        if not array:
            logger.error("Invalid 'messages': The array is empty or improperly formatted.")
            await self.send_message(ctx, "Your message must include text or valid attachments.")
            return
    
        logger.info(f"Final payload for processing: {json.dumps(array, indent=2)}")
    
        for item in array:
#            if self.config['openai_chat_completion'] and self.bot.user in ctx.message.mentions and not self.predicator.is_spawd(ctx):
#                async for chat_completion in self.handler.generate_chat_completion(
#                    custom_id=ctx.author.id, array=[item], sys_input=OPENAI_CHAT_SYS_INPUT
#                ):
#                    await self.handle_large_response(ctx, chat_completion)
    
            # Moderation Check
            if self.bot.user in ctx.message.mentions and self.config['openai_chat_moderation']:
                async for moderation_completion in create_moderation(input_array=[item]):
                    try:
                        full_response = json.loads(moderation_completion)
                        results = full_response.get('results', [])
                        
                        if results and results[0].get('flagged', False) and not self.predicator.is_spawd(ctx):
                            await self.handle_moderation(ctx.message)
                            return
                        
                    except Exception as e:
                        logger.error(traceback.format_exc())
                        if await self.predicator.is_at_home_func(ctx.guild.id):
                            print(f'An error occurred: {e}')
        if self.config['openai_chat_completion'] and self.bot.user in ctx.message.mentions:
#                        elif self.predicator.is_spawd(ctx) and self.bot.user in ctx.message.mentions:
            async for chat_completion in self.handler.generate_chat_completion(
                custom_id=ctx.author.id, array=array, sys_input=OPENAI_CHAT_SYS_INPUT
            ):
                await self.handle_large_response(ctx, chat_completion)
    
    async def handle_large_response(self, ctx: commands.Context, response: str):
        """Handles sending messages larger than 2000 characters by creating temp files."""
        if len(response) > 2000:
            unique_filename = f'temp_{uuid.uuid4()}.txt'
            
            try:
                with open(unique_filename, 'w') as f:
                    f.write(response)
                await self.send_file(ctx, file=discord.File(unique_filename))
            finally:
                os.remove(unique_filename)
        else:
            await self.send_message(ctx, response)

    async def handle_moderation(self, message: discord.Message):
        if await self.predicator.is_at_home_func(message.guild.id):
            if not await self.predicator.is_vegan_user(message.author):
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

    def handle_users(self, author: str):
        author_char = author[0].upper()  # Fixed: Get the first character in uppercase
        data = {letter: [] for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}
        users_file = join(DIR_HOME, '.users', 'users')
        if exists(users_file):
            with open(users_file, 'r+') as file:  # Fixed: 'r+' allows reading and writing
                try:
                    data = yaml.safe_load(file) or data  # Handle empty file case
                except yaml.YAMLError:
                    data = {letter: [] for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}  # Reset on error
                if author_char not in data:
                    data[author_char] = []  # Ensure the key exists
                if author not in data[author_char]:
                    data[author_char].append(author)  # Append user if not already there
                    file.seek(0)  # Move to beginning before writing
                    file.write(yaml.dump(data))
                    file.truncate()  # Truncate extra old content
        else:
            with open(users_file, 'w') as file:  # 'w' to create new file
                yaml.dump(data, file)

    async def send_message(self, ctx: commands.Context, print: str):
        if await self.predicator.is_at_home_func(ctx.guild.id):
            await ctx.reply(print)

    async def send_file(self, ctx: commands.Context, file: discord.File):
        if await self.predicator.is_at_home_func(ctx.guild.id):
            await ctx.reply(file=file)

#    async def ai_handler(self, ctx: commands.Context):
#        array = await self.handler.process_array(ctx.message.content, attachments=ctx.message.attachments)
#        if not array:
#            logger.error("Invalid 'messages': The array is empty or improperly formatted.")
#            await self.send_message(ctx, "Your message must include text or valid attachments.")
#            return
#        logger.info(f"Final payload for processing: {json.dumps(array, indent=2)}")
#        for item in array:
#            if self.config['openai_chat_completion'] and self.bot.user in ctx.message.mentions and not self.predicator.is_spawd(ctx):
#                async for chat_completion in self.handler.generate_chat_completion(
#                     custom_id=ctx.author.id, array=[item], sys_input=OPENAI_CHAT_SYS_INPUT
#                ):
#                    if len(chat_completion) > 2000:
#                        unique_filename = f'temp_{uuid.uuid4()}.txt'
#                        with open(unique_filename, 'w') as f:
#                             f.write(chat_completion)
#                        await self.send_file(ctx, file=discord.File(unique_filename))
#                        os.remove(unique_filename)
#                    else:
#                        await self.send_message(ctx, chat_completion)
#            # Moderation
#            if self.bot.user in ctx.message.mentions and self.config['openai_chat_moderation']:
#                async for moderation_completion in create_moderation(input_array=[item]):
#                    try:
#                        full_response = json.loads(moderation_completion)
#                        results = full_response.get('results', [])
#                        if results and results[0].get('flagged', False) and not self.predicator.is_spawd(ctx):
#                            await self.handle_moderation(ctx.message)
#                            return
#                        elif self.predicator.is_spawd(ctx) and self.bot.user in ctx.message.mentions:
#                            async for chat_completion in self.handler.generate_chat_completion(
#                                custom_id=ctx.author.id, array=[item], model='o1-mini'
#                            ):
#                                if len(chat_completion) > 2000:
#                                    unique_filename = f'temp_{uuid.uuid4()}.txt'
#                                    with open('part_1_' + unique_filename, 'w') as f:
#                                        f.write(chat_completion[:1000])
#                                    await send_file(ctx, file=discord.File('part_1_' + unique_filename))
#                                    os.remove('part_1_' + unique_filename)
#                                    with open('part_2_' + unique_filename, 'w') as f:
#                                       f.write(chat_completion[1000:])
#                                    await send_file(ctx, file=discord.File('part_2_' + unique_filename))
#                                    os.remove('part_2_' + unique_filename)
#                                else:
#                                    await self.send_message(ctx, chat_completion)
#                    except Exception as e:
#                        logger.error(traceback.format_exc())
#                        if await self.predicator.is_at_home_func(ctx.guild.id):
#                            print(f'An error occurred: {e}')
#            async for moderation_completion in self.handler.generate_moderation_completion(custom_id=ctx.author.id, array=array):
#                chat_moderation = json.loads(moderation_completion)
#                academic_dishonesty_results = chat_moderation.get('results', [])
#                academic_dishonesty_flagged = academic_dishonesty_results[0]['categories'].get('academic-dishonesty', False)
#                if academic_dishonesty_flagged:
#                    academic_dishonesty_score = academic_dishonesty_results[0]['category_scores'].get('academic-dishonesty', 0)
#                    await self.handle_moderation(message)
#                    NLPUtils.append_to_jsonl(PATH_TRAINING, academic_dishonesty_score, ctx.message.content, ctx.author.id)
#                    return
#

    @commands.Cog.listener()
    async def on_message(self, message):
        logger.info(f'Received message: {message.content}')
        try:
            if message.author == self.bot.user:
                return
            ctx = await self.bot.get_context(message)
            author = ctx.author.name
            self.handle_users(author)
            self.game.distribute_xp(ctx.author.id)
        except Exception as e:
            logger.error(traceback.format_exc())
            if await self.predicator.is_at_home_func(message.guild.id):
                print(f'An error occurred: {e}')
        finally:
            try:
                shutil.rmtree(DIR_TEMP)
                os.makedirs(DIR_TEMP, exist_ok=True)
                logger.info("Temporary files cleaned up successfully.")
            except Exception as cleanup_error:
                logger.error(f"Error cleaning up temporary files: {cleanup_error}")

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
