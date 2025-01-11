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
from lucy.utils.get_mol import get_mol
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

class Indica(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.conversations = bot.conversations
        self.db_pool = bot.db_pool
        self.handler = Message(self.config, self.conversations)
<<<<<<< HEAD
        self.predicator = Predicator(self.bot)
=======
        self.tag_manager = TagManager(self.bot.db_pool)
        self.daily_loop.start()
        self.channel_guild_map: Dict[int, int] = {
            798967615636504657: 730907954345279591,
            730907954877956179: 730907954345279591,
            1315735859848544378: 1300517536001036348,
        }
        self.guild_loops_index = defaultdict(int)


    def at_home(self):
        async def predicate(ctx):
            if ctx.guild is None:
                logger.warning('at_home called outside a guild context.')
                return False
            testing_guild_id = self.bot.config['discord_testing_guild_id']
            logger.info('Checking guild ID: {ctx.guild.id} against testing guild ID: {testing_guild_id}')
            return ctx.guild.id == testing_guild_id
        return commands.check(predicate)

    def release_mode(self):
        async def predicate(ctx):
            logger.info(f"Checking user ID: {ctx.author.id}")
            logger.info(f"Release mode setting: {self.bot.config['discord_release_mode']}")
            return ctx.author.id == 154749533429956608 or self.bot.config['discord_release_mode'] or isinstance(ctx.message.channel, discord.DMChannel)
        return commands.check(predicate)

    async def is_vegan(self, user: discord.User):
        guilds = [
            await self.bot.fetch_guild(self.config['discord_testing_guild_id']),
            await self.bot.fetch_guild(730907954345279591)
        ]
        for guild in guilds:
            vegan_role = get(guild.roles, name="Vegan")
            if vegan_role in user.roles:
                return True
        return False

    async def is_vegan_check(self, ctx):
        return await self.is_vegan(ctx.author)

    @tasks.loop(minutes=1)
    async def daily_loop(self):
        await self.bot.wait_until_ready()
        est_tz = pytz.timezone('US/Eastern')
        now_est = datetime.datetime.now(est_tz)
        print(now_est)
        if now_est.hour == 22 and now_est.minute == 0:
            guild_channels_map = {}
            for channel_id, guild_id in self.channel_guild_map.items():
                guild_channels_map.setdefault(guild_id, []).append(channel_id)
            for guild_id, channel_ids in guild_channels_map.items():
                loop_tags = await self.tag_manager.list_tags(
                    location_id=guild_id,
                    tag_type='loop'
                )
                loop_tags = [
                    t for t in loop_tags
                    if t.get('content') or t.get('attachment_url')
                ]
                if not loop_tags:
                    for cid in channel_ids:
                        channel = self.bot.get_channel(cid)
                        if channel:
                            await channel.send('No loop tags found for this guild.')
                    continue
                current_index = self.guild_loops_index[guild_id]
                for cid in channel_ids:
                    channel = self.bot.get_channel(cid)
                    if channel:
                        tag = loop_tags[current_index % len(loop_tags)]
                        msg = tag.get('content') or tag.get('attachment_url')
                        if msg:
                            await channel.send(msg)
                self.guild_loops_index[guild_id] = (current_index + 1) % len(loop_tags)

    @tasks.loop(hours=24)  # Adjust the interval as needed
    async def backup_task(self):
        try:
            backup_dir = setup_backup_directory('./backups')
            backup_file = perform_backup(
                db_user='postgres',
                db_name='lucy',
                db_host='localhost',
                backup_dir=backup_dir
            )

            logger.info(f'Backup completed successfully: {backup_file}')
        except Exception as e:
            logger.error(f'Error during database backup: {e}')

    @backup_task.before_loop
    async def before_backup(self):
        await self.bot.wait_until_ready()
>>>>>>> f442ff2 (93rd commit)

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if await self.predicator.is_at_home_func(before.guild.id):
            if before.content != after.content:
                ctx = await self.bot.get_context(after)
                if ctx.command:
                    await self.bot.invoke(ctx)

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
            ctx = await self.bot.get_context(message)
            if not array:
                logger.error("Invalid 'messages': The array is empty or improperly formatted.")
<<<<<<< HEAD
                if await self.predicator.is_at_home_func(message.guild.id):
                    print("Your message must include text or valid attachments.")
=======
                if await self.at_home().predicate(ctx):
                    await message.reply("Your message must include text or valid attachments.")
>>>>>>> f442ff2 (93rd commit)
                return
            logger.info(f"Final payload for processing: {json.dumps(array, indent=2)}")
            for item in array:
                # Moderation
                if self.config['openai_chat_moderation']:
                    async for moderation_completion in create_moderation(input_array=[item]):
                        try:
                            full_response = json.loads(moderation_completion)
                            results = full_response.get('results', [])
<<<<<<< HEAD
                            if results and results[0].get('flagged', False) and not self.predicator.is_spawd(ctx):
                                await self.handle_moderation(message)
                                return
                            elif self.predicator.is_spawd(ctx) and self.bot.user in message.mentions:
                                async for chat_completion in self.handler.generate_chat_completion(
                                    custom_id=message.author.id, array=[item], model='o1-mini'
                                ):
                                    if len(chat_completion) > 2000:
                                        unique_filename = f'temp_{uuid.uuid4()}.txt'
                                        with open('part_1_' + unique_filename, 'w') as f:
                                            f.write(chat_completion[:1000])
                                        await message.reply(file=discord.File('part_1_' + unique_filename))
                                        os.remove('part_1_' + unique_filename)
                                        with open('part_2_' + unique_filename, 'w') as f:
                                            f.write(chat_completion[1000:])
                                        await message.reply(file=discord.File('part_2_' + unique_filename))
                                        os.remove('part_2_' + unique_filename)
                                    else:
                                        await message.reply(chat_completion)
                            else:
                                async for moderation_completion in self.handler.generate_moderation_completion(custom_id=message.author.id, array=array):
                                    chat_moderation = json.loads(moderation_completion)
                                    academic_dishonesty_results = chat_moderation.get('results', [])
                                    academic_dishonesty_flagged = academic_dishonesty_results[0]['categories'].get('academic-dishonesty', False)
                                    if academic_dishonesty_flagged:
                                        academic_dishonesty_score = academic_dishonesty_results[0]['category_scores'].get('academic-dishonesty', 0)
                                        await self.handle_moderation(message)
                                        NLPUtils.append_to_jsonl(PATH_TRAINING, academic_dishonesty_score, message.content, message.author.id)
=======
                            if results and results[0].get('flagged', False):
                                if await self.at_home().predicate(ctx):
                                    await message.reply(
                                        f"Your file '{item.get('filename', 'unknown')}' was flagged for moderation."
                                    )
                                    await message.delete()
                                    return
                            else:
                                async for moderation_completion in self.handler.generate_moderation_completion(custom_id=message.author.id, array=array):
                                    carnism_flagged = results[0]['categories'].get('carnism', False)
                                    if carnism_flagged:
                                        carnism_score = results[0]['category_scores'].get('carnism', 0)
                                        if await self.at_home().predicate(ctx):
                                            await message.reply(
                                                f"Your file '{item.get('filename', 'unknown')}' was flagged for moderation."
                                            )
                                            await message.delete()
                                        NLPUtils.append_to_jsonl(PATH_TRAINING, carnism_score, message.content, message.author.id)
>>>>>>> f442ff2 (93rd commit)
                                        return
                                    if True:
                                        if self.config['openai_chat_completion'] and self.bot.user in message.mentions:
                                            async for chat_completion in self.handler.generate_chat_completion(
                                                custom_id=message.author.id, array=[item], sys_input=OPENAI_CHAT_SYS_INPUT
                                            ):
                                                if len(chat_completion) > 2000:
                                                    unique_filename = f'temp_{uuid.uuid4()}.txt'
                                                    with open(unique_filename, 'w') as f:
                                                        f.write(chat_completion)
                                                    await message.reply(file=discord.File(unique_filename))
                                                    os.remove(unique_filename)
                                                else:
                                                    await message.reply(chat_completion)
                        except Exception as e:
                            logger.error(traceback.format_exc())
<<<<<<< HEAD
                            if await self.predicator.is_at_home_func(message.guild.id):
                                print(f'An error occurred: {e}')
        except Exception as e:
            logger.error(traceback.format_exc())
            if await self.predicator.is_at_home_func(message.guild.id):
                print(f'An error occurred: {e}')
=======
                            if await self.at_home().predicate(ctx):
                                await message.reply(f'An error occurred: {e}')
                # Chat completion
                if await self.is_vegan(message.author) and message.guild.id is not 730907954345279591:
                    if self.config['openai_chat_completion'] and self.bot.user in message.mentions:
                        async for chat_completion in self.handler.generate_chat_completion(
                            custom_id=message.author.id, array=[item]
                        ):
                            await message.reply(chat_completion)
        except Exception as e:
            logger.error(traceback.format_exc())
            if await self.at_home().predicate(ctx):
                await message.reply(f'An error occurred: {e}')
>>>>>>> f442ff2 (93rd commit)
        finally:
            try:
                shutil.rmtree(DIR_TEMP)
                os.makedirs(DIR_TEMP, exist_ok=True)
                logger.info("Temporary files cleaned up successfully.")
            except Exception as cleanup_error:
                logger.error(f"Error cleaning up temporary files: {cleanup_error}")

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
