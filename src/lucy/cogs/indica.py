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
from lucy.utils.backup import perform_backup, setup_backup_directory
from lucy.utils.create_https_completion import Conversations
from lucy.utils.create_https_moderation import create_https_moderation
from lucy.utils.create_moderation import create_moderation
from lucy.utils.load_contents import load_contents
from lucy.utils.message import Message
from lucy.utils.nlp_utils import NLPUtils
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
from lucy.utils.helpers import *

class TagMenu(menus.ListPageSource):
    def __init__(self, tags):
        super().__init__(tags, per_page=1)  # One tag per page

    async def format_page(self, menu, tag):
        content = tag.get('content', '')
        attachment_url = tag.get('attachment_url', '')
        description = content or attachment_url or 'No content available.'
        embed = discord.Embed(
            title=f'Loop Tag: {tag['name']}',
            description=description,
            color=discord.Color.blurple()
        )
        return embed

class Indica(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.conversations = Conversations()
        self.handler = Message(self.config, self.conversations)
        self.tag_manager = TagManager(self.bot.db_pool)
        self.daily_loop.start()
        self.channel_guild_map: Dict[int, int] = {
            787738272616808509: 730907954345279591,
            730907954877956179: 730907954345279591,
            1315735859848544378: 1300517536001036348,
        }
        self.guild_loops_index = defaultdict(int)

    def at_home(self):
        async def predicate(ctx):
            return ctx.guild is not None and ctx.guild.id == self.bot.config['discord_testing_guild_id']
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
        if now_est.hour == 10 and now_est.minute == 0:
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
                for cid in channel_ids:
                    channel = self.bot.get_channel(cid)
                    if channel:
                        current_index = self.channel_loops_index[cid]
                        tag = loop_tags[current_index % len(loop_tags)]
                        msg = tag.get('content') or tag.get('attachment_url')
                        if msg:
                            await channel.send(msg)
                        self.channel_loops_index[cid] = (current_index + 1) % len(loop_tags)

    @tasks.loop(hours=24)
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

    @commands.Cog.listener()
    @commands.check(release_mode)
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)

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
                if await self.at_home().predicate(ctx):
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
                                if await self.at_home().predicate(ctx):
                                    if not await self.is_vegan(message.author):
                                        await message.reply(
                                            f"Your file '{item.get('filename', 'unknown')}' was flagged for moderation."
                                        )
                                        await message.delete()
                                    return
                            else:
                                async for moderation_completion in self.handler.generate_moderation_completion(custom_id=message.author.id, array=array):
                                    results = moderation_completion.get('results', [])
                                    carnism_flagged = results[0]['categories'].get('carnism', False)
                                    if carnism_flagged:
                                        carnism_score = results[0]['category_scores'].get('carnism', 0)
                                        if await self.at_home().predicate(ctx):
                                            if not await self.is_vegan(message.author):
                                                await message.reply(
                                                    f"Your file '{item.get('filename', 'unknown')}' was flagged for moderation."
                                                )
                                                await message.delete()
                                        NLPUtils.append_to_jsonl(PATH_TRAINING, carnism_score, message.content, message.author.id)
                                        return
                        except Exception as e:
                            logger.error(traceback.format_exc())
                            if await self.at_home().predicate(ctx):
                                await message.reply(f'An error occurred: {e}')
                # Chat completion
                if await self.is_vegan(message.author) and message.guild.id != 730907954345279591:
                    if self.config['openai_chat_completion'] and self.bot.user in message.mentions:
                        async for chat_completion in self.handler.generate_chat_completion(
                            custom_id=message.author.id, array=[item]
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
            if await self.at_home().predicate(ctx):
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
