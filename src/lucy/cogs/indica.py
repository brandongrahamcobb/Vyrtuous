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
            798967615636504657: 730907954345279591,
            730907954877956179: 730907954345279591,
            1315735859848544378: 1300517536001036348,
        }
        self.guild_loops_index = defaultdict(int)

    @classmethod
    def at_home(cls):
        async def predicate(ctx):
            return ctx.guild is not None and ctx.guild.id == cls.bot.testing_guild_id
        return commands.check(predicate)

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

    @commands.Cog.listener()
    @commands.check(at_home)
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
            if message.attachments:
                array = await self.handler.process_array(message.content, attachments=message.attachments)
            else:
                array = await self.handler.process_array(message.content)
#            if isinstance(message.channel, discord.DMChannel):
 #               await ai_message()
            if self.config['openai_chat_completion']:
                if self.bot.user in message.mentions:
                    guilds = [
                        await self.bot.fetch_guild(self.config['discord_testing_guild_id']),
                        await self.bot.fetch_guild(730907954345279591)
                    ]

               # Get roles by name
                    vegan_roles = [
                        get(guilds[0].roles, name="Vegan"),
                        get(guilds[1].roles, name="Vegan")
                    ]
                    if (vegan_roles[0] in message.author.roles or vegan_roles[1] in message.author.roles):
                        if (message.guild.id == self.config['discord_testing_guild_id'] or message.channel.id == 985926652041261117):
                            async for chat_completion in self.handler.generate_chat_completion(custom_id=message.author.id, array=array): #, sys_input=OPENAI_CHAT_SYS_INPUT):
                                if len(chat_completion) > 2000:
                                    with open(PATH_COMPLETION, "w") as file:
                                        file.write(chat_completion)
                                    with open(PATH_COMPLETION, "rb") as file:
                                        await message.reply(f'Your response exceeded {self.config['discord_character_limit']} characters:', file=discord.File(file))
                                else:
                                    await message.reply(chat_completion)

            # Moderate Text and Images
            if self.config['openai_chat_moderation']:
                async for moderation_completion in create_moderation(input_array=array):
                    try:
                        full_response = json.loads(moderation_completion)
                        results = full_response.get('results', [])
                        if results:
                            flagged = results[0].get('flagged', False)
                            if flagged:
                                guilds = [await self.bot.fetch_guild(self.config['discord_testing_guild_id']), await self.bot.fetch_guild(730907954345279591)]
                                vegan_roles = [guilds[0].get_role(self.config['discord_role_pass']), guilds[1].get_role(788114978020392982)]
                                if vegan_roles[0] in message.author.roles and flagged:
                                    await message.delete()
                            else:
                                async for moderation_completion in self.handler.generate_moderation_completion(custom_id=message.author.id, array=array):
                                    carnism_flagged = results[0]['categories'].get('carnism', False)
                                    if carnism_flagged:
                                        carnism_score = results[0]['category_scores'].get('carnism', 0)
                                        NLPUtils.append_to_jsonl(PATH_TRAINING, carnism_score, message.content, message.author.id)
                        else:
                            logger.info('No moderation results returned.')
                    except Exception as e:
                        logger.error(f'Error processing moderation response: {e}')
                        await message.reply('An error occurred while processing moderation.')
        except Exception as e:
            logger.error(traceback.format_exc())
            await message.reply(e)

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
