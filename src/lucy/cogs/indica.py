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
from os.path import abspath, dirname, exists, expanduser, join
from lucy.utils.create_https_completion import Conversations
from lucy.utils.create_https_moderation import create_https_moderation
from lucy.utils.create_moderation import create_moderation
from lucy.utils.nlp_utils import NLPUtils
from lucy.utils.load_contents import load_contents
from lucy.utils.message import Message

import asyncio
import datetime
import discord
import json
import os
import subprocess
import traceback
from lucy.utils.helpers import *

class Indica(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.conversations = Conversations()
#        self.hybrid = self.bot.get_cog('Hybrid')
 #       self.sativa = self.bot.get_cog('Sativa')
        self.hybrid = load_contents(PATH_HYBRID)
        self.indica = load_contents(PATH_INDICA)
        self.sativa = load_contents(PATH_SATIVA)
        self.add_watermark = load_contents(PATH_ADD_WATERMARK)
        self.adjust_hue_and_saturation = load_contents(PATH_ADJUST_HUE_AND_SATURATION)
        self.arpp = load_contents(PATH_ARPP)
        self.benchmark = load_contents(PATH_BENCHMARK)
        self.clear_screen = load_contents(PATH_CLEAR_SCREEN)
        self.combine = load_contents(PATH_COMBINE)
        self.create_batch_completion = load_contents(PATH_CREATE_BATCH_COMPLETION)
        self.create_https_completion = load_contents(PATH_CREATE_HTTPS_COMPLETION)
        self.create_moderation = load_contents(PATH_CREATE_MODERATION)
        self.discord = load_contents(PATH_DISCORD_UTILS)
        self.draw_fingerprint = load_contents(PATH_DRAW_FINGERPRINT)
        self.draw_watermarked_molecule = load_contents(PATH_DRAW_WATERMARKED_MOLECULE)
        self.fine_tuning = load_contents(PATH_FINE_TUNING)
        self.format_error_check = load_contents(PATH_FORMAT_ERROR_CHECK)
        self.get_molecule_name = load_contents(PATH_GET_MOLECULE_NAME)
        self.get_mol = load_contents(PATH_GET_MOL)
        self.get_proximity = load_contents(PATH_GET_PROXIMITY)
        self.google = load_contents(PATH_GOOGLE)
        self.gsrs = load_contents(PATH_GSRS)
        self.handler = Message(self.config, self.conversations)
        self.helpers = load_contents(PATH_HELPERS)
        self.increment_version = load_contents(PATH_INCREMENT_VERSION)
        self.load_contents = load_contents(PATH_LOAD_CONTENTS)
        self.load_yaml = load_contents(PATH_LOAD_YAML)
        self.prompt_for_values = load_contents(PATH_PROMPT_FOR_VALUES)
        self.script = load_contents(PATH_SCRIPT)
        self.setup_logging = load_contents(PATH_SETUP_LOGGING)
        self.tag = load_contents(PATH_TAG)
        self.unique_pairs = load_contents(PATH_UNIQUE_PAIRS)
        self.sum_of_paths = f'''
            {self.adjust_hue_and_saturation} and {self.arpp} and {self.benchmark} and {self.clear_screen} and {self.combine} and {self.create_batch_completion} and and {self.create_https_completion} and {self.create_moderation} and {self.discord} and {self.draw_fingerprint} and {self.draw_watermarked_molecule} and {self.fine_tuning} and {self.format_error_check} and {self.get_molecule_name} and {self.get_mol} and {self.get_proximity} and {self.google} and {self.gsrs} and {self.helpers} and {self.hybrid} and {self.increment_version} and {self.indica} and {self.load_contents} and {self.load_yaml} and {self.sativa} and {self.setup_logging} and {self.tag} and {self.unique_pairs}
        '''
        self.sys_input = f'''
            Your utilities are {self.sum_of_paths}.
        '''

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)

    @commands.Cog.listener()
    async def on_message(self, message):
        logger.info(f"Received message: {message.content}")
        if message.author == self.bot.user:
            return
        array = await self.handler.process_array(message.content, message.attachments)

        # Chat
        if self.config['openai_chat_completion']:
           async for chat_completion in self.handler.generate_chat_completion(custom_id=message.author.id, array=array):
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
                    if not self.config['discord_role_pass']:
                        await message.delete()
                    NLPUtils.append_to_other_jsonl('training.jsonl', carnism_score, message.content, message.author.id)

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
