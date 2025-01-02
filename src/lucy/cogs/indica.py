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
        self.conversations = Conversations()
        self.config = bot.config
        self.lock = asyncio.Lock()
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
        try:
            if self.bot.user == message.author:
                return
            async with self.lock:
                # Input Image/Text
                array = []
#                result = NLPUtils.combined_analysis(message.content)
 #               if result['sentiment']['label'].lower() == 'negative':
  #                  pass
   #             else:
    #                NLPUtils.append_to_jsonl('training.jsonl', result['sentiment'], message.content)
                input_text_dict = {
                    'type': 'text',
                    'text': message.content.replace('<@1318597210119864385>', '')
                }
                array.append(input_text_dict)
                for attachment in message.attachments:
                    if attachment.content_type and attachment.content_type.startswith('image/'):
                        input_image_dict = {
                            'type': 'image_url',
                            'image_url': {
                                'url': attachment.url
                            }
                        }
                        array.append(input_image_dict)
                    elif attachment.content_type and attachment.content_type.startswith('text/'):
                    # Read the content of the text file
                        try:
                            file_content = await attachment.read()  # Read the file content
                            text_content = file_content.decode('utf-8')  # Decode bytes to string
                            input_text_dict = {
                                'type': 'text',
                                'text': text_content  # Add the text content to the array
                            }
                            array.append(input_text_dict)
                        except Exception as e:
                            logger.error(f'Error reading text attachment {attachment.filename}: {e}')
                            await message.channel.send(f"Could not read the content of the text file: {attachment.filename}")
                async for moderation in create_https_moderation(message.author.id, array, model=OPENAI_MODERATION_MODEL):
                    results = moderation.get('results', [])
                    if results and results[0].get('flagged', False):
                        await message.delete()
                        channel = await message.author.create_dm()
                        await channel.send(self.config['openai_moderation_warning'])
                        break
      
#                if message.channel.id == 1317987851593584651:
#                if message.attachments:
 #                   async for moderation in create_moderation(array):
  #                      results = moderation.get('results', [])
   #                     if results and results[0].get('flagged', False):
    #                        await message.delete()
     #                       channel = await message.author.create_dm()
      #                      await channel.send(self.config['openai_moderation_warning'])
       #                     break
                # Chat Completion
#                if self.bot.user in message.raw_mentions and not isinstance(message.type, MessageType.reply):
 #                   self.conversations.clear()
                if self.bot.user in message.mentions:
                    async for response in self.conversations.create_https_completion(
                        completions=self.config['openai_chat_n'],
                        custom_id=message.author.id,
                        input_array=array,
                        max_tokens=self.config['openai_chat_max_tokens'],
                        model=self.config['openai_chat_model'],
                        response_format=self.config['openai_chat_response_format'],
                        stop=self.config['openai_chat_stop'],
                        store=self.config['openai_chat_store'],
                        stream=self.config['openai_chat_stream'],
                        sys_input=self.sys_input, #self.cwonfig['openai_chat_sys_input'], # self.sys_input, 
                        temperature=self.config['openai_chat_temperature'],
                        top_p=self.config['openai_chat_top_p'],
                        use_history=self.config['openai_chat_use_history'],
                        add_completion_to_history=self.config['openai_chat_add_completion_to_history']
                    ):
                        await message.reply(response)
#                # Chat Moderation
                if self.config['openai_chat_moderation']:
                    role = message.guild.get_role(1308689505158565918)
                    async for moderation in self.conversations.create_https_completion(
                        completions=OPENAI_CHAT_MODERATION_N,
                        custom_id=message.author.id,
                        input_array=array,
                        max_tokens=OPENAI_CHAT_MODERATION_MAX_TOKENS,
                        model=OPENAI_CHAT_MODERATION_MODEL,
                        response_format=OPENAI_CHAT_MODERATION_RESPONSE_FORMAT,
                        stop=OPENAI_CHAT_MODERATION_STOP,
                        store=OPENAI_CHAT_MODERATION_STORE,
                        stream=OPENAI_CHAT_MODERATION_STREAM,
                        sys_input=OPENAI_CHAT_MODERATION_SYS_INPUT,
                        temperature=OPENAI_CHAT_MODERATION_TEMPERATURE,
                        top_p=OPENAI_CHAT_MODERATION_TOP_P,
                        use_history=OPENAI_CHAT_MODERATION_USE_HISTORY,
                        add_completion_to_history=OPENAI_CHAT_MODERATION_ADD_COMPLETION_TO_HISTORY
                    ):
                        full_response = json.loads(moderation)
                        results = full_response.get('results', [])
                        flagged = results[0].get('flagged', False)
                        carnism_flagged = results[0]['categories'].get('carnism', False)
                        carnism_score = results[0]['category_scores'].get('carnism', 0)
                        total_carnism_score = sum(arg['category_scores'].get('carnism', 0) for arg in results)
                        if carnism_flagged or flagged:  # If carnism is flagged
                            if role not in message.author.roles:
                                if not self.config['discord_role_pass']:
                                    channel = await message.author.create_dm()
                                    await channel.send(f'Your message: {message.content} was flagged for promoting carnism.')
                                await message.delete()
                            NLPUtils.append_to_other_jsonl('training.jsonl', carnism_score, message.content, message.author.id) #results[0].get('flagged', False), message.content)
                            break
#                except Exception as e:
 #                   print(e)
        except Exception as e:
            print(e)

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
