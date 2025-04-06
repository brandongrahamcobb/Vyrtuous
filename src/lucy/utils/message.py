''' message.py  The purpose of this program is to handle messages to OpenAI.
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
import discord
from discord.ext import commands
import aiofiles
import base64
import json
import os
import shutil
import tiktoken
from lucy.utils.predicator import Predicator
from .helpers import *
from .setup_logging import logger
from os.path import exists
import yaml 
from lucy.utils.create_moderation import create_moderation
import uuid
os.makedirs(DIR_TEMP, exist_ok=True)

class Message:
    def __init__(self, bot, config, conversations, db_pool):
        self.bot = bot
        self.config = config
        self.conversations = conversations
        self.db_pool = db_pool
        self.user_messages = {}
        self.predicator = Predicator(self.bot)

    async def generate_chat_completion(
        self,
        custom_id,
        array,
        *,
        completions: int = Ellipsis,
        max_tokens: int = Ellipsis,
        model: str = Ellipsis,
        response_format: str = Ellipsis,
        stop: str = Ellipsis,
        store: bool = Ellipsis,
        stream: bool = Ellipsis,
        sys_input: str = Ellipsis,
        temperature: float = Ellipsis,
        top_p: int = Ellipsis,
        use_history: bool = Ellipsis,
        add_completion_to_history: bool = Ellipsis
    ):
        if completions is Ellipsis:
            completions = OPENAI_CHAT_N
        encoder = tiktoken.get_encoding('cl100k_base')  # Get encoding for the selected model
        total_input_tokens = sum(
            [len(encoder.encode(message.get('text', ''))) for message in array if 'text' in message]
        )
        if model is Ellipsis:
            model = self.config['openai_chat_model']
        if any(message.get("type") == "image_url" for message in array):
            model = "gpt-4o-mini"
        for message in array:
            if 'text' not in message:
                logger.warning(f'Missing "text" in message: {message}')
        available_tokens = OPENAI_MODEL_CONTEXT_LIMITS[model] - total_input_tokens
        max_tokens = min(available_tokens, OPENAI_MODEL_OUTPUT_LIMITS[model])
        total_tokens = sum([len(message.get('text', '').split()) for message in array if 'text' in message])
        if response_format is Ellipsis:
            response_format = OPENAI_CHAT_RESPONSE_FORMAT
        if stop is Ellipsis:
            stop = self.config['openai_chat_stop']
        if store is Ellipsis:
            store = self.config['openai_chat_store']
        if stream is Ellipsis:
            stream = self.config['openai_chat_stream']
        if sys_input is Ellipsis:
            sys_input = self.config['openai_chat_sys_input']
        if temperature is Ellipsis:
            temperature = self.config['openai_chat_temperature']
        if top_p is Ellipsis:
            top_p = self.config['openai_chat_top_p']
        if use_history is Ellipsis:
            use_history = self.config['openai_chat_use_history']
        if add_completion_to_history is Ellipsis:
            add_completion_to_history = self.config['openai_chat_add_completion_to_history']
        while total_tokens + max_tokens > OPENAI_MODEL_CONTEXT_LIMITS[self.config['openai_chat_model']]:
            removed_message = array.pop(0)
            total_tokens -= len(removed_message.get('text', '').split())
        print(model)
        async for chat_completion in self.conversations.create_https_completion(
            custom_id=custom_id,
            completions=completions,
            input_array=array,
            max_tokens=max_tokens,
            model=model,
            response_format=response_format,
            stop=stop,
            store=store,
            stream=stream,
            sys_input=sys_input,
            temperature=temperature,
            top_p=top_p,
            use_history=use_history,
            add_completion_to_history=add_completion_to_history
        ):
            yield chat_completion

    async def generate_moderation_completion(self, custom_id, array):
#        encoder = tiktoken.get_encoding(f'openai_public:{self.config['openai_chat_model']}')  # Get encoding for the selected model
        encoder = tiktoken.get_encoding('cl100k_base')  # Get encoding for the selected model
        total_input_tokens = sum(
            [len(encoder.encode(message.get('text', ''))) for message in array if 'text' in message]  # Use .get() to avoid KeyError
        )
        for message in array:
            if 'text' not in message and message.get('type') == 'text':
               logger.warning(f'Missing "text" in message: {message}')
        available_tokens = OPENAI_MODEL_CONTEXT_LIMITS[OPENAI_CHAT_MODERATION_MODEL] - total_input_tokens
        max_tokens = min(available_tokens, OPENAI_MODEL_OUTPUT_LIMITS[self.config['openai_chat_model']])
        async for moderation_completion in self.conversations.create_https_completion(
            completions=OPENAI_CHAT_MODERATION_N,
            custom_id=custom_id,
            input_array=array,
            max_tokens=max_tokens,
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
            yield moderation_completion


    async def process_array(self, content, *, attachments=None):
        array = []
        if content.strip():
            array = await self.process_text_message(content)
        if attachments:
            processed_attachments = await self.process_attachments(attachments)
            array.extend(processed_attachments)
        if not array:
            logger.warning("Generated 'array' is empty. No valid text or attachments found.")
        return array

    async def process_attachments(self, attachments):
        processed_attachments = []
        for attachment in attachments:
            file_path = os.path.join(DIR_TEMP, attachment.filename)
            try:
                async with aiofiles.open(file_path, 'wb') as f:
                    await f.write(await attachment.read())
                if attachment.content_type.startswith('image/'):
                    async with aiofiles.open(file_path, 'rb') as f:
                        image_base64 = base64.b64encode(await f.read()).decode('utf-8')
                    mime_type = attachment.content_type
                    processed_attachments.append({
                        "type": "image_url",
                        "image_url": {
                            "url": f"data:{mime_type};base64,{image_base64}"
                        }
                    })
                elif attachment.content_type.startswith('text/'):
                    async with aiofiles.open(file_path, 'r', encoding='utf-8') as f:
                        text_content = await f.read()
                    processed_attachments.append({
                        "type": "text",
                        "text": text_content
                    })
                else:
                    logger.info(f"Skipping unsupported attachment type: {attachment.content_type}")
            except Exception as e:
                logger.error(f"Error processing file {attachment.filename}: {e}")
                continue
        return processed_attachments

    async def process_text_message(self, content):
        return [{
            'type': 'text',
            'text': content.replace(f'<@1318597210119864385>', '')
        }]

    def validate_array(self, array):
        valid = True
        for item in array:
            if item.get('type') == 'image_base64':
                if not item.get('image_data') or not item.get('content_type'):
                    logger.error(f"Invalid Base64 image data: {item}")
                    valid = False
            elif item.get('type') == 'text':
                if not item.get('text', '').strip():
                    logger.error(f"Invalid text content: {item}")
                    valid = False
        return valid

    async def send_dm(self, member: discord.Member, *, content: str = None, file: discord.File = None, embed: discord.Embed = None):
        channel = await member.create_dm()
        await self._send_message(channel.send, content=content, file=file, embed=embed)

    async def send_message(self, ctx: commands.Context, *, content: str = None, file: discord.File = None, embed: discord.Embed = None):
        can_send = (
            ctx.guild
            and isinstance(ctx.channel, discord.abc.GuildChannel)
            and ctx.channel.permissions_for(ctx.guild.me).send_messages
        )
        async with ctx.typing():
            if can_send:
                try:
                    await self._send_message(ctx.reply, content=content, file=file, embed=embed)
                except discord.HTTPException as e:
                    if e.code == 50035:  # Invalid Form Body due to message_reference
                        await self._send_message(ctx.send, content=content, file=file, embed=embed)
                    else:
                        raise
            else:
                await self.send_dm(ctx.author, content=content, file=file, embed=embed)

    async def _send_message(self, send_func, *, content: str = None, file: discord.File = None, embed: discord.Embed = None):
        kwargs = {}
        if content:
            kwargs["content"] = content
        if file:
            kwargs["file"] = file
        if embed:
            kwargs["embed"] = embed
        await send_func(**kwargs)


    async def ai_handler(self, ctx: commands.Context):
        array = await self.process_array(ctx.message.content, attachments=ctx.message.attachments)
        if not array:
            logger.error("Invalid 'messages': The array is empty or improperly formatted.")
            await self.send_message(ctx, content="Your message must include text or valid attachments.")
            return
        logger.info(f"Final payload for processing: {json.dumps(array, indent=2)}")
        for item in array:
            if self.config['openai_chat_moderation']:
                async for moderation_completion in create_moderation(input_array=[item]):
                    try:
                        full_response = json.loads(moderation_completion)
                        results = full_response.get('results', [])
                        if results and results[0].get('flagged', False) and not self.predicator.is_spawd(ctx):
                            result = results[0]
                            flagged = result.get('flagged', False)
                            categories = result.get('categories', {})
                            reasons = [category.replace("/", " → ").replace("-", " ").capitalize() for category, value in categories.items() if value is True]
                            if reasons:
                                reason_str = ", ".join(reasons)
                            else:
                                reason_str = "Unspecified moderation issue"
                            await self.handle_moderation(ctx.message, reason_str)
                            return
                    except Exception as e:
                        logger.error(traceback.format_exc())
                        print(f'An error occurred: {e}')
        if self.config['openai_chat_completion'] and self.bot.user in ctx.message.mentions:
            async for chat_completion in self.generate_chat_completion(
                custom_id=ctx.author.id, array=array, sys_input=OPENAI_CHAT_SYS_INPUT,
            ):
                await self.handle_large_response(ctx, chat_completion)

    def handle_users(self, author: str):
        author_char = author[0].upper()
        data = {letter: [] for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}
        users_file = join(DIR_HOME, '.users', 'users')
        if exists(users_file):
            with open(users_file, 'r+') as file:
                try:
                    data = yaml.safe_load(file) or data
                except yaml.YAMLError:
                    data = {letter: [] for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}
                if author_char not in data:
                    data[author_char] = []
                if author not in data[author_char]:
                    data[author_char].append(author)
                    file.seek(0)
                    file.write(yaml.dump(data))
                    file.truncate()
        else:
            with open(users_file, 'w') as file:
                yaml.dump(data, file)

    async def handle_large_response(self, ctx: commands.Context, response: str):
        if len(response) > 2000:
            unique_filename = f'temp_{uuid.uuid4()}.txt'
            try:
                with open(unique_filename, 'w') as f:
                    f.write(response)
                await self.send_file(ctx, file=discord.File(unique_filename))
            finally:
                os.remove(unique_filename)
        else:
            await self.send_message(ctx, content=response)

    async def handle_moderation(self, message: discord.Message, reason_str: str = "Unspecified moderation issue"):
        unfiltered_role = get(message.guild.roles, name=DISCORD_ROLE_PASS)
        if unfiltered_role in message.author.roles:
            return
        user_id = message.author.id
        async with self.db_pool.acquire() as connection:
            async with connection.transaction():
                row = await connection.fetchrow(
                    "SELECT flagged_count FROM moderation_counts WHERE user_id = $1", user_id
                )
                if row:
                    flagged_count = row["flagged_count"] + 1
                    await connection.execute(
                        "UPDATE moderation_counts SET flagged_count = $1 WHERE user_id = $2",
                        flagged_count, user_id
                    )
                    logger.info(f"Updated flagged count for user {user_id}: {flagged_count}")
                else:
                    flagged_count = 1
                    await connection.execute(
                        "INSERT INTO moderation_counts (user_id, flagged_count) VALUES ($1, $2)",
                        user_id, flagged_count
                    )
                    logger.info(f"Inserted new flagged count for user {user_id}: {flagged_count}")
        await message.delete()
        if flagged_count == 1:
            await self.send_message(
                message.author,
                content=f"{self.config['discord_moderation_warning']}. Your message was flagged for: {reason_str}"
            )
        elif flagged_count in [2, 3, 4]:
            if flagged_count == 4:
                await self.send_message(
                    message.author,
                    content=f"{self.config['discord_moderation_warning']}. Your message was flagged for: {reason_str}"
                )
        elif flagged_count >= 5:
            await self.send_message(
                message.author,
                content=f"{self.config['discord_moderation_warning']}. Your message was flagged for: {reason_str}"
            )
            await self.send_message(
                message.author,
                content="You have been timed out for 5 minutes due to repeated violations."
            )
            await message.author.timeout(datetime.timedelta(seconds=300))
            async with self.db_pool.acquire() as connection:
                await connection.execute(
                    "UPDATE moderation_counts SET flagged_count = 0 WHERE user_id = $1", user_id
                )

    def _handle_large_response(self, content: str) -> str:
        if len(content) > 2000:
            buffer = io.StringIO(content)
            file = discord.File(fp=buffer, filename="output.txt")
            return file
        return content

