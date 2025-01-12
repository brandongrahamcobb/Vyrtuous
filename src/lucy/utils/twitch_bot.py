# bots/twitch_bot.py
from lucy.utils.create_https_completion import Conversations
from lucy.utils.message import Message
from lucy.utils.nlp_utils import NLPUtils
from lucy.utils.setup_logging import logger
from twitchio.ext import commands

import logging
import json


class TwitchBot(commands.Bot):
    def __init__(self, config, db_pool, conversations, lock, oauth_token):
        super().__init__(
            token=oauth_token,
            client_id=config['api_keys']['Twitch']['client_id'],
            nick='Lucy_',
            prefix="!",
            initial_channels=['spawdspawd']
        )
        self.conversations = conversations
        self.config = config
        self.db_pool = db_pool
        self.lock = lock
        self.handler = Message(self.config, self.conversations)

    async def event_ready(self):
        logger.info(f"Twitch Bot is ready! Logged in as {self.nick}")
        logger.info(f"Connected to channels: {', '.join(self.initial_channels)}")

    async def event_message(self, message):
        logger.info(f'Received message: {message.content}')
        try:
            if message.author.bot:
                return
            array = await self.handler.process_array(
                message.content, attachments=message.attachments
            )
            if not array:
                logger.error("Invalid 'messages': The array is empty or improperly formatted.")
                await message.send("Your message must include text or valid attachments.")
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
                                await message.send(
                                    f"Your file '{item.get('filename', 'unknown')}' was flagged for moderation."
                                )
                                await message.delete()
                            else:
                                async for moderation_completion in self.handler.generate_moderation_completion(custom_id=message.author.id, array=array):
                                    full_response = json.loads(moderation_completion)
                                    results = full_response.get('results', [])
                                    carnism_flagged = results[0]['categories'].get('carnism', False)
                                    if carnism_flagged:
                                        carnism_score = results[0]['category_scores'].get('carnism', 0)
                                        await message.send(
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
                if await self.is_donor(message.author):
                    if self.config['openai_chat_completion'] and self.bot.user in message.mentions:
                        async for chat_completion in self.handler.generate_chat_completion(
                            custom_id=message.author.id, array=[item]
                        ):
                            await message.send(chat_completion)
        except Exception as e:
            logger.error(traceback.format_exc())
            await message.send(f'An error occurred: {e}')
