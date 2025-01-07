# bots/twitch_bot.py

import logging
from twitchio.ext import commands
from .create_https_completion import Conversations
from .message import Message
import json

from .setup_logging import logger

class TwitchBot(commands.Bot):
    def __init__(self, config, db_pool, conversations, lock, oauth_token):
        super().__init__(
            token=oauth_token,
            client_id=config['api_keys']['Twitch']['client_id'],
            nick='Lucy_',
            prefix="!",
            initial_channels=config.get('twitch_initial_channels', ['spawdspawd'])
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
        if message.echo:
            return
        logger.info(f"Received Twitch message: {message.content}")
        array = await self.handler.process_array(message.content)

        # Chat Completion
        if self.config.get('openai_chat_completion', False):
            async for response in self.handler.generate_moderation_completion(
                custom_id=message.author.id, array=array
            ):
                await message.reply(response)

        # Moderate Text and Images
        if self.config.get('openai_chat_moderation', False):
            async for moderation_completion in self.handler.generate_moderation_completion(
                custom_id=message.author.id, array=array
            ):
                full_response = json.loads(moderation_completion)
                results = full_response.get('results', [])
                flagged = results[0].get('flagged', False)
                carnism_flagged = results[0]['categories'].get('carnism', False)
                carnism_score = results[0]['category_scores'].get('carnism', 0)
                total_carnism_score = sum(arg['category_scores'].get('carnism', 0) for arg in results)
                if carnism_flagged or flagged:
                    await message.delete()
                    # Assuming NLPUtils is defined elsewhere
                    NLPUtils.append_to_other_jsonl('training.jsonl', carnism_score, message.content, message.author.id)
