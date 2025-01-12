# bots/twitch_bot.py
from lucy.utils.create_https_completion import Conversations
from lucy.utils.create_moderation import create_moderation
from lucy.utils.message import Message
from lucy.utils.nlp_utils import NLPUtils
from lucy.utils.setup_logging import logger
from twitchio.ext import commands
import logging
import json
import traceback  # Importing traceback for error handling

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
