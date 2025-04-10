# bots/twitch_bot.py
from lucy.utils.inc.setup_logging import logger
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

    async def event_message(self, message):
        logger.info(f'Received message: {message.content}')
