from datetime import datetime
from lucy.utils.inc.setup_logging import logger

import aiohttp
import asyncio
import logging

class LinkedInBot:
    def __init__(self, config, db_pool, conversations, lock, oauth_token):
        self.config = config
        self.db_pool = db_pool
        self.conversations = conversations
        self.lock = lock
        self.oauth_token = oauth_token
        self.base_url = "https://api.linkedin.com/v2"

    async def event_ready(self):
        logger.info("LinkedIn Bot is ready!")
        logger.info("Listening for LinkedIn messages and events...")

    async def run(self):
        try:
            await self.event_ready()
        except Exception as e:
            logger.error(f"Error in LinkedInBot run loop: {e}")
