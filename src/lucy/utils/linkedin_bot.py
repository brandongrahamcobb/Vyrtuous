import aiohttp
import asyncio
import logging
from datetime import datetime
from .setup_logging import logger
from .linkedin_oauth import LinkedInOAuth

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

    async def reply_to_message(self, recipient_id, message):
        """Send a reply to a LinkedIn message."""
        url = f"{self.base_url}/messages"
        headers = {
            "Authorization": f"Bearer {self.oauth_token}",
            "Content-Type": "application/json",
            "X-Restli-Protocol-Version": "2.0.0",
        }
        payload = {
            "recipients": [f"urn:li:person:{recipient_id}"],
            "message": {"text": message},
        }

        async with aiohttp.ClientSession() as session:
            async with session.post(url, headers=headers, json=payload) as response:
                if response.status == 201:
                    logger.info(f"Replied to {recipient_id}: {message}")
                else:
                    logger.error(f"Failed to send LinkedIn message. Status: {response.status}, Response: {await response.text()}")

    async def run(self):
        """Main loop to fetch and process LinkedIn messages."""
        await self.event_ready()
        while True:
            try:
                messages = await self.fetch_messages()
                if messages:
                    await self.process_messages(messages)
            except Exception as e:
                logger.error(f"Error in LinkedInBot run loop: {e}")
            await asyncio.sleep(30)  # Adjust polling interval as needed
