# bots/twitch_bot.py
from .create_https_completion import Conversations
from .create_moderation import create_moderation
from .message import Message
from .nlp_utils import NLPUtils
from .setup_logging import logger
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
        try:
            if message.author.bot:
                return
            
            array = await self.handler.process_array(
                message.content, attachments=message.attachments
            )
            if not array:
                logger.error("Invalid 'messages': The array is empty or improperly formatted.")
                await message.channel.send("Your message must include text or valid attachments.")
                return
            
            logger.info(f"Final payload for processing: {json.dumps(array, indent=2)}")
            for item in array:
                # Moderation
                if self.config.get('openai_chat_moderation', False):
                    moderation_completion = await create_moderation(input_array=[item])  # Ensure create_moderation is defined
                    try:
                        full_response = json.loads(moderation_completion)
                        results = full_response.get('results', [])
                        if results and results[0].get('flagged', False):
                            await message.channel.send(
                                f"Your file '{item.get('filename', 'unknown')}' was flagged for moderation."
                            )
                            await message.delete()
                        else:
                            moderation_results = await self.handler.generate_moderation_completion(custom_id=message.author.id, array=array)
                            full_response = json.loads(moderation_results)
                            results = full_response.get('results', [])
                            carnism_flagged = results[0]['categories'].get('carnism', False)
                            if carnism_flagged:
                                carnism_score = results[0]['category_scores'].get('carnism', 0)
                                await message.channel.send(
                                    f"Your file '{item.get('filename', 'unknown')}' was flagged for moderation."
                                )
                                await message.delete()
                                NLPUtils.append_to_jsonl(PATH_TRAINING, carnism_score, message.content, message.author.id)
                                return
                    except Exception as e:
                        logger.error(traceback.format_exc())
                        await message.channel.send(f'An error occurred during moderation: {e}')
                
                # Chat completion
                if await self.is_donor(message.author):
                    if self.config.get('openai_chat_completion', False) and self.nick in message.mentions:
                        chat_completion = await self.handler.generate_chat_completion(
                            custom_id=message.author.id, array=[item]
                        )
                        await message.channel.send(chat_completion)

        except Exception as e:
            logger.error(traceback.format_exc())
            await message.channel.send(f'An error occurred: {e}')
    
    async def is_donor(self, user):
        # Example implementation for checking donor status
        donor_list = self.config.get('donors', [])
        return user.name in donor_list

