import json

from .helpers import *
from .setup_logging import logger

class Message:
    def __init__(self, config, conversations):
        self.config = config
        self.conversations = conversations

    async def generate_chat_completion(self, custom_id, array):
        async for chat_completion in self.conversations.create_https_completion(
            custom_id=custom_id,
            completions=self.config['openai_chat_n'],
            input_array=array,
            max_tokens=self.config['openai_chat_max_tokens'],
            model=self.config['openai_chat_model'],
            response_format=self.config['openai_chat_response_format'],
            stop=self.config['openai_chat_stop'],
            store=self.config['openai_chat_store'],
            stream=self.config['openai_chat_stream'],
            sys_input=self.config['openai_chat_sys_input'],
            temperature=self.config['openai_chat_temperature'],
            top_p=self.config['openai_chat_top_p'],
            use_history=self.config['openai_chat_use_history'],
            add_completion_to_history=self.config['openai_chat_add_completion_to_history']
        ):
            yield chat_completion

    async def generate_moderation_completion(self, custom_id, array):
        async for moderation_completion in self.conversations.create_https_completion(
            completions=OPENAI_CHAT_MODERATION_N,
            custom_id=custom_id,
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
            yield moderation_completion


    async def process_array(self, content, attachments):
        array = await self.process_text_message(content)
        array.extend(await self.process_attachments(attachments))
        return array

    async def process_attachments(self, attachments):
        array = []
        for attachment in attachments:
            if attachment.content_type and attachment.content_type.startswith('image/'):
                array.append({
                    'type': 'image_url',
                    'image_url': {'url': attachment.url}
                })
            elif attachment.content_type and attachment.content_type.startswith('text/'):
                try:
                    file_content = await attachment.read()
                    text_content = file_content.decode('utf-8')
                    array.append({'type': 'text', 'text': text_content})
                except Exception as e:
                    logger.error(f"Error reading attachment {attachment.filename}: {e}")
        return array

    async def process_text_message(self, content):
        return [{
            'type': 'text',
            'text': content.replace(f'<@1318597210119864385>', '')
        }]

