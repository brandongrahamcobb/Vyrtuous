
import aiofiles
import base64
import json
import tiktoken
import os
import shutil

from .helpers import *
from .setup_logging import logger

os.makedirs(DIR_TEMP, exist_ok=True)

class Message:
    def __init__(self, config, conversations):
        self.config = config
        self.conversations = conversations

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
        for message in array:
            if 'text' not in message:
                logger.warning(f'Missing "text" in message: {message}')
        available_tokens = OPENAI_MODEL_CONTEXT_LIMITS[self.config['openai_chat_model']] - total_input_tokens
        max_tokens = min(available_tokens, OPENAI_MODEL_OUTPUT_LIMITS[self.config['openai_chat_model']])
        total_tokens = sum([len(message.get('text', '').split()) for message in array if 'text' in message])
        if model is Ellipsis:
            model = self.config['openai_chat_model']
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
    
        # Add text content if present
        if content.strip():
            array = await self.process_text_message(content)
    
        # Process attachments if present
        if attachments:
            processed_attachments = await self.process_attachments(attachments)
            array.extend(processed_attachments)
    
        # Log a warning if the array is empty
        if not array:
            logger.warning("Generated 'array' is empty. No valid text or attachments found.")
    
        return array

    async def process_attachments(self, attachments):
        """
        Processes attachments:
        - Saves files in DIR_TEMP/temp/.
        - Converts files to Base64 or raw text for further processing.
        """
        processed_attachments = []
    
        for attachment in attachments:
            file_path = os.path.join(DIR_TEMP, attachment.filename)
            try:
                # Save the file locally
                async with aiofiles.open(file_path, 'wb') as f:
                    await f.write(await attachment.read())
    
                if attachment.content_type.startswith('image/'):
                    # Convert image to Base64
                    async with aiofiles.open(file_path, 'rb') as f:
                        image_base64 = base64.b64encode(await f.read()).decode('utf-8')
                    processed_attachments.append({
                        'type': 'image_base64',
                        'image_data': image_base64,
                        'filename': attachment.filename
                    })
    
                elif attachment.content_type.startswith('text/'):
                    # Read text content
                    async with aiofiles.open(file_path, 'r', encoding='utf-8') as f:
                        text_content = await f.read()
                    processed_attachments.append({
                        'type': 'text',
                        'text': text_content,
                        'filename': attachment.filename
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

