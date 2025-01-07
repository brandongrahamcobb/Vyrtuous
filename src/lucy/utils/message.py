import json

from .helpers import *
from .setup_logging import logger

class Message:
    def __init__(self, config, conversations):
        self.config = config
        self.conversations = conversations

    async def generate_chat_completion(
        self,
        custom_id,
        array,
        *,
        completions: int = Ellipsis,  # Use Ellipsis as a sentinel
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
        # Handle defaults, keeping None if explicitly set
        if completions is Ellipsis:
            completions = OPENAI_CHAT_N
        if max_tokens is Ellipsis:
            max_tokens = OPENAI_MODEL_CONTEXT_LIMITS[self.config['openai_chat_model']] - OPENAI_MODEL_OUTPUT_LIMITS[self.config['openai_chat_model']]
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

        # Continue with the main logic
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
        async for moderation_completion in self.conversations.create_https_completion(
            completions=OPENAI_CHAT_MODERATION_N,
            custom_id=custom_id,
            input_array=array,
            max_tokens=OPENAI_MODEL_CONTEXT_LIMITS[self.config['openai_chat_moderation_model']] - OPENAI_MODEL_OUTPUT_LIMITS[self.config['openai_chat_moderation_model']],
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


    async def process_array(self, content, *, attachments = None):
        array = await self.process_text_message(content)
        if attachments:
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

