
import json
import tiktoken

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
        array = await self.process_text_message(content)
        image_exceeded = False
        if attachments:
            processed_array, image_exceeded = await self.process_attachments(attachments)
            array.extend(processed_array)
        return array, image_exceeded

    async def process_attachments(self, attachments):
        """
        Processes attachments by:
        - Keeping the text content for the first attachment.
        - For subsequent image attachments, dropping the text content.
        - Yielding each attachment one by one.
        
        Yields:
            tuple: (processed_attachment, image_exceeded_flag)
        """
        image_count = 0  # Counter to track the number of images processed
        image_exceeded = False  # Flag to indicate if image limit was exceeded
        first_attachment = True  # Flag to check if processing the first attachment

        for attachment in attachments:
            processed_attachment = {}

            if attachment.content_type and attachment.content_type.startswith('image/'):
                image = {
                    'type': 'image_url',
                    'image_url': {'url': attachment.url}
                }

                if first_attachment:
                    # Keep any text associated with the first image if applicable
                    processed_attachment = image
                    # If there's associated text, include it
                    # Assuming attachment has 'description' or similar field
                    # Modify as per your actual data structure
                    if hasattr(attachment, 'description') and attachment.description:
                        processed_attachment['text'] = attachment.description
                    image_count += 1
                    first_attachment = False
                else:
                    # Subsequent images: Append without text
                    processed_attachment = image
                    image_count += 1

            elif attachment.content_type and attachment.content_type.startswith('text/'):
                if first_attachment:
                    try:
                        file_content = await attachment.read()
                        text_content = file_content.decode('utf-8')
                        processed_attachment = {
                            'type': 'text',
                            'text': text_content
                        }
                        first_attachment = False
                    except Exception as e:
                        logger.error(f"Error reading attachment {attachment.filename}: {e}")
                        continue  # Skip this attachment on error
                else:
                    # For subsequent text attachments, decide whether to include or skip
                    # Here, we'll skip additional text attachments
                    logger.info(f"Skipping additional text attachment: {attachment.filename}")
                    continue
            else:
                # Handle other types of attachments if necessary
                logger.info(f"Skipping unsupported attachment type: {attachment.content_type}")
                continue

            # Check if image limit is exceeded (example: more than 5 images)
            MAX_IMAGES = 5  # Define as per your requirements
            if image_count > MAX_IMAGES:
                image_exceeded = True
                logger.warning(f'More than {MAX_IMAGES} images detected in message.')

            # Yield the processed attachment and the image_exceeded flag
            if self.config.get('openai_chat_model') in OPENAI_CHAT_MODELS.get('current', ['o1-mini', 'o1-preview']):
                yield processed_attachment, image_exceeded

        # After processing all attachments, you can perform any final actions if needed
        # For example, notifying if image limit was exceeded
        # This part depends on how you intend to use the yielded data

    async def process_text_message(self, content):
        return [{
            'type': 'text',
            'text': content.replace(f'<@1318597210119864385>', '')
        }]

