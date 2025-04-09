''' ai.py  The purpose of this program is to provide generic Natural Language Processing functionality.
    Copyright (C) 2024  github.com/brandongrahamcobb

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from lucy.utils.helpers import *
from lucy.utils.load_yaml import load_yaml
from lucy.utils.setup_logging import logger
from openai import AsyncOpenAI
import aiohttp
from typing import List, Optional, Any, Dict, Union

import aiohttp
import asyncio
import datetime
import json
import openai
import os
import traceback

config = load_yaml(PATH_CONFIG_YAML)
OPENAI_API_KEY = config['api_keys']['OpenAI']['api_key']

class NLPUtils:

    @staticmethod
    def append_to_jsonl(file_path, carnism, completion, custom_id):
        entry = {
            'messages': [
                {'role': 'system', 'content': 'You are Vyrtuous, an AI. Vyrtuous is programmed to generate a training file which contains non-vegan messages.'},
                {'role': 'user', 'content': json.dumps({'carnism': carnism})},
                {'role': 'assistant', 'content': completion}
            ],
            'metadata': {
                'user': str(custom_id),
                'timestamp': str(datetime.datetime.now(datetime.timezone.utc))
            }
        }
        try:
            with open(file_path, 'a') as file:
                json.dump(entry, file)
                file.write('\n')
        except Exception as e:
            logger.error(f'Error occurred while appending entry: {e}')
            return {'error': str(e)}

class BatchProcessor:

    def __init__(self, bot):
        self.bot = bot
        self.results = {}
        if os.path.exists(PATH_OPENAI_RESULTS):
            with open(PATH_OPENAI_RESULTS, 'r') as f:
                self.results = json.load(f)

    async def upload_file(self):
        url = 'https://api.openai.com/v1/files'
        headers = {
            'Authorization': f'Bearer {OPENAI_API_KEY}',
        }
        if not os.path.exists(PATH_OPENAI_REQUESTS) or os.stat(PATH_OPENAI_REQUESTS).st_size == 0:
            return 'No batch requests to process.'
        form_data = aiohttp.FormData()
        form_data.add_field(
            'file',
            open(PATH_OPENAI_REQUESTS, 'rb'),
            filename=os.path.basename(PATH_OPENAI_REQUESTS),
            content_type='application/jsonl'
        )
        form_data.add_field('purpose', 'batch')
        async with aiohttp.ClientSession() as session:
            async with session.post(url, headers=headers, data=form_data) as response:
                resp_json = await response.json()
                if response.status == 200:
                    return resp_json['id']  # Return the uploaded file ID
                else:
                    return f'File upload failed: {resp_json}'

    async def create_batch(self, input_file_id):
        url = 'https://api.openai.com/v1/batches'
        headers = {'Authorization': f'Bearer {OPENAI_API_KEY}', 'Content-Type': 'application/json'}
        data = {'input_file_id': input_file_id, 'endpoint': '/v1/chat/completions', 'completion_window': '24h'}
        async with aiohttp.ClientSession() as session:
            async with session.post(url, headers=headers, json=data) as response:
                return await response.json()

    async def retrieve_batch_results(self, batch_id):
        url = f'https://api.openai.com/v1/batches/{batch_id}'
        headers = {'Authorization': f'Bearer {OPENAI_API_KEY}'}
        async with aiohttp.ClientSession() as session:
            while True:
                async with session.get(url, headers=headers) as response:
                    batch_status = await response.json()
                    if batch_status['status'] == 'completed':
                        output_file_id = batch_status['output_file_id']
                        return await self.download_results(output_file_id)
                    elif batch_status['status'] in ['failed', 'cancelled']:
                        return f'Batch processing failed: {batch_status}'
                    print('Batch still processing... Checking again in 5 minutes.')
                    await asyncio.sleep(300)

    async def download_results(self, output_file_id):
        url = f'https://api.openai.com/v1/files/{output_file_id}/content'
        headers = {'Authorization': f'Bearer {OPENAI_API_KEY}'}
        async with aiohttp.ClientSession() as session:
            async with session.get(url, headers=headers) as response:
                results = await response.json()
                for entry in results:
                    user_id = entry['custom_id'].split('-')[0]
                    user = self.bot.get_user(int(user_id))
                    if not user:
                        continue
                    username = user.name
                    first_letter = username[0].upper()
                    if first_letter not in self.results:
                        self.results[first_letter] = {}
                    if username not in self.results[first_letter]:
                        self.results[first_letter][username] = []
                    self.results[first_letter][username].append(entry['response'])
                with open(PATH_OPENAI_RESULTS, 'w') as file:
                    json.dump(self.results, file, indent=4)
                return 'Batch processing completed! Results saved.'

    async def process_batches(self):
        if not os.path.exists(PATH_OPENAI_REQUESTS) or os.stat(PATH_OPENAI_REQUESTS).st_size == 0:
            return 'No batch requests to process.'
        file_id = await self.upload_file()
        batch_response = await self.create_batch(file_id)
        batch_id = batch_response['id']
        print(f'Batch created: {batch_id}. Waiting for results...')
        result_message = await self.retrieve_batch_results(batch_id)
        open(PATH_OPENAI_REQUESTS, 'w').close()
        return result_message

    def get_user_responses(self, user: discord.User):
        first_letter = user.name[0].upper()
        if first_letter in self.results and user.name in self.results[first_letter]:
            responses = self.results[first_letter][user.name]
            del self.results[first_letter][user.name]
            if not self.results[first_letter]:  # Remove empty letters
                del self.results[first_letter]
            with open(PATH_OPENAI_RESULTS, 'w') as file:
                json.dump(self.results, file, indent=4)
            return responses
        return None

class Conversations:
    def __init__(self):
        self.conversations = defaultdict(list)

    def trim_conversation_history(self, model, custom_id):
        max_context_length = OPENAI_MODEL_CONTEXT_LIMITS.get(model, 4096)
        total_tokens = sum(len(msg['content']) for msg in self.conversations[custom_id])
        while total_tokens > max_context_length:
            removed_message = self.conversations[custom_id].pop(0)
            total_tokens -= len(removed_message['content'])

    async def create_https_completion(
        self,
        completions,
        custom_id,
        input_array,
        max_tokens,
        model,
        response_format,
        stop,
        store,
        stream,
        sys_input,
        temperature,
        top_p,
        use_history=True,
        add_completion_to_history=True
    ):
        try:
            config = load_yaml(PATH_CONFIG_YAML)
            api_key = config['api_keys']['OpenAI']['api_key']
            ai_client = AsyncOpenAI(api_key=api_key)
            headers = {'Authorization': f'Bearer {api_key}'}
            logger.info('Headers prepared for the request.')
            messages = []
            if sys_input:
                messages.append({'role': 'system', 'content': sys_input})
                logger.info('System input added at the beginning.')
            for message in input_array:
                if 'text' in message and message['text']:
                    messages.append({'role': 'user', 'content': message['text']})
                elif 'files' in message:
                    for file in message['files']:
                        file_content = await self.process_file(file)
                        messages.append({'role': 'user', 'content': file_content})
            if use_history and custom_id in self.conversations:
                messages = self.conversations[custom_id] + messages
            self.conversations[custom_id].append({'role': 'user', 'content': input_array})
            self.trim_conversation_history(model, custom_id)
            request_data = {
                'messages': messages,
                'model': model,
                'temperature': float(temperature),
                'top_p': float(top_p),
                'n': int(completions),
                'stop': stop,
                'store': bool(store),
                'stream': bool(stream),
            }
            if response_format:
                request_data['response_format'] = response_format
            if model in {'chatgpt-4o-latest', 'o1-mini', 'o1-preview'}:
                request_data['max_completion_tokens'] = int(max_tokens)
                request_data['temperature'] = 1.0
            else:
                request_data['max_tokens'] = int(max_tokens)
            if store:
                request_data.update({
                    'metadata': {'user': str(custom_id), 'timestamp': str(datetime.datetime.now(datetime.timezone.utc))}
                })
            async with aiohttp.ClientSession() as session:
                try:
                    async with session.post(url=OPENAI_ENDPOINT_URLS['chat'], headers=headers, json=request_data) as response:
                        logger.info(f'Received response with status: {response.status}.')
                        full_response = ''
                        if stream:
                            if response.status != 200:
                                return
                            async for line in response.content:
                                decoded_line = line.decode('utf-8').strip()
                                if not decoded_line.startswith('data: ') or len(decoded_line) <= 6:
                                    continue
                                try:
                                    data_chunk = json.loads(decoded_line[6:])
                                    if 'choices' in data_chunk:
                                        for choice in data_chunk['choices']:
                                            content = choice['delta'].get('content', '')
                                            full_response += content
                                            if choice.get('finish_reason') == 'stop':
                                                logger.info('Completion streaming stopped.')
                                                break
                                except json.JSONDecodeError:
                                    continue
                        else:
                            full_response_json = await response.json()
                            full_response = full_response_json['choices'][0]['message']['content']
                        if add_completion_to_history:
                            self.conversations[custom_id].append({'role': 'assistant', 'content': full_response})
                        for chunk in self.split_long_response(full_response, DISCORD_CHARACTER_LIMIT):
                            yield chunk
                except Exception as e:
                    yield traceback.format_exc()
        except Exception as e:
            logger.error('Error in create_https_completion.', exc_info=True)
            yield traceback.format_exc()

    def split_long_response(self, response, limit):
        parts = response.split('```')
        output_chunks = []
        for i, part in enumerate(parts):
            if i % 2 == 0:
                non_code_chunks = [part[j:j + limit] for j in range(0, len(part), limit)]
                output_chunks.extend(non_code_chunks)
            else:
                output_chunks.append(f'```{part}```')
        return output_chunks

async def create_completion(input_array):
    try:
        config = load_yaml(PATH_CONFIG_YAML)
        api_key = config['api_keys']['OpenAI']['api_key']
        ai_client = AsyncOpenAI(api_key=api_key)
        response = await ai_client.chat.completions.create(
            model='gpt-4o-mini',
            messages=input_array,
            response_format=OPENAI_CHAT_COLORIZE_RESPONSE_FORMAT
        )
        yield response.choices[0].message.content
    except Exception as e:
        yield {'error': traceback.format_exc()}

async def create_https_moderation(custom_id, input_array, model):
    try:
        config = load_yaml(PATH_CONFIG_YAML)
        api_key = config['api_keys']['OpenAI']['api_key']
        ai_client = AsyncOpenAI(api_key=api_key)
        headers = {'Authorization': f'Bearer {api_key}'}
        request_data = {
            'input': input_array,
            'model': model,
        }
        async with aiohttp.ClientSession() as session:
            try:
                async with session.post(url=OPENAI_ENDPOINT_URLS['moderations'], headers=headers, json=request_data) as moderation_object:
                    if moderation_object.status == 200:
                        response_data = await moderation_object.json()
                        yield response_data
                    else:
                        error_message = await moderation_object.text()
                        yield {'error': error_message}
            except Exception as e:
                logger.error('An error occurred while making the HTTP request.')
                logger.error(traceback.format_exc())
                yield traceback.format_exc()
    except Exception as e:
        logger.error('An error occurred in create_https_moderation.')
        logger.error(traceback.format_exc())
        yield traceback.format_exc()

async def create_moderation(input_array):
    try:
        config = load_yaml(PATH_CONFIG_YAML)
        api_key = config['api_keys']['OpenAI']['api_key']
        ai_client = AsyncOpenAI(api_key=api_key)
        for item in input_array:
            response = await ai_client.moderations.create(
                model='omni-moderation-latest',
                input=input_array,
            )
            moderation_response = response.json()
            yield moderation_response
    except Exception as e:
        error_details = traceback.format_exc()
        logger.error(f'An error occurred during moderation: {error_details}')
        yield {'error': error_details}

@dataclass
class ModerationResult:
    object: str
    input_tokens: int
    num_model_requests: int
    project_id: Optional[str]
    user_id: Optional[str]
    api_key_id: Optional[str]
    model: Optional[str]

@dataclass
class ModerationBucket:
    object: str
    start_time: int
    end_time: int
    results: List[ModerationResult]

@dataclass
class CompletionResult:
    object: str
    input_tokens: int
    output_tokens: int
    input_cached_tokens: int
    input_audio_tokens: int
    output_audio_tokens: int
    num_model_requests: int
    project_id: Optional[str]
    user_id: Optional[str]
    api_key_id: Optional[str]
    model: Optional[str]
    batch: Optional[bool]

@dataclass
class CompletionBucket:
    object: str
    start_time: int
    end_time: int
    results: List[CompletionResult]

@dataclass
class UsagePage:
    object: str
    data: List[Union[ModerationBucket, CompletionBucket]]
    has_more: bool
    next_page: Optional[str]

class OpenAIUsageClient:
    BASE_URL = 'https://api.openai.com/v1/organization/usage'

    def __init__(self, api_key: str, organization_id: Optional[str] = None):
        self.api_key = api_key
        self.organization_id = organization_id
        self.headers = {
            'Authorization': f'Bearer {self.api_key}',
            'Content-Type': 'application/json'
        }

    async def _get_usage(self, endpoint: str, params: Dict[str, Any]) -> UsagePage:
        url = f'{self.BASE_URL}/{endpoint}'
        if self.organization_id:
            params['organization_id'] = self.organization_id

        async with aiohttp.ClientSession() as session:
            async with session.get(url, headers=self.headers, params=params) as response:
                if response.status != 200:
                    text = await response.text()
                    raise Exception(f'API request failed with status {response.status}: {text}')
                data = await response.json()
        usage_page = UsagePage(
            object=data.get('object'),
            data=[],
            has_more=data.get('has_more'),
            next_page=data.get('next_page')
        )
        for bucket in data.get('data', []):
            if endpoint == 'moderations':
                results = [ModerationResult(**result) for result in bucket.get('results', [])]
                usage_bucket = ModerationBucket(
                    object=bucket.get('object'),
                    start_time=bucket.get('start_time'),
                    end_time=bucket.get('end_time'),
                    results=results
                )
            elif endpoint == 'completions':
                results = [CompletionResult(**result) for result in bucket.get('results', [])]
                usage_bucket = CompletionBucket(
                    object=bucket.get('object'),
                    start_time=bucket.get('start_time'),
                    end_time=bucket.get('end_time'),
                    results=results
                )
            else:
                continue
            usage_page.data.append(usage_bucket)
        return usage_page

    async def get_moderations_usage(self, **params) -> UsagePage:
        return await self._get_usage('moderations', params)

    async def get_completions_usage(self, **params) -> UsagePage:
        return await self._get_usage('completions', params)

async def main():
    api_key = self.config['api_keys']['OpenAI']['api_key']
    client = OpenAIUsageClient(api_key=api_key)
    start_time = 1730419200  # Example Unix timestamp
    limit = 1
    # Fetch Moderations usage

# To run the example, uncomment the following lines:
if __name__ == '__main__':
     asyncio.run(main())
