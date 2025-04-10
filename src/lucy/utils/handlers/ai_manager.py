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
from lucy.utils.inc.load_yaml import load_yaml
from lucy.utils.inc.setup_logging import logger
from openai import AsyncOpenAI
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

import aiohttp
import asyncio
import time
import datetime
import os
class Benchmark:
    def __init__(self):
        self.config = load_yaml(PATH_CONFIG_YAML)
        self.conversations = Conversations()
        self.handler = Message(self.config, self.conversations)
        self.initiation_time = datetime.datetime.now()
        self.api_key = self.config['api_keys']['OpenAI']['api_key']
        self.model = self.config['openai_chat_model']
        self.endpoint = OPENAI_ENDPOINT_URLS['chat']
        self.num_requests = OPENAI_CHAT_N
        self.concurrency = ''
        self.prompt = ''
        self.latencies = []
        self.responses = []
        self.success_requests = 0
        self.failed_requests = 0

    async def fetch(self, session, prompt, request_id):
        headers = {
            'Content-Type': 'application/json',
            'Authorization': f'Bearer {self.api_key}',
        }
        payload = {
            'model': self.model,
            'messages': [{'role': 'user', 'content': prompt}],
            'temperature': 0.7,
        }
        start_time = time.time()
        try:
            async with session.post(self.endpoint, json=payload, headers=headers) as response:
                resp = await response.json()
                end_time = time.time()
                latency = end_time - start_time
                if response.status == 200:
                    self.success_requests += 1
                    self.latencies.append(latency)
                    self.responses.append(resp)
                    logger.info(f'Request {request_id} succeeded in {latency:.4f} seconds.')
                else:
                    self.failed_requests += 1
                    self.latencies.append(latency)
                    logger.error(f'Request {request_id} failed with status {response.status}: {resp}')
        except Exception as e:
            end_time = time.time()
            latency = end_time - start_time
            self.failed_requests += 1
            self.latencies.append(latency)
            logger.error(f'Request {request_id} encountered an exception after {latency:.4f} seconds. Error: {e}')

    async def run_benchmark(self):
        connector = aiohttp.TCPConnector(limit=self.concurrency)
        timeout = aiohttp.ClientTimeout(total=60)  # Adjust timeout as needed
        async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
            tasks = [
                asyncio.create_task(self.fetch(session, self.prompt, i))
                for i in range(1, self.num_requests + 1)
            ]
            await asyncio.gather(*tasks)

    def analyze_response_quality(self):
        pass

    async def main(self):
        start_benchmark = time.time()
        logger.info('Starting benchmark...')
        await self.run_benchmark()
        end_benchmark = time.time()
        total_time = end_benchmark - start_benchmark
        avg_latency = sum(self.latencies) / len(self.latencies) if self.latencies else 0
        throughput = self.success_requests / total_time if total_time > 0 else 0
        error_rate = (self.failed_requests / self.num_requests) * 100 if self.num_requests > 0 else 0
        logger.info('----- Benchmark Results -----')
        logger.info(f'Total Requests: {self.num_requests}')
        logger.info(f'Successful Requests: {self.success_requests}')
        logger.info(f'Failed Requests: {self.failed_requests}')
        logger.info(f'Total Time: {total_time:.2f} seconds')
        logger.info(f'Average Latency: {avg_latency:.4f} seconds')
        logger.info(f'Throughput: {throughput:.2f} requests/second')
        logger.info(f'Error Rate: {error_rate:.2f}%')
        logger.info('------------------------------')
        self.analyze_response_quality()

if __name__ == '__main__':
    benchmark = Benchmark()
    asyncio.run(benchmark.main())

"""
API REQUEST PARALLEL PROCESSOR

Using the OpenAI API to process lots of text quickly takes some care.
If you trickle in a million API requests one by one, they'll take days to complete.
If you flood a million API requests in parallel, they'll exceed the rate limits and fail with errors.
To maximize throughput, parallel requests need to be throttled to stay under rate limits.

This script parallelizes requests to the OpenAI API while throttling to stay under rate limits.

Features:
- Streams requests from file, to avoid running out of memory for giant jobs
- Makes requests concurrently, to maximize throughput
- Throttles request and token usage, to stay under rate limits
- Retries failed requests up to {max_attempts} times, to avoid missing data
- Logs errors, to diagnose problems with requests

Example command to call script:
```
python examples/api_request_parallel_processor.py \
  --requests_filepath examples/data/example_requests_to_parallel_process.jsonl \
  --save_filepath examples/data/example_requests_to_parallel_process_results.jsonl \
  --request_url https://api.openai.com/v1/embeddings \
  --max_requests_per_minute 1500 \
  --max_tokens_per_minute 6250000 \
  --token_encoding_name cl100k_base \
  --max_attempts 5 \
  --logging_level 20
```

Inputs:
- requests_filepath : str
    - path to the file containing the requests to be processed
    - file should be a jsonl file, where each line is a json object with API parameters and an optional metadata field
    - e.g., {"model": "text-embedding-3-small", "input": "embed me", "metadata": {"row_id": 1}}
    - as with all jsonl files, take care that newlines in the content are properly escaped (json.dumps does this automatically)
    - an example file is provided at examples/data/example_requests_to_parallel_process.jsonl
    - the code to generate the example file is appended to the bottom of this script
- save_filepath : str, optional
    - path to the file where the results will be saved
    - file will be a jsonl file, where each line is an array with the original request plus the API response
    - e.g., [{"model": "text-embedding-3-small", "input": "embed me"}, {...}]
    - if omitted, results will be saved to {requests_filename}_results.jsonl
- request_url : str, optional
    - URL of the API endpoint to call
    - if omitted, will default to "https://api.openai.com/v1/embeddings"
- api_key : str, optional
    - API key to use
    - if omitted, the script will attempt to read it from an environment variable {os.getenv("OPENAI_API_KEY")}
- max_requests_per_minute : float, optional
    - target number of requests to make per minute (will make less if limited by tokens)
    - leave headroom by setting this to 50% or 75% of your limit
    - if requests are limiting you, try batching multiple embeddings or completions into one request
    - if omitted, will default to 1,500
- max_tokens_per_minute : float, optional
    - target number of tokens to use per minute (will use less if limited by requests)
    - leave headroom by setting this to 50% or 75% of your limit
    - if omitted, will default to 125,000
- token_encoding_name : str, optional
    - name of the token encoding used, as defined in the `tiktoken` package
    - if omitted, will default to "cl100k_base" (used by `text-embedding-3-small`)
- max_attempts : int, optional
    - number of times to retry a failed request before giving up
    - if omitted, will default to 5
- logging_level : int, optional
    - level of logging to use; higher numbers will log fewer messages
    - 40 = ERROR; will log only when requests fail after all retries
    - 30 = WARNING; will log when requests his rate limits or other errors
    - 20 = INFO; will log when requests start and the status at finish
    - 10 = DEBUG; will log various things as the loop runs to see when they occur
    - if omitted, will default to 20 (INFO).

The script is structured as follows:
    - Imports
    - Define main()
        - Initialize things
        - In main loop:
            - Get next request if one is not already waiting for capacity
            - Update available token & request capacity
            - If enough capacity available, call API
            - The loop pauses if a rate limit error is hit
            - The loop breaks when no tasks remain
    - Define dataclasses
        - StatusTracker (stores script metadata counters; only one instance is created)
        - APIRequest (stores API inputs, outputs, metadata; one method to call API)
    - Define functions
        - api_endpoint_from_url (extracts API endpoint from request URL)
        - append_to_jsonl (writes to results file)
        - num_tokens_consumed_from_request (bigger function to infer token usage from request)
        - task_id_generator_function (yields 0, 1, 2, ...)
    - Run main()
"""

# imports
import aiohttp  # for making API calls concurrently
import argparse  # for running script from command line
import asyncio  # for running API calls concurrently
import json  # for saving results to a jsonl file
import logging  # for logging rate limit warnings and other messages
import os  # for reading API key
import re  # for matching endpoint from request URL
import tiktoken  # for counting tokens
import time  # for sleeping after rate limit is hit
from dataclasses import (
    dataclass,
    field,
)  # for storing API inputs, outputs, and metadata


async def process_api_requests_from_file(
    requests_filepath: str,
    save_filepath: str,
    request_url: str,
    api_key: str,
    max_requests_per_minute: float,
    max_tokens_per_minute: float,
    token_encoding_name: str,
    max_attempts: int,
    logging_level: int,
    status_tracker=None
):
    if status_tracker is None:
        status_tracker = StatusTracker()
    """Processes API requests in parallel, throttling to stay under rate limits."""
    # constants
    seconds_to_pause_after_rate_limit_error = 15
    seconds_to_sleep_each_loop = (
        0.001  # 1 ms limits max throughput to 1,000 requests per second
    )

    # initialize logging
    logging.basicConfig(level=logging_level)
    logging.debug(f"Logging initialized at level {logging_level}")

    # infer API endpoint and construct request header
    api_endpoint = api_endpoint_from_url(request_url)
    request_header = {"Authorization": f"Bearer {api_key}"}
    # use api-key header for Azure deployments
    if "/deployments" in request_url:
        request_header = {"api-key": f"{api_key}"}

    # initialize trackers
    queue_of_requests_to_retry = asyncio.Queue()
    task_id_generator = (
        task_id_generator_function()
    )  # generates integer IDs of 0, 1, 2, ...
    status_tracker = (
        StatusTracker()
    )  # single instance to track a collection of variables
    next_request = None  # variable to hold the next request to call

    # initialize available capacity counts
    available_request_capacity = max_requests_per_minute
    available_token_capacity = max_tokens_per_minute
    last_update_time = time.time()

    # initialize flags
    file_not_finished = True  # after file is empty, we'll skip reading it
    logging.debug(f"Initialization complete.")

    # initialize file reading
    with open(requests_filepath) as file:
        # `requests` will provide requests one at a time
        requests = file.__iter__()
        logging.debug(f"File opened. Entering main loop")
        async with aiohttp.ClientSession() as session:  # Initialize ClientSession here
            while True:
                # get next request (if one is not already waiting for capacity)
                if next_request is None:
                    if not queue_of_requests_to_retry.empty():
                        next_request = queue_of_requests_to_retry.get_nowait()
                        logging.debug(
                            f"Retrying request {next_request.task_id}: {next_request}"
                        )
                    elif file_not_finished:
                        try:
                            # get new request
                            request_json = json.loads(next(requests))
                            next_request = APIRequest(
                                task_id=next(task_id_generator),
                                request_json=request_json,
                                token_consumption=num_tokens_consumed_from_request(
                                    request_json, api_endpoint, token_encoding_name
                                ),
                                attempts_left=max_attempts,
                                metadata=request_json.pop("metadata", None),
                            )
                            status_tracker.num_tasks_started += 1
                            status_tracker.num_tasks_in_progress += 1
                            logging.debug(
                                f"Reading request {next_request.task_id}: {next_request}"
                            )
                        except StopIteration:
                            # if file runs out, set flag to stop reading it
                            logging.debug("Read file exhausted")
                            file_not_finished = False

                # update available capacity
                current_time = time.time()
                seconds_since_update = current_time - last_update_time
                available_request_capacity = min(
                    available_request_capacity
                    + max_requests_per_minute * seconds_since_update / 60.0,
                    max_requests_per_minute,
                )
                available_token_capacity = min(
                    available_token_capacity
                    + max_tokens_per_minute * seconds_since_update / 60.0,
                    max_tokens_per_minute,
                )
                last_update_time = current_time

                # if enough capacity available, call API
                if next_request:
                    next_request_tokens = next_request.token_consumption
                    if (
                        available_request_capacity >= 1
                        and available_token_capacity >= next_request_tokens
                    ):
                        # update counters
                        available_request_capacity -= 1
                        available_token_capacity -= next_request_tokens
                        next_request.attempts_left -= 1

                        # call API
                        asyncio.create_task(
                            next_request.call_api(
                                session=session,
                                request_url=request_url,
                                request_header=request_header,
                                retry_queue=queue_of_requests_to_retry,
                                save_filepath=save_filepath,
                                status_tracker=status_tracker,
                            )
                        )
                        next_request = None  # reset next_request to empty

                # if all tasks are finished, break
                if status_tracker.num_tasks_in_progress == 0:
                    break

                # main loop sleeps briefly so concurrent tasks can run
                await asyncio.sleep(seconds_to_sleep_each_loop)

                # if a rate limit error was hit recently, pause to cool down
                seconds_since_rate_limit_error = (
                    time.time() - status_tracker.time_of_last_rate_limit_error
                )
                if (
                    seconds_since_rate_limit_error
                    < seconds_to_pause_after_rate_limit_error
                ):
                    remaining_seconds_to_pause = (
                        seconds_to_pause_after_rate_limit_error
                        - seconds_since_rate_limit_error
                    )
                    await asyncio.sleep(remaining_seconds_to_pause)
                    # ^e.g., if pause is 15 seconds and final limit was hit 5 seconds ago
                    logging.warn(
                        f"Pausing to cool down until {time.ctime(status_tracker.time_of_last_rate_limit_error + seconds_to_pause_after_rate_limit_error)}"
                    )

        # after finishing, log final status
        logging.info(
            f"""Parallel processing complete. Results saved to {save_filepath}"""
        )
        if status_tracker.num_tasks_failed > 0:
            logging.warning(
                f"{status_tracker.num_tasks_failed} / {status_tracker.num_tasks_started} requests failed. Errors logged to {save_filepath}."
            )
        if status_tracker.num_rate_limit_errors > 0:
            logging.warning(
                f"{status_tracker.num_rate_limit_errors} rate limit errors received. Consider running at a lower rate."
            )


# dataclasses


@dataclass
class StatusTracker:
    """Stores metadata about the script's progress. Only one instance is created."""

    num_tasks_started: int = 0
    num_tasks_in_progress: int = 0  # script ends when this reaches 0
    num_tasks_succeeded: int = 0
    num_tasks_failed: int = 0
    num_rate_limit_errors: int = 0
    num_api_errors: int = 0  # excluding rate limit errors, counted above
    num_other_errors: int = 0
    time_of_last_rate_limit_error: int = 0  # used to cool off after hitting rate limits


@dataclass
class APIRequest:
    """Stores an API request's inputs, outputs, and other metadata. Contains a method to make an API call."""

    task_id: int
    request_json: dict
    token_consumption: int
    attempts_left: int
    metadata: dict
    result: list = field(default_factory=list)

    async def call_api(
        self,
        session: aiohttp.ClientSession,
        request_url: str,
        request_header: dict,
        retry_queue: asyncio.Queue,
        save_filepath: str,
        status_tracker: StatusTracker,
    ):
        """Calls the OpenAI API and saves results."""
        logging.info(f"Starting request #{self.task_id}")
        error = None
        try:
            async with session.post(
                url=request_url, headers=request_header, json=self.request_json
            ) as response:
                response = await response.json()
            if "error" in response:
                logging.warning(
                    f"Request {self.task_id} failed with error {response['error']}"
                )
                status_tracker.num_api_errors += 1
                error = response
                if "rate limit" in response["error"].get("message", "").lower():
                    status_tracker.time_of_last_rate_limit_error = time.time()
                    status_tracker.num_rate_limit_errors += 1
                    status_tracker.num_api_errors -= (
                        1  # rate limit errors are counted separately
                    )

        except (
            Exception
        ) as e:  # catching naked exceptions is bad practice, but in this case we'll log & save them
            logging.warning(f"Request {self.task_id} failed with Exception {e}")
            status_tracker.num_other_errors += 1
            error = e
        if error:
            self.result.append(error)
            if self.attempts_left:
                retry_queue.put_nowait(self)
            else:
                logging.error(
                    f"Request {self.request_json} failed after all attempts. Saving errors: {self.result}"
                )
                data = (
                    [self.request_json, [str(e) for e in self.result], self.metadata]
                    if self.metadata
                    else [self.request_json, [str(e) for e in self.result]]
                )
                append_to_jsonl(data, save_filepath)
                status_tracker.num_tasks_in_progress -= 1
                status_tracker.num_tasks_failed += 1
        else:
            data = (
                [self.request_json, response, self.metadata]
                if self.metadata
                else [self.request_json, response]
            )
            append_to_jsonl(data, save_filepath)
            status_tracker.num_tasks_in_progress -= 1
            status_tracker.num_tasks_succeeded += 1
            logging.debug(f"Request {self.task_id} saved to {save_filepath}")


# functions


def api_endpoint_from_url(request_url):
    """Extract the API endpoint from the request URL."""
    match = re.search("^https://[^/]+/v\\d+/(.+)$", request_url)
    if match is None:
        # for Azure OpenAI deployment urls
        match = re.search(
            r"^https://[^/]+/openai/deployments/[^/]+/(.+?)(\?|$)", request_url
        )
    return match[1]


def append_to_jsonl(data, filename: str) -> None:
    """Append a json payload to the end of a jsonl file."""
    json_string = json.dumps(data)
    with open(filename, "a") as f:
        f.write(json_string + "\n")


def num_tokens_consumed_from_request(
    request_json: dict,
    api_endpoint: str,
    token_encoding_name: str,
):
    """Count the number of tokens in the request. Only supports completion and embedding requests."""
    encoding = tiktoken.get_encoding(token_encoding_name)
    # if completions request, tokens = prompt + n * max_tokens
    if api_endpoint.endswith("completions"):
        max_tokens = request_json.get("max_tokens", 15)
        n = request_json.get("n", 1)
        completion_tokens = n * max_tokens

        # chat completions
        if api_endpoint.startswith("chat/"):
            num_tokens = 0
            for message in request_json["messages"]:
                num_tokens += 4  # every message follows <im_start>{role/name}\n{content}<im_end>\n
                for key, value in message.items():
                    num_tokens += len(encoding.encode(value))
                    if key == "name":  # if there's a name, the role is omitted
                        num_tokens -= 1  # role is always required and always 1 token
            num_tokens += 2  # every reply is primed with <im_start>assistant
            return num_tokens + completion_tokens
        # normal completions
        else:
            prompt = request_json["prompt"]
            if isinstance(prompt, str):  # single prompt
                prompt_tokens = len(encoding.encode(prompt))
                num_tokens = prompt_tokens + completion_tokens
                return num_tokens
            elif isinstance(prompt, list):  # multiple prompts
                prompt_tokens = sum([len(encoding.encode(p)) for p in prompt])
                num_tokens = prompt_tokens + completion_tokens * len(prompt)
                return num_tokens
            else:
                raise TypeError(
                    'Expecting either string or list of strings for "prompt" field in completion request'
                )
    # if embeddings request, tokens = input tokens
    elif api_endpoint == "embeddings":
        input = request_json["input"]
        if isinstance(input, str):  # single input
            num_tokens = len(encoding.encode(input))
            return num_tokens
        elif isinstance(input, list):  # multiple inputs
            num_tokens = sum([len(encoding.encode(i)) for i in input])
            return num_tokens
        else:
            raise TypeError(
                'Expecting either string or list of strings for "inputs" field in embedding request'
            )
    # more logic needed to support other API calls (e.g., edits, inserts, DALL-E)
    else:
        raise NotImplementedError(
            f'API endpoint "{api_endpoint}" not implemented in this script'
        )


def task_id_generator_function():
    """Generate integers 0, 1, 2, and so on."""
    task_id = 0
    while True:
        yield task_id
        task_id += 1


# run script


if __name__ == "__main__":
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--requests_filepath")
    parser.add_argument("--save_filepath", default=None)
    parser.add_argument("--request_url", default="https://api.openai.com/v1/embeddings")
    parser.add_argument("--api_key", default=os.getenv("OPENAI_API_KEY"))
    parser.add_argument("--max_requests_per_minute", type=int, default=3_000 * 0.5)
    parser.add_argument("--max_tokens_per_minute", type=int, default=250_000 * 0.5)
    parser.add_argument("--token_encoding_name", default="cl100k_base")
    parser.add_argument("--max_attempts", type=int, default=5)
    parser.add_argument("--logging_level", default=logging.INFO)
    args = parser.parse_args()

    if args.save_filepath is None:
        args.save_filepath = args.requests_filepath.replace(".jsonl", "_results.jsonl")

    # run script
    asyncio.run(
        process_api_requests_from_file(
            requests_filepath=args.requests_filepath,
            save_filepath=args.save_filepath,
            request_url=args.request_url,
            api_key=args.api_key,
            max_requests_per_minute=float(args.max_requests_per_minute),
            max_tokens_per_minute=float(args.max_tokens_per_minute),
            token_encoding_name=args.token_encoding_name,
            max_attempts=int(args.max_attempts),
            logging_level=int(args.logging_level),
        )
    )


"""
APPENDIX

The example requests file at openai-cookbook/examples/data/example_requests_to_parallel_process.jsonl contains 10,000 requests to text-embedding-3-small.

It was generated with the following code:

```python
import json

filename = "data/example_requests_to_parallel_process.jsonl"
n_requests = 10_000
jobs = [{"model": "text-embedding-3-small", "input": str(x) + "\n"} for x in range(n_requests)]
with open(filename, "w") as f:
    for job in jobs:
        json_string = json.dumps(job)
        f.write(json_string + "\n")
```

As with all jsonl files, take care that newlines in the content are properly escaped (json.dumps does this automatically).
"""
from openai import AsyncOpenAI

import asyncio

config = load_yaml(PATH_CONFIG_YAML)
api_key = config['api_keys']['OpenAI']['api_key']
ai_client = AsyncOpenAI(api_key=api_key)

async def cancel():
    await ai_client.fine_tuning.jobs.cancel('ftjob-VBRw83PIls4zA25bypQBcCHH')

async def main():
    await ai_client.fine_tuning.jobs.create(
        training_file='file-LvuzigtnKkifPazQptC7Mz',
        model='gpt-4o-mini-2024-07-18',
        suffix='vyrtuous'
    )

if __name__ == '__main__':
    asyncio.run(main())

from collections import defaultdict

import json
import numpy as np
import tiktoken

encoding = tiktoken.encoding_for_model('gpt-4o-mini-07-18') #get_encoding('cl100k_base')

# not exact!
# simplified from https://github.com/openai/openai-cookbook/blob/main/examples/How_to_count_tokens_with_tiktoken.ipynb
def num_tokens_from_messages(messages, tokens_per_message=3, tokens_per_name=1):
    num_tokens = 0
    for message in messages:
        num_tokens += tokens_per_message
        for key, value in message.items():
            if isinstance(value, str):  # Check if the value is a string
                num_tokens += len(encoding.encode(value))
            else:
                print(f'Warning: Non-string value {value} of type {type(value)} encountered. Skipping.')
            if key == 'name':
                num_tokens += tokens_per_name
    num_tokens += 3  # For end-of-sequence token
    return num_tokens

def num_assistant_tokens_from_messages(messages):
    num_tokens = 0
    for message in messages:
        if message['role'] == 'assistant':
            content = message.get('content', '')
            if isinstance(content, str):
                num_tokens += len(encoding.encode(content))
            else:
                print(f'Warning: Non-string assistant content {content} of type {type(content)} encountered. Skipping.')
    return num_tokens

def print_distribution(values, name):
    print(f'\n#### Distribution of {name}:')
    print(f'min / max: {min(values)}, {max(values)}')
    print(f'mean / median: {np.mean(values)}, {np.median(values)}')
    print(f'p5 / p95: {np.quantile(values, 0.1)}, {np.quantile(values, 0.9)}')

def format_error_check():
    format_errors = defaultdict(int)
    for ex in dataset:
        if not isinstance(ex, dict):
            format_errors['data_type'] += 1
            continue
        messages = ex.get('messages', None)
        if not messages:
            format_errors['missing_messages_list'] += 1
            continue
        for message in messages:
            if 'role' not in message or 'content' not in message:
                format_errors['message_missing_key'] += 1
            if any(k not in ('role', 'content', 'name', 'function_call', 'weight') for k in message):
                format_errors['message_unrecognized_key'] += 1
            if message.get('role', None) not in ('system', 'user', 'assistant', 'function'):
                format_errors['unrecognized_role'] += 1
            content = message.get('content', None)
            function_call = message.get('function_call', None)
            if (not content and not function_call) or not isinstance(content, str):
                format_errors['missing_content'] += 1
        if not any(message.get('role', None) == 'assistant' for message in messages):
            format_errors['example_missing_assistant_message'] += 1
    if format_errors:
        print('Found errors:')
        for k, v in format_errors.items():
            print(f'{k}: {v}')
    else:
        print('No errors found')

if __name__ == '__main__':
    data_path = 'training_temp.jsonl'
    with open(data_path, 'r', encoding='utf-8') as f:
       dataset = [json.loads(line) for line in f]
#    print('Num examples:', len(dataset))
 #   print('First example:')
  #  for message in dataset[0]['messages']:
   #     print(message)
    # Warnings and tokens counts
    n_missing_system = 0
    n_missing_user = 0
    n_messages = []
    convo_lens = []
    assistant_message_lens = []
    
    for ex in dataset:
        messages = ex['messages']
        if not any(message['role'] == 'system' for message in messages):
            n_missing_system += 1
        if not any(message['role'] == 'user' for message in messages):
            n_missing_user += 1
        n_messages.append(len(messages))
        convo_lens.append(num_tokens_from_messages(messages))
        assistant_message_lens.append(num_assistant_tokens_from_messages(messages))
        
    print('Num examples missing system message:', n_missing_system)
    print('Num examples missing user message:', n_missing_user)
    print_distribution(n_messages, 'num_messages_per_example')
    print_distribution(convo_lens, 'num_total_tokens_per_example')
    print_distribution(assistant_message_lens, 'num_assistant_tokens_per_example')
    n_too_long = sum(l > 16385 for l in convo_lens)
    print(f'\n{n_too_long} examples may be over the 16,385 token limit, they will be truncated during fine-tuning') 
    # Pricing and default n_epochs estimate
    MAX_TOKENS_PER_EXAMPLE = 16385
    
    TARGET_EPOCHS = 3
    MIN_TARGET_EXAMPLES = 100
    MAX_TARGET_EXAMPLES = 25000
    MIN_DEFAULT_EPOCHS = 1
    MAX_DEFAULT_EPOCHS = 25
    
    n_epochs = TARGET_EPOCHS
    n_train_examples = len(dataset)
    if n_train_examples * TARGET_EPOCHS < MIN_TARGET_EXAMPLES:
        n_epochs = min(MAX_DEFAULT_EPOCHS, MIN_TARGET_EXAMPLES // n_train_examples)
    elif n_train_examples * TARGET_EPOCHS > MAX_TARGET_EXAMPLES:
        n_epochs = max(MIN_DEFAULT_EPOCHS, MAX_TARGET_EXAMPLES // n_train_examples)
    
    n_billing_tokens_in_dataset = sum(min(MAX_TOKENS_PER_EXAMPLE, length) for length in convo_lens)
    print(f'Dataset has ~{n_billing_tokens_in_dataset} tokens that will be charged for during training')
    print(f'By default, you\'ll train for {n_epochs} epochs on this dataset')
    print(f'By default, you\'ll be charged for ~{n_epochs * n_billing_tokens_in_dataset} tokens')
