''' benchmark.py  This is not a functioning program.
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
import aiohttp
import asyncio
import time
import datetime
import os
from lucy.utils.create_https_completion import Conversations
from lucy.utils.helpers import load_yaml, PATH_CONFIG_YAML
from lucy.utils.message import Message
from lucy.utils.setup_logging import logger

class Benchmark:
    def __init__(self):
        self.config = load_yaml(PATH_CONFIG_YAML)
        self.conversations = Conversations()
        self.handler = Message(self.config, self.conversations)
        self.initiation_time = datetime.datetime.now()
        self.api_key = self.config['api_keys']['OpenAI]['api_key']
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
        # Implement your quality analysis logic here
        # For example, check if responses meet certain criteria
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

        # Perform any additional analysis
        self.analyze_response_quality()

if __name__ == '__main__':
    benchmark = Benchmark()
    asyncio.run(benchmark.main())
