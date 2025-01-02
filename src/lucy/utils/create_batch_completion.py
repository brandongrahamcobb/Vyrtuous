''' create_batch_completionn.py  The purpose of this program is to load a better queue for 50% cheaper API requests during low traffic times from cd ../.
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
from openai import AsyncOpenAI
from utils.setup_logging import logger

import load_yaml
import openai

async def create_batch_completion(batch_input_filename):
    config = load_yaml(helpers.path_config_yaml)
    api_key = config['api_keys']['api_key_1']
    ai_client = AsyncOpenAI(api_key=api_key)
    try:
        with open(batch_input_filename, 'r') as file:
            requests = [json.loads(line) for line in file]
    except FileNotFoundError:
        return False, f'File {batch_input_filename} not found. Please create the file and try again.'
    except json.JSONDecodeError:
        return False, f'Error decoding {batch_input_filename}. Ensure it's in proper JSONL format.'
    try:
        batch_input_file = await ai_client.files.create(file=open(batch_input_filename, 'rb'), purpose='batch')
        batch_input_file_id = batch_input_file.id
    except Exception as e:
        return False, f'Error uploading batch file: {e}'
    try:
        batch = await ai_client.batches.create(
            input_file_id=batch_input_file_id,
            endpoint='/v1/chat/completions',
            completion_window='24h',
            metadata={'description': 'Discord batch processing'}
        )
        batch_id = batch.id  # Use attribute-style access here
    except Exception as e:
        return False, f'Error creating batch: {e}'
    while True:
        batch_status = await ai_client.batches.retrieve(batch_id)
        print(batch_status.errors)
        if batch_status.status == 'completed':
            break
        elif batch_status.status in ['failed', 'cancelled', 'expired']:
            return False, f'Batch failed with status: {batch_status.status}'
    try:
        output_file_id = batch_status.output_file_id  # Use attribute-style access
        output_file = await ai_client.files.content(output_file_id)
        output_data = [json.loads(line) for line in output_file.splitlines()]
        results = '\n'.join(
            f'Request {i + 1}: {res['response']['body']['choices'][0]['message']['content']}'
            for i, res in enumerate(output_data)
        )
        return True, results
    except Exception as e:
        return False, f'Error retrieving batch results: {e}'

