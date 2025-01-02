''' create_https_moderation.py Still, OpenAI's python SDK is way better than this. However,
                               I'm proud of this program. It's my current work of art.
                               Functioning from cd ../
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
from datetime import datetime
from openai import AsyncOpenAI
from utils.load_yaml import load_yaml
from utils.nlp_utils import NLPUtils
from utils.setup_logging import logger

import aiohttp
import json
import openai
import traceback
import utils.helpers as helpers

async def create_https_moderation(custom_id, input_array, model):
    try:
        logger.info('Loading configuration file.')
        # Load the configuration file
        config = load_yaml(helpers.PATH_CONFIG_YAML)
        api_key = config['api_keys']['api_key_1']['api_key']
        logger.info('API key loaded successfully.')

        # Set up the OpenAI client
        ai_client = AsyncOpenAI(api_key=api_key)
        headers = {'Authorization': f'Bearer {api_key}'}
        logger.info('Headers set for the request.')

        # Prepare the request data
        request_data = {
            'input': input_array,
            'model': model,
        }
        logger.info('Request data prepared.')

        # Send the request to the moderation endpoint
        async with aiohttp.ClientSession() as session:
            try:
                logger.info('Sending request to OpenAI moderation endpoint.')
                async with session.post(url=helpers.OPENAI_ENDPOINT_URLS['moderations'], headers=headers, json=request_data) as moderation_object:
                    logger.info(f'Received response with status: {moderation_object.status}')

                    if moderation_object.status == 200:
                        # If the response is successful, parse and return the response data
                        response_data = await moderation_object.json()
                        logger.info('Request successful. Returning response data.')
                        yield response_data
                    else:
                        # If the response is not successful, log the error message
                        error_message = await moderation_object.text()
                        logger.error(f'Error in response: {error_message}')
                        yield {'error': error_message}

            except Exception as e:
                # If there is an exception during the request, log the exception
                logger.error('An error occurred while making the HTTP request.')
                logger.error(traceback.format_exc())
                yield traceback.format_exc()

    except Exception as e:
        # Log any error that happens in the main function
        logger.error('An error occurred in create_https_moderation.')
        logger.error(traceback.format_exc())
        yield traceback.format_exc()


