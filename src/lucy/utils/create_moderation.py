''' create_moderation.py  The purpose of this program is to be a simpler implementation of create_https_moderation.py from cd ../.
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
from utils.load_yaml import load_yaml
from utils.setup_logging import logger

import openai
import traceback
import utils.helpers as helpers

async def create_moderation(input_text):
    try:
        logger.info('Starting moderation process...')
        config = load_yaml(helpers.PATH_CONFIG_YAML)
        logger.debug('Configuration loaded successfully.')

        api_key = config['api_keys']['api_key_1']
        logger.debug('API key retrieved successfully.')

        ai_client = AsyncOpenAI(api_key=api_key)
        logger.info('AI client initialized.')

        response = await ai_client.moderations.create(
            model='omni-moderation-latest',
            input=input_text,
        )
        logger.info('Moderation API call completed.')

        moderation_response = await response.json()
        logger.debug(f'Moderation response: {moderation_response}')

        yield moderation_response

    except Exception as e:
        error_details = traceback.format_exc()
        logger.error(f'An error occurred during moderation: {error_details}')
        yield {'error': error_details}
