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
from lucy.utils.helpers import *
from lucy.utils.load_yaml import load_yaml
from lucy.utils.setup_logging import logger
from openai import AsyncOpenAI

import json
import openai
import traceback

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
