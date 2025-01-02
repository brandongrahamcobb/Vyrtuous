''' create_completion.py  The purpose of this program is to be a simpler implementation of create_https_completion.py from cd ../.
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

async def create_completion(input_array):
    try:
        config = load_yaml(helpers.PATH_CONFIG_YAML)
        api_key = config['api_keys']['api_key_1']['api_key']
        ai_client = AsyncOpenAI(api_key=api_key)
        response = await ai_client.chat.completions.create(
            model='gpt-4o-mini',
            messages=input_array,
            response_format=helpers.OPENAI_CHAT_RESPONSE_FORMAT
        )
        yield response.choices[0].message.content
    except Exception as e:
        yield {'error': traceback.format_exc()}
