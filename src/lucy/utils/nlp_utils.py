''' nlp_utils.py  The purpose of this program is to provide generic Natural Language Processing functionality from cd ../.
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
from .setup_logging import logger

import datetime
import json
import os

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
            logger.info(f'Preparing to append entry to JSONL file: {file_path}')
            logger.debug(f'Entry content: {entry}')
            with open(file_path, 'a') as file:
                json.dump(entry, file)
                file.write('\n')
                logger.info(f'Successfully appended entry to JSONL file: {file_path}')
        except Exception as e:
            logger.error(f'Error occurred while appending entry: {e}')
            return {'error': str(e)}
