''' fine_tuning.py  The purpose of this program is to access OpenAI's v1/fine-tuning endpoint using their python SDK.
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
