from openai import AsyncOpenAI
from utils.load_yaml import load_yaml
from utils.setup_logging import logger

import asyncio
import utils.helpers as helpers

config = load_yaml(helpers.PATH_CONFIG_YAML)
api_key = config['api_keys']['api_key_1']['api_key']
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
    asyncio.run(cancel())
