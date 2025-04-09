''' create_completion.py  The purpose of this program is to be a simpler implementation of create_https_completion.py.
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
from lucy.utils.load_yaml import load_yaml
from lucy.utils.helpers import *
from lucy.utils.setup_logging import logger
from openai import AsyncOpenAI, OpenAI

import aiohttp
import discord
import io
import openai
import os
import traceback

async def create_image(prompt):
    try:
        config = load_yaml(PATH_CONFIG_YAML)
        api_key = config['api_keys']['OpenAI']['api_key']
        client = AsyncOpenAI(api_key=api_key)
        response = await client.images.generate(
            model='dall-e-3',
            prompt=prompt,
            size='1024x1024',
            quality='standard',
            n=1,
        )
        image_url = response.data[0].url
        return discord.File(await download_image(image_url), filename='generated.png')
    except Exception as e:
        return f'Error: {traceback.format_exc()}'


async def create_image_variation(image_file, prompt):
    try:
        config = load_yaml(PATH_CONFIG_YAML)
        api_key = config['api_keys']['OpenAI']['api_key']
        client = OpenAI(api_key=api_key)
        file_path = os.path.join(DIR_TEMP, 'uploaded_image.png')
        with open(file_path, 'wb') as f:
            f.write(image_file.fp.read())
        with open(file_path, 'rb') as image:
            response = client.images.create_variation(
                model='dall-e-2',
                image=image,
                n=1,
            )
        image_url = response.data[0].url
        return discord.File(await download_image(image_url), filename='generated_variation.png')
    except Exception as e:
        error_msg = f'Error creating image variation: {traceback.format_exc()}'
        logger.error(error_msg)
        return error_msg

async def edit_image(image_file, mask_file, prompt):
    try:
        config = load_yaml(PATH_CONFIG_YAML)
        api_key = config['api_keys']['OpenAI']['api_key']
        image_bytes = image_file.fp.read()
        mask_bytes = await mask_file.read() if mask_file else None
        data = aiohttp.FormData()
        data.add_field('image', image_bytes, filename='image.png', content_type='image/png')
        if mask_bytes:
            data.add_field('mask', mask_bytes, filename='mask.png', content_type='image/png')
        data.add_field('prompt', prompt)
        data.add_field('n', '1')
        data.add_field('size', '1024x1024')
        headers = {
            'Authorization': f'Bearer {api_key}',
        }
        async with aiohttp.ClientSession() as session:
            async with session.post('https://api.openai.com/v1/images/edits', headers=headers, data=data) as response:
                if response.status == 200:
                    logger.info('Received successful response from OpenAI API.')
                    response_json = await response.json()
                    image_url = response_json['data'][0]['url']
                    return discord.File(await download_image(image_url), filename='edited_variation.png')
                else:
                    error_message = await response.text()
                    return None
    except Exception as e:
        logger.error(f'Error during image edit: {e}')
        return None

async def download_image(url):
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as resp:
                if resp.status != 200:
                    return None
                img_data = io.BytesIO(await resp.read())
                img_data.seek(0)
                return img_data
    except Exception as e:
        error_msg = traceback.format_exc()
        logger.error(f'Error downloading image: {error_msg}')  # ✅ Log error
        return None
