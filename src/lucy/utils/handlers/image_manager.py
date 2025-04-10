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
from gradio_client import Client
from itertools import zip_longest
from io import BytesIO
from lucy.utils.inc.load_yaml import load_yaml
from lucy.utils.helpers import *
from lucy.utils.inc.setup_logging import logger
from math import ceil, sqrt
from openai import AsyncOpenAI, OpenAI
from os.path import join
from PIL import Image, ImageDraw, ImageEnhance, ImageFont
from random import randint

import aiohttp
import colorsys
import discord
import io
import math
import openai
import os
import traceback

def add_watermark(image: BytesIO, watermark_text: str = 'Vyrtuous', bottom: bool = True) -> BytesIO:
    if not bottom:
        normalized_text = normalize_text(watermark_text)
    else:
        normalized_text = watermark_text
    try:
        base_image = Image.open(image).convert('RGBA')
        width, height = base_image.size
        diagonal = math.sqrt(width**2 + height**2)
        font_size = int(diagonal / 15)
        try:
            font = ImageFont.truetype(PATH_FONT, font_size)
        except IOError:
            logger.warning('Roboto-Regular.ttf not found. Falling back to default font.')
            font = ImageFont.load_default()
        max_text_width = int(width * 0.8)
        min_font_size = 30
        while True:
            draw = ImageDraw.Draw(Image.new('RGBA', (1, 1)))  # Dummy image to get text size
            bbox = draw.textbbox((0, 0), normalized_text, font=font)
            text_width = bbox[2] - bbox[0]
            if text_width <= max_text_width:
                 break
            font_size -= 1
            font = ImageFont.truetype(PATH_FONT, font_size)
        text_height = bbox[3] - bbox[1]
        text_x = (width - text_width) / 2
        if bottom:
            text_y = height - (2 * text_height)
        else:
            text_y = text_height  # Position near the top
        watermark_layer = Image.new('RGBA', base_image.size, (255, 255, 255, 0))
        draw = ImageDraw.Draw(watermark_layer)
        draw.text((text_x, text_y), normalized_text, font=font, fill=(255, 255, 255, 128))  # White text with transparency
        watermarked_image = Image.alpha_composite(base_image, watermark_layer)
        output = BytesIO()
        watermarked_image.save(output, format='PNG')
        output.seek(0)
        return output
    except Exception as e:
        logger.error('An error occurred during the watermarking process.', exc_info=True)
        raise

def adjust_hue_and_saturation(image, hue_shift, saturation_shift) -> BytesIO:
    try:
        image = image.convert('RGB')
        pixels = list(image.getdata())
        adjusted_pixels = []
        for r, g, b in pixels:
            h, s, v = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
            h = (h + hue_shift / 360.0) % 1.0
            s = min(max(s + saturation_shift / 100.0, 0), 1)
            r, g, b = colorsys.hsv_to_rgb(h, s, v)
            adjusted_pixels.append((int(r * 255), int(g * 255), int(b * 255)))
        new_image = Image.new('RGB', image.size)
        new_image.putdata(adjusted_pixels)
        output = BytesIO()
        new_image.save(output, format='PNG')
        output.seek(0)
        return output
    except Exception as e:
        logger.error(f'An error occurred during hue and saturation adjustment: {e}')
        raise

def combine_gallery(images: list, names: list, title: str, quantity: int = 1, linearity: bool = False) -> BytesIO:
    if len(images) <= 2:
        linearity = True
    processed_images = []
    repeated_images = list(images * quantity)
    for img_bytes in repeated_images:
        img_bytes.seek(0)
        img = Image.open(img_bytes).convert('RGBA')
        processed_images.append(img)
    if linearity:
        combined_width = sum(img.size[0] for img in processed_images)
        combined_height = max(img.size[1] for img in processed_images)
        footer_ratio = 0.15
        footer_height = int(combined_height * footer_ratio)
        new_height = combined_height + footer_height
        combined_img = Image.new('RGB', (combined_width, new_height), (0, 0, 0))
        x_offset = 0
        for img in processed_images:
            y_offset = 0
            img = img.resize((img.size[0], combined_height), Image.LANCZOS)
            combined_img.paste(img, (x_offset, y_offset), img)
            x_offset += img.size[0]
    else:
        num_images = len(processed_images)
        grid_size = ceil(sqrt(num_images))
        image_size = max(img.size[0] for img in processed_images)
        canvas_size = grid_size * image_size
        footer_ratio = 0.15
        footer_height = int(canvas_size * footer_ratio)
        new_height = canvas_size + footer_height
        combined_img = Image.new('RGB', (canvas_size, new_height), (0, 0, 0))
        for index, img in enumerate(processed_images):
            row = index // grid_size
            col = index % grid_size
            x_offset = col * image_size
            y_offset = row * image_size
            img = img.resize((image_size, image_size), Image.LANCZOS)
            combined_img.paste(img, (x_offset, y_offset), img)
    combined_img_buffer = BytesIO()
    combined_img.save(combined_img_buffer, format='PNG')
    combined_img_buffer.seek(0)
    final_image_buffer = add_watermark(combined_img_buffer, title, True)
    final_image_buffer.seek(0)
    return final_image_buffer

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

def normalize_text(text: str) -> str:
    letters_only = ''.join(filter(str.isalpha, text))
    if letters_only.isupper():
        return text
    else:
        return text.lower().capitalize()

def stable_cascade(prompt):
    try:
        client = Client('multimodalart/stable-cascade')
        result = client.predict(
            prompt=prompt,
            negative_prompt='',
            seed=randint(0, 2147483647),
            width=1024,
            height=1024,
            prior_num_inference_steps=20,
            prior_guidance_scale=4,
            decoder_num_inference_steps=10,
            decoder_guidance_scale=0,
            num_images_per_prompt=1,
            api_name='/run',
        )
        return discord.File(result, 'image.webp')
    except ConnectionError as conn_err:
        print(f'Connection error: {conn_err}')
        return 'Failed to connect to the server. Please try again later.'
    except Exception as e:
        print(f'An unexpected error occurred: {e}')
        return f'An error occurred: {e}'
