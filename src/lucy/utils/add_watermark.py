''' add_watermark.py  The purpose of this program is to provide a watermark to PIL Image objects.
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
from io import BytesIO
from lucy.utils.helpers import *
from lucy.utils.setup_logging import logger
from PIL import Image, ImageDraw, ImageFont
from os.path import join

import math

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

def normalize_text(text: str) -> str:
    letters_only = ''.join(filter(str.isalpha, text))
    if letters_only.isupper():
        return text
    else:
        return text.lower().capitalize()
