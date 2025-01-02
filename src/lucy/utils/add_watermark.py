''' add_watermark.py  The purpose of this program is to provide a watermark to PIL Image objects from cd ../.
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
from PIL import Image, ImageDraw, ImageFont
import math
from utils.setup_logging import logger

def add_watermark(image: BytesIO, watermark_text: str = 'Discord') -> BytesIO:
    logger.info('Starting the watermarking process.')
    
    try:
        RGB_image = Image.open(image)
        RGBA_image = RGB_image.convert('RGBA')
        logger.info('Image loaded and converted to RGBA.')

        draw = ImageDraw.Draw(RGBA_image)
        width, height = RGBA_image.size
        logger.info(f'Image dimensions: width={width}, height={height}.')

        diagonal = math.sqrt(width**2 + height**2)
        font_size = int(diagonal / 15)
        logger.info(f'Calculated initial font size: {font_size}.')

        try:
            font = ImageFont.truetype('Roboto-Regular.ttf', font_size)
            logger.info('Loaded Roboto-Regular.ttf font.')
        except IOError:
            logger.warning('Roboto-Regular.ttf not found. Falling back to default font.')
            font = ImageFont.load_default()

        while True:
            bbox = draw.textbbox((0, 0), watermark_text, font=font)
            text_width = bbox[2] - bbox[0]
            if text_width <= 512 or font_size <= 1:
                break
            font_size -= 1
            font = ImageFont.truetype('Roboto-Regular.ttf', font_size)
            logger.debug(f'Font size reduced to: {font_size}.')

        text_height = bbox[3] - bbox[1]
        text_x = (width - text_width) / 2
        text_y = height - (2 * text_height)
        logger.info(f'Text position calculated: x={text_x}, y={text_y}.')

        watermark_image = Image.new('RGBA', RGBA_image.size, (0, 0, 0, 0))
        watermark_draw = ImageDraw.Draw(watermark_image)
        watermark_draw.text((text_x, text_y), watermark_text, font=font, fill=(255, 255, 255, 64))
        logger.info('Watermark text added to the watermark image.')

        mask = watermark_image.split()[3]
        RGBA_image.paste(watermark_image, (0, 0), mask)
        logger.info('Watermark image pasted onto the original image.')

        output = BytesIO()
        RGBA_image.save(output, format='PNG')
        output.seek(0)
        logger.info('Watermarked image saved to output stream.')

        return output

    except Exception as e:
        logger.error('An error occurred during the watermarking process.', exc_info=True)
        raise

