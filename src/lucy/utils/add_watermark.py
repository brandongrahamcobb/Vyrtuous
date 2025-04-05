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

#def add_watermark(image: BytesIO, watermark_text: str = 'Discord') -> BytesIO:
#
#    logger.info('Starting the watermarking process.')
#
#    try:
#        RGB_image = Image.open(image)
#        RGBA_image = RGB_image.convert('RGBA')
#        logger.info('Image loaded and converted to RGBA.')
#        draw = ImageDraw.Draw(RGBA_image)
#        width, height = RGBA_image.size
#        logger.info(f'Image dimensions: width={width}, height={height}.')
#        diagonal = math.sqrt(width**2 + height**2)
#        font_size = int(diagonal / 15)
#        logger.info(f'Calculated initial font size: {font_size}.')
#        try:
#            font = ImageFont.truetype(PATH_FONT, font_size)
#            logger.info('Loaded Roboto-Regular.ttf font.')
#        except IOError:
#            logger.warning('Roboto-Regular.ttf not found. Falling back to default font.')
#            font = ImageFont.load_default()
#        min_font_size = 30
#        max_text_width = 512
#        
#        while True:
#            bbox = draw.textbbox((0, 0), watermark_text, font=font)
#            text_width = bbox[2] - bbox[0]
#            
#            if text_width <= max_text_width or font_size <= min_font_size:
#                break
#        
#            font_size -= 1
#            if font_size < min_font_size:
#                font_size = min_font_size
#                break  # stop shrinking if it reaches the minimum font size
#        
#            font = ImageFont.truetype(PATH_FONT, font_size)
#            logger.debug(f'Font size reduced to: {font_size}.')
##        while True:
##            bbox = draw.textbbox((0, 0), watermark_text, font=font)
##            text_width = bbox[2] - bbox[0]
##            if text_width <= 512 or font_size <= 1:
##                break
##            font_size -= 1
##            font = ImageFont.truetype('Roboto-Regular.ttf', font_size)
##            logger.debug(f'Font size reduced to: {font_size}.')
#        text_height = bbox[3] - bbox[1]
#        text_x = (width - text_width) / 2
#        text_y = height - (2 * text_height)
#        logger.info(f'Text position calculated: x={text_x}, y={text_y}.')
#        watermark_image = Image.new('RGBA', RGBA_image.size, (0, 0, 0, 0))
#        watermark_draw = ImageDraw.Draw(watermark_image)
#        watermark_draw.text((text_x, text_y), watermark_text, font=font, fill=(255, 255, 255, 64))
#        logger.info('Watermark text added to the watermark image.')
#        mask = watermark_image.split()[3]
#        RGBA_image.paste(watermark_image, (0, 0), mask)
#        logger.info('Watermark image pasted onto the original image.')
#        logger.info('Watermarked image saved to output stream.')
#        output = BytesIO()
#        RGBA_image.save(output, format='PNG')
#        output.seek(0)
#        return output
#
#    except Exception as e:
#        logger.error('An error occurred during the watermarking process.', exc_info=True)
#        raise

def add_watermark(image: BytesIO, watermark_text: str = 'Unknown', bottom: bool = True) -> BytesIO:
    logger.info('Starting the watermarking process.')

    if not bottom:
        normalized_text = normalize_text(watermark_text)
    else:
        normalized_text = watermark_text
    try:
        # Open image and ensure it's in RGBA mode (preserving transparency)
        base_image = Image.open(image).convert('RGBA')
        logger.info('Image loaded and converted to RGBA.')

        width, height = base_image.size
        logger.info(f'Image dimensions: width={width}, height={height}.')

        diagonal = math.sqrt(width**2 + height**2)
        font_size = int(diagonal / 15)
        logger.info(f'Calculated initial font size: {font_size}.')

        try:
            font = ImageFont.truetype(PATH_FONT, font_size)
            logger.info('Loaded Roboto-Regular.ttf font.')
        except IOError:
            logger.warning('Roboto-Regular.ttf not found. Falling back to default font.')
            font = ImageFont.load_default()

        max_text_width = int(width * 0.8)
        min_font_size = 30

        # Adjust font size to fit
        while True:
            draw = ImageDraw.Draw(Image.new("RGBA", (1, 1)))  # Dummy image to get text size
            bbox = draw.textbbox((0, 0), normalized_text, font=font)
            text_width = bbox[2] - bbox[0]

            if text_width <= max_text_width:
                 break

            font_size -= 1
            font = ImageFont.truetype(PATH_FONT, font_size)
            logger.debug(f'Font size reduced to: {font_size}.')

        text_height = bbox[3] - bbox[1]
        text_x = (width - text_width) / 2
        if bottom:
            text_y = height - (2 * text_height)
        else:
            text_y = text_height  # Position near the top
        logger.info(f'Text position calculated: x={text_x}, y={text_y}.')

        # Create a transparent overlay for the watermark
        watermark_layer = Image.new("RGBA", base_image.size, (255, 255, 255, 0))
        draw = ImageDraw.Draw(watermark_layer)
        draw.text((text_x, text_y), normalized_text, font=font, fill=(255, 255, 255, 128))  # White text with transparency

        # Composite the watermark with the original image
        watermarked_image = Image.alpha_composite(base_image, watermark_layer)

        # Save output
        output = BytesIO()
        watermarked_image.save(output, format='PNG')
        output.seek(0)
        return output

    except Exception as e:
        logger.error('An error occurred during the watermarking process.', exc_info=True)
        raise
#def add_watermark(image: BytesIO, watermark_text: str = 'Discord') -> BytesIO:
#    logger.info('Starting the watermarking process.')
#
#    try:
#        # Open image and ensure it's in RGB mode
#        RGB_image = Image.open(image).convert('RGBA')
#        logger.info('Image loaded and converted to RGB.')
#        
#        draw = ImageDraw.Draw(RGB_image)
#        width, height = RGB_image.size
#        logger.info(f'Image dimensions: width={width}, height={height}.')
#        
#        diagonal = math.sqrt(width**2 + height**2)
#        font_size = int(diagonal / 15)
#        logger.info(f'Calculated initial font size: {font_size}.')
#
#        try:
#            font = ImageFont.truetype(PATH_FONT, font_size)
#            logger.info('Loaded Roboto-Regular.ttf font.')
#        except IOError:
#            logger.warning('Roboto-Regular.ttf not found. Falling back to default font.')
#            font = ImageFont.load_default()
#
#        min_font_size = 30
#        max_text_width = 512
#
#        # Adjust font size to fit
#        while True:
#            bbox = draw.textbbox((0, 0), watermark_text, font=font)
#            text_width = bbox[2] - bbox[0]
#
#            if text_width <= max_text_width or font_size <= min_font_size:
#                break
#
#            font_size -= 1
#            font = ImageFont.truetype(PATH_FONT, font_size)
#            logger.debug(f'Font size reduced to: {font_size}.')
#
#        text_height = bbox[3] - bbox[1]
#        text_x = (width - text_width) / 2
#        text_y = height - (2 * text_height)
#        logger.info(f'Text position calculated: x={text_x}, y={text_y}.')
#
#        # Draw watermark directly on RGB image (no transparency)
#        draw.text((text_x, text_y), watermark_text, font=font, fill=(255, 255, 255, 100))
#        logger.info('Watermark text added to the image.')
#
#        # Save output
#        output = BytesIO()
#        RGB_image.save(output, format='PNG')
#        output.seek(0)
#        return output
#
#    except Exception as e:
#        logger.error('An error occurred during the watermarking process.', exc_info=True)
#        raise

def normalize_text(text: str) -> str:
    # Extract only alphabetic characters for the all-uppercase check
    letters_only = ''.join(filter(str.isalpha, text))

    if letters_only.isupper():
        return text
    else:
        return text.lower().capitalize()
