''' combine.py  The purpose of this program is to slap two baddie Image objects next to eachother from from cd ../.
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
from PIL import Image
from utils.add_watermark import add_watermark
from utils.adjust_hue_and_saturation import adjust_hue_and_saturation
from utils.setup_logging import logger

import math

def combine(images: list, names: list) -> BytesIO:
    logger.info('Starting the image combination process.')
    
    combined_images = []
    
    for index, (bytes_io, name) in enumerate(zip(images, names)):
        logger.info(f'Processing image {index + 1}/{len(images)}: {name}')
        
        img = Image.open(bytes_io)
        logger.info(f'Image {name} opened successfully. Size: {img.size}')
        
        inverted_image = Image.eval(img, lambda x: 255 - x)
        logger.info(f'Inverted image {name}.')
        
        img_bytes = BytesIO()
        inverted_image.save(img_bytes, format='PNG')
        logger.info(f'Saved the inverted image {name} to BytesIO.')
        
        img_bytes_final = add_watermark(img_bytes, watermark_text=name)
        logger.info(f'Watermark added to image {name}.')
        
        img_final = Image.open(img_bytes_final)
        combined_images.append(img_final)
        logger.info(f'Image {name} added to the combined images list.')
    
    widths, heights = zip(*(img.size for img in combined_images))
    total_width = sum(widths)
    max_height = max(heights)
    logger.info(f'Combined image dimensions calculated. Total width: {total_width}, Max height: {max_height}.')
    
    combined_img = Image.new('RGB', (total_width, max_height))
    x_offset = 0
    logger.info(f'Created a new blank image for pasting the combined images.')
    
    for img in combined_images:
        combined_img.paste(img, (x_offset, 0))
        x_offset += img.width
        logger.info(f'Pasted an image at offset {x_offset}.')
    
    logger.info('Applying hue and saturation adjustments.')
    output = adjust_hue_and_saturation(combined_img, hue_shift=-180, saturation_shift=160)
    
    output.seek(0)
    logger.info('Image combination process completed successfully.')
    
    return output
