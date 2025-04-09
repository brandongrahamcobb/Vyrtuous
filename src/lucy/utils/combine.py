''' combine.py  The purpose of this program is to slap two baddie Image objects next to eachother.
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
from PIL import Image, ImageDraw, ImageEnhance, ImageFont
from itertools import zip_longest
from lucy.utils.add_watermark import add_watermark
from lucy.utils.adjust_hue_and_saturation import adjust_hue_and_saturation
from lucy.utils.setup_logging import logger
from math import ceil, sqrt

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
