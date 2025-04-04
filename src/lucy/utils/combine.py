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
from PIL import Image, ImageEnhance
from itertools import zip_longest
from lucy.utils.add_watermark import add_watermark
from lucy.utils.adjust_hue_and_saturation import adjust_hue_and_saturation
from lucy.utils.setup_logging import logger
from math import ceil, sqrt

def combine_gallery(images: list, names: list, title: str, quantity: int = 1, linearity: bool = False) -> BytesIO:
    if len(images) <= 2:
        linearity = True

    processed_images = []

    # Repeat images based on quantity
    repeated_images = list(images * quantity)  # Ensuring the quantity requirement is met

    # Process each image and store it in the processed_images list
    for img_bytes in repeated_images:
        img_bytes.seek(0)  # Reset BytesIO pointer
        img = Image.open(img_bytes).convert("RGBA")  # Open and convert to RGBA
        processed_images.append(img)

    if linearity == True:
        # Create a linear (horizontal) arrangement of images
        combined_width = sum(img.size[0] for img in processed_images)  # Total width
        combined_height = max(img.size[1] for img in processed_images)  # Max height

        combined_img = Image.new('RGBA', (combined_width, combined_height), (255, 255, 255, 0))  # Transparent background
        x_offset = 0

        for img in processed_images:
            combined_img.paste(img, (x_offset, 0), img)  # Use alpha mask for transparency
            x_offset += img.size[0]  # Update horizontal offset

    else:
        # Create a grid layout
        num_images = len(processed_images)
        grid_size = ceil(sqrt(num_images))  # Determine grid size
        image_size = max(img.size[0] for img in processed_images)  # Use largest width
        canvas_size = grid_size * image_size  # Square canvas
        footer_ratio = 0.15  # Footer is 15% of the square canvas
        footer_height = int(canvas_size * footer_ratio)
        new_height = canvas_size + footer_height


        combined_img = Image.new('RGB', (canvas_size, new_height), (0, 0, 0))  # Transparent background

        for index, img in enumerate(processed_images):
            row = index // grid_size
            col = index % grid_size
            x_offset = col * image_size
            y_offset = row * image_size

            img = img.resize((image_size, image_size), Image.LANCZOS)
            combined_img.paste(img, (x_offset, y_offset), img)  # Use alpha mask for transparency

    # Save image to buffer
    combined_img_buffer = BytesIO()
    combined_img.save(combined_img_buffer, format='PNG')
    combined_img_buffer.seek(0)

    # Apply watermark if multiple images are combined
    if len(images) > 1:
        final_image_buffer = add_watermark(combined_img_buffer, title, True)
    else:
        final_image_buffer = combined_img_buffer

    final_image_buffer.seek(0)
    return final_image_buffer

#def combine_gallery(images: list, names: list, title: str, quantity: int = 1, linear: bool = False) -> BytesIO:
#    if linear:
#        
#    num_images = len(images)*quantity
#    grid_size = ceil(sqrt(num_images))  # Determine grid size
#    processed_images = []
#
#    for img_bytes, name in zip_longest(images, names, fillvalue=''):
#        img_bytes.seek(0)  # Reset BytesIO pointer
#        img = Image.open(img_bytes).convert("RGBA")  # Open and convert to RGBA
###        inverted_img = Image.eval(img, lambda x: 255 - x)  # Invert colors
##        adjusted_output = adjust_hue_and_saturation(inverted_img, hue_shift=-180, saturation_shift=160)  
##        adjusted_output.seek(0)
##        watermarked_image_buffer = add_watermark(adjusted_output, name)
##        watermarked_image = Image.open(watermarked_image_buffer).convert("RGBA")
##        processed_images.append(watermarked_image)
#        processed_images.append(img)
#
#    # Determine square dimensions
#    image_size = max(img.size[0] for img in processed_images)  # Use the largest width
#    canvas_size = grid_size * image_size  # Create square canvas
#
#    # Create a blank white canvas with RGBA mode (to support transparency)
#    combined_img = Image.new('RGBA', (canvas_size, canvas_size), (255, 255, 255, 0))  # Transparent background
#
#    for index, img in enumerate(processed_images):
#        row = index // grid_size
#        col = index % grid_size
#        x_offset = col * image_size
#        y_offset = row * image_size
#
#        # Ensure consistent size before pasting
#        img = img.resize((image_size, image_size), Image.LANCZOS)
#
#        # Paste the image onto the transparent background
#        combined_img.paste(img, (x_offset, y_offset), img)  # Use alpha mask for transparency
#
#    combined_img_buffer = BytesIO()
#    combined_img.save(combined_img_buffer, format='PNG')
#    combined_img_buffer.seek(0)
#    if len(images) != 1:
#        final_image_buffer = add_watermark(combined_img_buffer, title)
#        final_image_buffer.seek(0)
#    else:
#        final_image_buffer = combined_img_buffer
#        final_image_buffer.seek(0)
#    return final_image_buffer
