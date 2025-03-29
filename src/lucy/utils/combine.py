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

def combine_gallery(images: list, names: list, title: str) -> BytesIO:
    num_images = len(images)
    grid_size = ceil(sqrt(num_images))  # Determine grid size
    processed_images = []

    for img_bytes, name in zip_longest(images, names, fillvalue=''):
        img_bytes.seek(0)  # Reset BytesIO pointer
        img = Image.open(img_bytes).convert("RGBA")  # Open and convert to RGBA
##        inverted_img = Image.eval(img, lambda x: 255 - x)  # Invert colors
#        adjusted_output = adjust_hue_and_saturation(inverted_img, hue_shift=-180, saturation_shift=160)  
#        adjusted_output.seek(0)
#        watermarked_image_buffer = add_watermark(adjusted_output, name)
#        watermarked_image = Image.open(watermarked_image_buffer).convert("RGBA")
#        processed_images.append(watermarked_image)
        processed_images.append(img)

    # Determine square dimensions
    image_size = max(img.size[0] for img in processed_images)  # Use the largest width
    canvas_size = grid_size * image_size  # Create square canvas

    # Create a blank white canvas with RGBA mode (to support transparency)
    combined_img = Image.new('RGBA', (canvas_size, canvas_size), (255, 255, 255, 0))  # Transparent background

    for index, img in enumerate(processed_images):
        row = index // grid_size
        col = index % grid_size
        x_offset = col * image_size
        y_offset = row * image_size

        # Ensure consistent size before pasting
        img = img.resize((image_size, image_size), Image.LANCZOS)

        # Paste the image onto the transparent background
        combined_img.paste(img, (x_offset, y_offset), img)  # Use alpha mask for transparency

    combined_img_buffer = BytesIO()
    combined_img.save(combined_img_buffer, format='PNG')
    combined_img_buffer.seek(0)
    if len(images) != 1:
        final_image_buffer = add_watermark(combined_img_buffer, title)
        final_image_buffer.seek(0)
    else:
        final_image_buffer = combined_img_buffer
        final_image_buffer.seek(0)
    return final_image_buffer
#from itertools import zip_longest
#from lucy.utils.add_watermark import add_watermark
#from lucy.utils.adjust_hue_and_saturation import adjust_hue_and_saturation
#from lucy.utils.setup_logging import logger
#from PIL import Image
#from io import BytesIO
#from math import ceil, sqrt
##
#def combine_gallery(images: list, names: list, title: str) -> BytesIO:
#    num_images = len(images)
#    grid_size = ceil(sqrt(num_images))  # Determine grid size
#    processed_images = []
#
#    # Load and ensure all images are valid
#    # Load and ensure all images are valid
##    for img_bytes in images:
##        img_bytes.seek(0)  # Reset BytesIO pointer
##        img = Image.open(img_bytes).convert("RGBA")  # Ensure RGBA mode
##        processed_images.append(img)
#    for img_bytes, name in zip_longest(images, names, fillvalue=''):
#        img_bytes.seek(0)  # Reset BytesIO pointer
#        inverted_img = Image.eval(img_bytes, lambda x: 255 - x)
#        adjusted_output = adjust_hue_and_saturation(inverted_img, hue_shift=-180, saturation_shift=160)  
#        adjusted_output.seek(0)
#        watermarked_image_buffer = add_watermark(adjusted_output, name)
#        watermarked_image = Image.open(watermarked_image_buffer).convert("RGBA")
#        processed_images.append(watermarked_image)
#
#    # Determine square dimensions
#    image_size = max(img.size[0] for img in processed_images)  # Use the largest width
#    canvas_size = grid_size * image_size  # Create square canvas
#
#    # Create a blank white canvas with RGB mode (no transparency)
#    combined_img = Image.new('RGB', (canvas_size, canvas_size), (255, 255, 255))  # White background
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
#        # Paste the image onto the white background (no transparency mask)
#        combined_img.paste(img, (x_offset, y_offset))
#
#    combined_img_buffer = BytesIO()
#    combined_img.save(combined_img_buffer, format='PNG')
#    combined_img_buffer.seek(0)
#    final_image_buffer = add_watermark(combined_img_buffer, title)
#    final_image_buffer.seek(0)
#    return final_image_buffer
#
#    # Save the final combined image to a BytesIO object
#    output = BytesIO()
#    combined_img.save(output, format='PNG')
#    output.seek(0)
#
#    pil_image = Image.open(output).convert("RGB")
#
#    # Apply inversion
#    inverted_img = Image.eval(pil_image, lambda x: 255 - x)
#
#    adjusted_output = adjust_hue_and_saturation(inverted_img, hue_shift=-180, saturation_shift=160)  
#    adjusted_output.seek(0)
#
#    # Convert adjusted output back to PIL Image for inversion
#    adjusted_img = Image.open(adjusted_output).convert("RGB")
#
#    # Save final inverted image
#    final_output = BytesIO()
#    adjusted_img.save(final_output, format='PNG')
#    final_output.seek(0)
#    return final_output

#def addwatermark(image: Image.Image, text: str) -> Image.Image:
#    """Adds a watermark with the given text to an image."""
#    draw = ImageDraw.Draw(image)
#    font = ImageFont.load_default()  # Load a default font
#    text_size = draw.textsize(text, font=font)
#    position = (image.width - text_size[0] - 10, image.height - text_size[1] - 10)  # Bottom-right corner
#    draw.text(position, text, fill=(255, 255, 255, 128), font=font)  # White text with transparency
#    return image

#def combine_gallery(images: list, names: list, title: str) -> BytesIO:
#    num_images = len(images)
#    grid_size = ceil(sqrt(num_images))  # Determine grid size
#    processed_images = []
#
#    # Load and ensure all images are valid
#    for img_bytes in images:
#        img_bytes.seek(0)  # Reset BytesIO pointer
#        img = Image.open(img_bytes).convert("RGBA")  # Ensure RGBA mode
#        img = addwatermark(img, name) 
#        processed_images.append(img)
#
#    # Determine square dimensions
#    image_size = max(img.size[0] for img in processed_images)  # Use the largest width
#    canvas_size = grid_size * image_size  # Create square canvas
#
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
#        combined_img.paste(img, (x_offset, y_offset), img)  # Use alpha mask for transparency
#
#    output = BytesIO()
#    inverted_img = Image.eval(combined_img, lambda x: 255 - x)
#    inverted_img.save(output, format='PNG')
#    output = adjust_hue_and_saturation(inverted_img, hue_shift=-180, saturation_shift=160)
#    output.seek(0)
#    return output

#def combine_gallery(images: list, names: list, title: str) -> BytesIO:
#    num_images = len(images)
#    grid_size = ceil(sqrt(num_images))  # Determine grid size
#    processed_images = []
#
#    # Load and process individual images
#    for img_bytes, name in zip(images, names):
#        img_bytes.seek(0)  # Reset BytesIO pointer
#        img = Image.open(img_bytes).convert("RGBA")  # Ensure RGBA mode
#
#        # Convert back to BytesIO and add watermark using add_watermark
#        img_buffer = BytesIO()
#        img.save(img_buffer, format='PNG')
#        img_buffer.seek(0)
#        watermarked_img = add_watermark(img_buffer, name)  # Use existing add_watermark
#        processed_images.append(Image.open(watermarked_img))
#
#    # Determine square dimensions
#    image_size = max(img.size[0] for img in processed_images)  # Use the largest width
#    canvas_size = grid_size * image_size  # Create square canvas
#
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
#        combined_img.paste(img, (x_offset, y_offset), img)  # Use alpha mask for transparency
#
#    # Convert combined image to BytesIO and add watermark with title
#    combined_img_buffer = BytesIO()
#    combined_img.save(combined_img_buffer, format='PNG')
#    combined_img_buffer.seek(0)
#    inverted_img = Image.eval(combined_img, lambda x: 255 - x)
#    adjusted_output = adjust_hue_and_saturation(inverted_img, hue_shift=-180, saturation_shift=160)  
#    adjusted_output.seek(0)
#    final_image_buffer = add_watermark(adjusted_output, title)
#    final_image_buffer.seek(0)
#    return final_image_buffer
