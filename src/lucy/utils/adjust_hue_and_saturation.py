''' adjust_hue_and_saturation.py  The purpose of this program is reverse the colors of a PIL Image objects from cd ../.
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
from utils.setup_logging import logger

import colorsys

def adjust_hue_and_saturation(image, hue_shift, saturation_shift) -> BytesIO:
    try:
        logger.info("Starting hue and saturation adjustment.")
        logger.debug(f"Hue shift: {hue_shift}, Saturation shift: {saturation_shift}")

        # Convert image to RGB mode
        image = image.convert('RGB')
        logger.debug("Image converted to RGB mode.")

        # Get pixel data and apply adjustments
        pixels = list(image.getdata())
        adjusted_pixels = []
        logger.debug(f"Total pixels to process: {len(pixels)}")

        for r, g, b in pixels:
            h, s, v = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
            h = (h + hue_shift / 360.0) % 1.0
            s = min(max(s + saturation_shift / 100.0, 0), 1)
            r, g, b = colorsys.hsv_to_rgb(h, s, v)
            adjusted_pixels.append((int(r * 255), int(g * 255), int(b * 255)))

        logger.info("Pixel adjustments completed.")

        # Create a new image with adjusted pixels
        new_image = Image.new('RGB', image.size)
        new_image.putdata(adjusted_pixels)
        logger.debug("New image created with adjusted pixels.")

        # Save the new image to a BytesIO object
        output = BytesIO()
        new_image.save(output, format='PNG')
        output.seek(0)
        logger.info("Image saved to BytesIO object successfully.")

        return output
    except Exception as e:
        logger.error(f"An error occurred during hue and saturation adjustment: {e}")
        raise
