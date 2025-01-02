''' gsrs.py  The purpose of this program is to fetch a GSRS molecule image from from cd ../.
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
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from utils.add_watermark import add_watermark
from utils.setup_logging import logger
from webdriver_manager.chrome import ChromeDriverManager

import math
import os
import requests

path = os.path.dirname(os.path.abspath(__file__))

def gsrs(arg):
    logger.info(f'Starting GSRS search for argument: {arg}')

    # Set up Chrome options for headless operation
    chrome_options = Options()
    chrome_options.add_argument('--headless')  # Run headless Chrome (no UI)
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')

    executable_path = os.path.join(path, 'chromedriver')
    logger.debug(f'Using executable path: {executable_path}')

    driver = webdriver.Chrome(service=Service(executable_path=executable_path), options=chrome_options)
    try:
        # Navigate to the GSRS URL
        search_url = f'https://gsrs.ncats.nih.gov/ginas/app/beta/browse-substance?search={arg}'
        logger.info(f'Navigating to URL: {search_url}')
        driver.get(search_url)
        driver.implicitly_wait(10)  # Adjust the wait time as needed

        # Locate the image element
        logger.debug('Searching for the image element.')
        img_element = driver.find_element(By.CSS_SELECTOR, 'body > app-root > app-base > app-substances-browse > div > div.substance-cards > app-substance-summary-card > mat-card > mat-card-title > a')

        if img_element:
            logger.info('Image element found.')
            img_src = img_element.get_attribute('href')

            if img_src:
                logger.info(f'Image source URL retrieved: {img_src}')

                stripped = img_src.split('/', -1)[-1:]
                link = f'https://gsrs.ncats.nih.gov/api/v1/substances/render({stripped[0]})?format=png&size=512&stereo=true'
                logger.debug(f'Constructed image link: {link}')

                # Fetch the image
                response = requests.get(link)
                response.raise_for_status()
                logger.info('Image successfully fetched from GSRS.')

                image_bytes = response.content
                image = Image.open(BytesIO(image_bytes)).convert('RGBA')
                draw = ImageDraw.Draw(image)

                # Add watermark
                width, height = image.size
                diagonal = math.sqrt(width**2 + height**2)
                font_size = int(diagonal / 15)
                logger.debug(f'Calculated font size for watermark: {font_size}')

                try:
                    font = ImageFont.truetype('Roboto-Regular.ttf', font_size)
                    logger.info('Loaded custom font for watermark.')
                except IOError:
                    font = ImageFont.load_default()
                    logger.warning('Custom font not found. Using default font.')

                bbox = draw.textbbox((0, 0), arg, font=font)
                text_width = bbox[2] - bbox[0]
                text_height = bbox[3] - bbox[1]
                text_x = (width - text_width) / 2
                text_y = (height - text_height) / 2

                watermark_image = Image.new('RGBA', image.size, (0, 0, 0, 0))
                watermark_draw = ImageDraw.Draw(watermark_image)
                watermark_draw.text((text_x, text_y), arg, font=font, fill=(255, 255, 255, 64))
                mask = watermark_image.split()[3]
                image.paste(watermark_image, (0, 0), mask)

                logger.info('Watermark successfully added to the image.')
                return image
            else:
                logger.warning('No `src` attribute found in the <img> element.')
                return 'No src attribute found in the <img> element'
        else:
            logger.warning('No <img> element found with the specified CSS path.')
            return 'No <img> element found with the specified CSS path'
    except Exception as e:
        logger.error(f'An error occurred during the GSRS process: {e}')
        raise
    finally:
        driver.quit()
        logger.info('WebDriver session closed.')
