''' google.py  The purpose of this program is to search using the Google Custom Search Restricted API  from cd ../.
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
from bs4 import BeautifulSoup
from utils.setup_logging import logger

import requests

def google(query: str, num_results: int = 5):
    logger.info(f'Starting Google search for query: `{query}` with {num_results} results.')

    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }
    search_url = 'https://www.google.com/search'
    params = {'q': query, 'num': num_results}

    try:
        # Send the search request
        logger.debug(f'Sending request to Google search URL: {search_url} with params: {params}')
        response = requests.get(search_url, headers=headers, params=params)
        response.raise_for_status()  # Raise an error for bad HTTP responses
        logger.info('Received response from Google search.')

        # Parse the search results
        soup = BeautifulSoup(response.text, 'html.parser')
        results = []
        for g in soup.find_all('div', class_='g'):
            title = g.find('h3').text if g.find('h3') else 'No title'
            link = g.find('a')['href'] if g.find('a') else 'No link'
            results.append({'title': title, 'link': link})
            logger.debug(f'Extracted result: Title: {title}, Link: {link}')

            if len(results) >= num_results:
                break

        logger.info(f'Successfully extracted {len(results)} search results.')
        return results

    except requests.exceptions.RequestException as e:
        logger.error(f'Error during the web request: {e}')
        return []

    except Exception as e:
        logger.error(f'An unexpected error occurred: {e}')
        return []
