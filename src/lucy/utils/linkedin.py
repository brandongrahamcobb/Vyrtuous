''' linkedin.py   endpoint using their python SDK is
                                much more efficient than this program. This is a complicated
                                way to get a completion from OpenAI from cd ../.
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

from collections import defaultdict
from PIL import Image
from utils.create_https_completion import create_https_completion
from utils.load_yaml import load_yaml
from utils.setup_logging import logger

import asyncio
import io
import requests
import shlex
import time
import traceback
import utils.helpers as helpers
import yaml


class LinkedInBot:
    def __init__(self, config_path, conversation_log="training.jsonl"):
        self.config = load_yaml(config_path)
        self.access_token = self.config['api_keys'].get('api_key_2').get('api_key')
        self.base_url = "https://api.linkedin.com/v2"
        self.headers = {
            "Authorization": f"Bearer {self.access_token}" if self.access_token else "",
            "Content-Type": "application/json"
        }
        self.conversations = defaultdict(list)
        self.conversation_log = conversation_log
        self.profile_urn = self.get_profile_urn()

    def get_authorization_url(self):
        base_url = "https://www.linkedin.com/oauth/v2/authorization"
        redirect_uri = self.config['api_keys']['api_key_2']['redirect_uri']
        client_id = self.config['api_keys']['api_key_2']['client_id']
        scopes = "profile w_member_social"
        return (
            f"{base_url}?response_type=code"
            f"&client_id={client_id}"
            f"&redirect_uri={redirect_uri}"
            f"&scope={scopes}"
        )

    def exchange_authorization_code(self, authorization_code):
        url = "https://www.linkedin.com/oauth/v2/accessToken"
        payload = {
            "grant_type": "authorization_code",
            "code": authorization_code,
            "redirect_uri": self.config['api_keys']['api_key_2']['redirect_uri'],
            "client_id": self.config['api_keys']['api_key_2']['client_id'],
            "client_secret": self.config['api_keys']['api_key_2']['client_secret']
        }
        headers = {"Content-Type": "application/x-www-form-urlencoded"}
        response = requests.post(url, data=payload, headers=headers)
        if response.status_code == 200:
            token_data = response.json()
            self.access_token = token_data.get("access_token")
            self.headers["Authorization"] = f"Bearer {self.access_token}"
            print("Access token acquired successfully!")
        else:
            print(f"Error exchanging authorization code: {response.status_code}, {response.text}")

    def get_profile_urn(self):
        url = f"{self.base_url}/me"
        response = requests.get(url, headers=self.headers)
        if response.status_code == 200:
            profile_data = response.json()
            return profile_data.get("id")
        elif response.status_code == 403:
            print("Error: Insufficient permissions. Ensure 'profile' scope is included.")
        else:
            print(f"Error fetching profile URN: {response.status_code}, {response.text}")
        return None

    # Remaining methods are unchanged (create_first_post, respond_to_messages, etc.)

    async def main(self):
        if not self.access_token or not self.profile_urn:
            print("Authorization is required. Generate an authorization URL:")
            print(self.get_authorization_url())
            return

        await self.create_first_post()
        while True:
            self.respond_to_messages()
            time.sleep(60)


if __name__ == "__main__":
    CONFIG_PATH = helpers.PATH_CONFIG_YAML
    bot = LinkedInBot(CONFIG_PATH)
    asyncio.run(bot.main())
