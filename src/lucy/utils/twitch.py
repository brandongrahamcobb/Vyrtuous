''' twitch.py The purpose of this program is to handle all twitch bot related functions.
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
# twitch.py

import aiohttp
from quart import Quart, request, session, redirect
from datetime import datetime, timedelta
import os
import yaml

# Import your logging utility
# e.g., from utils.setup_logging import setup_logging
# and initialize the logger for this file (using __name__ helps identify the source)
# logger = setup_logging(__name__)

# If your utils.setup_logging directly provides a logger, you might do:
from utils.create_https_moderation import create_https_moderation
from utils.create_https_moderation import create_https_moderation
from utils.load_yaml import load_yaml
from utils.setup_logging import logger
import utils.helpers as helpers

# Load your configuration YAML or however you handle config
CONFIG = load_yaml(helpers.PATH_CONFIG_YAML)

app = Quart(__name__)
app.secret_key = CONFIG['twitch_token']

# Twitch OAuth constants
CLIENT_ID = 'l4dmn34apx38yj70jcwchqihuaieby'
CLIENT_SECRET = '1dzuy41ztp2ybfcqwouy2pz9zxte0v'
REDIRECT_URI = 'http://localhost:5000/callback'
TOKEN_URL = 'https://id.twitch.tv/oauth2/token'
API_URL = 'https://api.twitch.tv/helix/'
AUTH_URL = 'https://id.twitch.tv/oauth2/authorize'
SCOPES = ['chat:read', 'chat:edit']

@app.route("/authorize")
async def authorize():
    """
    Build the Twitch OAuth authorization URL and redirect the user.
    """
    logger.info("Starting the authorization flow...")

    # Generate a random state token to mitigate CSRF
    state = os.urandom(16).hex()
    session["oauth_state"] = state

    # Build the Twitch OAuth authorization URL
    params = {
        "client_id": CLIENT_ID,
        "redirect_uri": REDIRECT_URI,
        "response_type": "code",
        "scope": " ".join(SCOPES),  # or use '+'. e.g.: "chat:read+chat:edit"
        "state": state
    }

    # Construct the full URL
    url = (
        f"{AUTH_URL}"
        f"?client_id={params['client_id']}"
        f"&redirect_uri={params['redirect_uri']}"
        f"&response_type={params['response_type']}"
        f"&scope={params['scope'].replace(' ', '+')}"
        f"&state={params['state']}"
    )

    logger.debug(f"OAuth URL constructed: {url}")

    # Redirect user to Twitch for authentication
    return redirect(url)


@app.route("/callback", methods=['GET'])
async def oauth_callback():
    """
    Handles the redirect from Twitch, with 'code' in query parameters.
    Exchanges the code for an access token that belongs to the user.
    """
    logger.info("Received Twitch OAuth callback.")

    # Check for state mismatch (basic security check)
    if request.args.get("state") != session.get("oauth_state"):
        logger.warning("Invalid state token: possible CSRF detected.")
        return "Invalid state token", 400

    # Get the 'code' from the query parameters
    code = request.args.get("code")
    if not code:
        logger.error("Missing authorization code in callback.")
        return "Missing authorization code", 400

    logger.debug(f"Authorization code received: {code}")

    # Exchange the authorization code for a user access token
    data = {
        "client_id": CLIENT_ID,
        "client_secret": CLIENT_SECRET,
        "code": code,
        "grant_type": "authorization_code",
        "redirect_uri": REDIRECT_URI,
    }

    async with aiohttp.ClientSession() as http_session:
        async with http_session.post(TOKEN_URL, data=data) as resp:
            token_data = await resp.json()

    if "access_token" not in token_data:
        logger.error(f"Token exchange failed: {token_data}")
        return f"Token exchange failed: {token_data}", 400

    # At this point, we have a user-based token
    user_access_token = token_data["access_token"]
    refresh_token = token_data.get("refresh_token")
    expires_in = token_data["expires_in"]

    # Store tokens in session (or DB, as needed)
    session["user_access_token"] = user_access_token
    session["refresh_token"] = refresh_token
    expiration_time = datetime.utcnow() + timedelta(seconds=expires_in)
    session["expires_at"] = expiration_time.isoformat()

    logger.info("Successfully exchanged authorization code for an access token.")

    # (Optional) Get the user’s Twitch info to confirm identity
    async with aiohttp.ClientSession() as http_session:
        headers = {
            "Client-ID": CLIENT_ID,
            "Authorization": f"Bearer {user_access_token}"
        }
        async with http_session.get("https://api.twitch.tv/helix/users", headers=headers) as user_resp:
            user_data = await user_resp.json()

    if not user_data.get("data"):
        logger.error("Could not fetch user data from Twitch.")
        return "Could not fetch user data", 400

    user = user_data["data"][0]
    user_id = user["id"]
    user_login = user["login"]
    display_name = user["display_name"]

    logger.info(
        f"Got user token for {display_name} (ID {user_id}). "
        f"Expires in {expires_in} seconds."
    )

    # (Optional) Store user info in session
    session["user_id"] = user_id
    session["user_login"] = user_login
    session["display_name"] = display_name

    # Here, you can redirect them to a "success" page, or run your bot, etc.
    with open("token.txt", "w") as f:
        f.write(user_access_token)

    return f"Success! Token acquired. You are {display_name}."


# ---------------------------------------------------------
# BOT CODE
# ---------------------------------------------------------
from discord.ext import commands as discord_commands
from twitchio.ext import commands as twitch_commands
from utils.create_https_completion import Conversations
import json

class Vyrtuous(twitch_commands.Bot):
    def __init__(self, bot: discord_commands.Bot, access_token: str):
        super().__init__(
            token=access_token,
            client_id=CLIENT_ID,
            nick='Lucy_',
            prefix="!",
            initial_channels=['spawdspawd']
        )
        self.bot = bot
        self.conversations = Conversations()
        self.config = CONFIG

    async def event_ready(self):
        logger.info("Hello World!")
        logger.info(f"Bot is ready! Logged in as Lucy_")
        logger.info(f"Connected to channel: spawdspawd")

    async def event_message(self, message):
        """Handle every message in the Twitch chat."""
        # You might log incoming messages at info or debug level
        logger.info(f"Received message: {message.content}")

        if message.author.name.lower() == 'Lucy_':
            # Ignore the bot's own messages
            return

        logger.info(f"Message from {message.author.name}: {message.content}")

        # Prepare OpenAI API request
        array = []
        input_text_dict = {
            'type': 'text',
            'text': message.content
        }
        array.append(input_text_dict)

        # Make your OpenAI API call here. This is hypothetical code.
        async for response in self.conversations.create_https_completion(
            completions=self.config['openai_chat_n'],
            custom_id=message.author.id,
            input_array=array,
            max_tokens=self.config['openai_chat_max_tokens'],
            model=self.config['openai_chat_model'],
            response_format=self.config['openai_chat_response_format'],
            stop=self.config['openai_chat_stop'],
            store=self.config['openai_chat_store'],
            stream=self.config['openai_chat_stream'],
            sys_input=self.config['openai_chat_sys_input'],
            temperature=self.config['openai_chat_temperature'],
            top_p=self.config['openai_chat_top_p'],
            use_history=self.config['openai_chat_use_history'],
            add_completion_to_history=self.config['openai_chat_add_completion_to_history']
        ):
            await message.channel.send(response)
            logger.debug(f"Sent message: {response}")
#        channel = await self.bot.fetch_channel(1315735859848544378)
 #       await channel.send(message.content)

        async for moderation in create_https_moderation(message.author.id, array, model=helpers.OPENAI_MODERATION_MODEL):
            results = moderation.get('results', [])
            if results and results[0].get('flagged', False):
                await message.delete()

        if self.config['openai_chat_moderation']:
            async for moderation in self.conversations.create_https_completion(
                completions=helpers.OPENAI_CHAT_MODERATION_N,
                custom_id=message.author.id,
                input_array=array,
                max_tokens=helpers.OPENAI_CHAT_MODERATION_MAX_TOKENS,
                model=helpers.OPENAI_CHAT_MODERATION_MODEL,
                response_format=helpers.OPENAI_CHAT_MODERATION_RESPONSE_FORMAT,
                stop=helpers.OPENAI_CHAT_MODERATION_STOP,
                store=helpers.OPENAI_CHAT_MODERATION_STORE,
                stream=helpers.OPENAI_CHAT_MODERATION_STREAM,
                sys_input=helpers.OPENAI_CHAT_MODERATION_SYS_INPUT,
                temperature=helpers.OPENAI_CHAT_MODERATION_TEMPERATURE,
                top_p=helpers.OPENAI_CHAT_MODERATION_TOP_P,
                use_history=helpers.OPENAI_CHAT_MODERATION_USE_HISTORY,
                add_completion_to_history=helpers.OPENAI_CHAT_MODERATION_ADD_COMPLETION_TO_HISTORY
             ):
                full_response = json.loads(moderation)
                results = full_response.get('results', [])
                flagged = results[0].get('flagged', False)
                carnism_flagged = results[0]['categories'].get('carnism', False)
                carnism_score = results[0]['category_scores'].get('carnism', 0)
                total_carnism_score = sum(arg['category_scores'].get('carnism', 0) for arg in results)
                if carnism_flagged or flagged:  # If carnism is flagged
                    if not self.config['discord_role_pass']:
                        await message.delete()
                    NLPUtils.append_to_other_jsonl('training.jsonl', carnism_score, message.content, message.author.id) #results[0].get('flagged', False), message.content)
