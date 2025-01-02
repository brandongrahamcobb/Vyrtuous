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
import aiohttp
from quart import Quart, request, session, redirect
from datetime import datetime, timedelta
import os
import yaml

# Import your logging utility
from .create_https_moderation import create_https_moderation
from .load_yaml import load_yaml
from .setup_logging import logger
from .helpers import *

# Load your configuration YAML
CONFIG = load_yaml(PATH_CONFIG_YAML)

app = Quart(__name__)
app.secret_key = os.urandom(24)  # Use a secure random key

# Twitch OAuth constants
CLIENT_ID = 'l4dmn34apx38yj70jcwchqihuaieby'
CLIENT_SECRET = '1dzuy41ztp2ybfcqwouy2pz9zxte0v'
REDIRECT_URI = 'http://localhost:5000/callback'
TOKEN_URL = 'https://id.twitch.tv/oauth2/token'
API_URL = 'https://api.twitch.tv/helix/'
SCOPES = ['chat:read', 'chat:edit']

# In-memory token storage (consider using persistent storage)
oauth_data = {
    "user_access_token": None,
    "refresh_token": None,
    "expires_at": None,
}

async def exchange_token(code):
    """
    Exchange the authorization code for an access token.
    """
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
        return None

    # Update in-memory token storage
    oauth_data["user_access_token"] = token_data["access_token"]
    oauth_data["refresh_token"] = token_data.get("refresh_token")
    oauth_data["expires_at"] = datetime.utcnow() + timedelta(seconds=token_data["expires_in"])

    logger.info("Successfully exchanged authorization code for an access token.")
    return token_data

async def refresh_token():
    """
    Refresh the OAuth token if expired.
    """
    if not oauth_data["refresh_token"]:
        logger.error("No refresh token available.")
        return None

    data = {
        "client_id": CLIENT_ID,
        "client_secret": CLIENT_SECRET,
        "refresh_token": oauth_data["refresh_token"],
        "grant_type": "refresh_token",
    }

    async with aiohttp.ClientSession() as http_session:
        async with http_session.post(TOKEN_URL, data=data) as resp:
            token_data = await resp.json()

    if "access_token" not in token_data:
        logger.error(f"Token refresh failed: {token_data}")
        return None

    # Update in-memory token storage
    oauth_data["user_access_token"] = token_data["access_token"]
    oauth_data["refresh_token"] = token_data.get("refresh_token")
    oauth_data["expires_at"] = datetime.utcnow() + timedelta(seconds=token_data["expires_in"])

    logger.info("Token refreshed successfully.")
    return token_data

async def ensure_token():
    """
    Ensure a valid token is available (exchange or refresh as needed).
    """
    if oauth_data["user_access_token"] and datetime.utcnow() < oauth_data["expires_at"]:
        logger.debug(f"Valid token found: {oauth_data['user_access_token']}")
        return oauth_data["user_access_token"]

    # If no token or expired token, initiate refresh or authorization
    if oauth_data["refresh_token"]:
        logger.debug("Attempting to refresh token...")
        refreshed_token = await refresh_token()
        if refreshed_token:
            return oauth_data["user_access_token"]
        else:
            logger.error("Token refresh failed.")
            return None

    logger.error("No valid token available; authorization is required.")
    return None

@app.route("/authorize")
async def authorize():
    """
    Build the Twitch OAuth authorization URL and automatically fetch the token.
    """
    logger.info("Starting the automated authorization flow...")

    # Construct the authorization URL
    auth_url = (
        f"https://id.twitch.tv/oauth2/authorize"
        f"?client_id={CLIENT_ID}"
        f"&redirect_uri={REDIRECT_URI}"
        f"&response_type=code"
        f"&scope={' '.join(SCOPES)}"
    )
    logger.debug(f"Redirecting to Twitch OAuth URL: {auth_url}")
    return redirect(auth_url)

@app.route("/callback")
async def oauth_callback():
    """
    Handles the redirect from Twitch and exchanges the authorization code for an access token.
    """
    code = request.args.get("code")
    if not code:
        logger.error("Missing authorization code in callback.")
        return "Missing authorization code", 400

    logger.debug(f"Authorization code received: {code}")
    token_data = await exchange_token(code)
    if not token_data:
        return "Token exchange failed", 400

    # Optionally, you can initialize the Twitch bot here or signal that the token is ready
    return "Token exchange complete! You can now start using the bot."

# Example endpoint to demonstrate token usage
@app.route("/validate_token")
async def validate_token():
    """
    Validate and return the current access token.
    """
    token = await ensure_token()
    if not token:
        return "No valid token available. Reauthorization required.", 401

    return f"Access Token: {token}"

# Initialize the Twitch bot after ensuring the token is available
@app.before_serving
async def startup():
    """
    Handle startup tasks, such as initializing the Twitch bot.
    """
    # Optionally, you can check if a token exists and is valid
    token = await ensure_token()
    if not token:
        logger.warning("No valid token available on startup. Please authorize the application.")
        return

    # Initialize your Twitch bot here with the valid token
    # Example:
    twitch_bot = Vyrtuous(bot, token)
    await twitch_bot.start()

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)

# ---------------------------------------------------------
# BOT CODE
# ---------------------------------------------------------
from discord.ext import commands as discord_commands
from twitchio.ext import commands as twitch_commands
from .create_https_completion import Conversations
import json

class Vyrtuous(twitch_commands.Bot):
    def __init__(self, bot: discord_commands.Bot, access_token):
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

        async for moderation in create_https_moderation(message.author.id, array, model=OPENAI_MODERATION_MODEL):
            results = moderation.get('results', [])
            if results and results[0].get('flagged', False):
                await message.delete()

        if self.config['openai_chat_moderation']:
            async for moderation in self.conversations.create_https_completion(
                completions=OPENAI_CHAT_MODERATION_N,
                custom_id=message.author.id,
                input_array=array,
                max_tokens=OPENAI_CHAT_MODERATION_MAX_TOKENS,
                model=OPENAI_CHAT_MODERATION_MODEL,
                response_format=OPENAI_CHAT_MODERATION_RESPONSE_FORMAT,
                stop=OPENAI_CHAT_MODERATION_STOP,
                store=OPENAI_CHAT_MODERATION_STORE,
                stream=OPENAI_CHAT_MODERATION_STREAM,
                sys_input=OPENAI_CHAT_MODERATION_SYS_INPUT,
                temperature=OPENAI_CHAT_MODERATION_TEMPERATURE,
                top_p=OPENAI_CHAT_MODERATION_TOP_P,
                use_history=OPENAI_CHAT_MODERATION_USE_HISTORY,
                add_completion_to_history=OPENAI_CHAT_MODERATION_ADD_COMPLETION_TO_HISTORY
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
