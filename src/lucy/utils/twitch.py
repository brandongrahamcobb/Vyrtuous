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

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)

# ---------------------------------------------------------
# BOT CODE
# ---------------------------------------------------------
from discord.ext import commands as discord_commands
from twitchio.ext import commands as twitch_commands
from .create_https_completion import Conversations
from .message import Message
import json

class Vyrtuous(twitch_commands.Bot):
    def __init__(self, config, conversations, db_pool, lock, token):
        super().__init__(
            token=token,
            client_id=CLIENT_ID,
            nick='Lucy_',
            prefix="!",
            initial_channels=['spawdspawd']
        )
        self.conversations = Conversations()
        self.config = config
        self.lock = lock
        self.handler = Message(self.config, self.conversations)

    async def event_ready(self):
        logger.info(f"Bot is ready! Logged in as Lucy_")
        logger.info(f"Connected to channel: spawdspawd")

    async def event_message(self, message):
        logger.info(f"Received message: {message.content}")
        array = await self.handler.process_array(message.content)

        # Chat
        if self.config['openai_chat_completion']:
           async for chat_completion in self.handler.generate_moderation_completion(custom_id=message.author.id, array=array):
               await message.reply(response)

        # Moderate Text and Images
        if self.config['openai_chat_moderation']:
            async for moderation_completion in self.handler.generate_moderation_completion(custom_id=message.author.id, array=array):
                full_response = json.loads(moderation_completion)
                results = full_response.get('results', [])
                flagged = results[0].get('flagged', False)
                carnism_flagged = results[0]['categories'].get('carnism', False)
                carnism_score = results[0]['category_scores'].get('carnism', 0)
                total_carnism_score = sum(arg['category_scores'].get('carnism', 0) for arg in results)
                if carnism_flagged or flagged:
                    if not self.config['discord_role_pass']:
                        await message.delete()
                    NLPUtils.append_to_other_jsonl('training.jsonl', carnism_score, message.content, message.author.id)
