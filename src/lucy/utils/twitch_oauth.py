import aiohttp
from quart import Quart, request, redirect
from datetime import datetime, timedelta
import asyncio
import logging

from .setup_logging import logger

twitch_app = Quart(__name__)

class TwitchOAuth:
    def __init__(self, config):
        self.config = config
        self.access_token = None
        self.refresh_token = None
        self.expires_at = None
        self.token_event = asyncio.Event()
        self.client_id = self.config['api_keys']['Twitch']['client_id']
        self.client_secret = self.config['api_keys']['Twitch']['client_secret']
        self.redirect_uri = self.config['api_keys']['Twitch']['redirect_uri']

    def get_authorization_url(self):
        return (
            f"https://id.twitch.tv/oauth2/authorize"
            f"?client_id={self.client_id}"
            f"&redirect_uri={self.redirect_uri}"
            f"&response_type=code"
            f"&scope={' '.join(self.config.get('scopes', ['chat:read', 'chat:edit']))}"
        )

    async def exchange_token(self, code):
        data = {
            "client_id": self.client_id,
            "client_secret": self.client_secret,
            "code": code,
            "grant_type": "authorization_code",
            "redirect_uri": self.redirect_uri,
        }

        async with aiohttp.ClientSession() as session:
            async with session.post("https://id.twitch.tv/oauth2/token", data=data) as resp:
                token_data = await resp.json()

        if "access_token" not in token_data:
            logger.error(f"Twitch token exchange failed: {token_data}")
            return False

        self.access_token = token_data["access_token"]
        self.refresh_token = token_data.get("refresh_token")
        self.expires_at = datetime.utcnow() + timedelta(seconds=token_data["expires_in"])

        logger.info("Successfully exchanged Twitch authorization code for an access token.")
        self.token_event.set()
        return True

    async def refresh_token_func(self):
        if not self.refresh_token:
            logger.error("No Twitch refresh token available.")
            return False

        data = {
            "client_id": self.client_id,
            "client_secret": self.client_secret,
            "refresh_token": self.refresh_token,
            "grant_type": "refresh_token",
        }

        async with aiohttp.ClientSession() as session:
            async with session.post("https://id.twitch.tv/oauth2/token", data=data) as resp:
                token_data = await resp.json()

        if "access_token" not in token_data:
            logger.error(f"Twitch token refresh failed: {token_data}")
            return False

        self.access_token = token_data["access_token"]
        self.refresh_token = token_data.get("refresh_token")
        self.expires_at = datetime.utcnow() + timedelta(seconds=token_data["expires_in"])

        logger.info("Twitch token refreshed successfully.")
        return True

    async def ensure_token(self):
        if self.access_token and datetime.utcnow() < self.expires_at:
            logger.debug("Valid Twitch token available.")
            return self.access_token

        if self.refresh_token:
            logger.debug("Refreshing Twitch token...")
            refreshed = await self.refresh_token_func()
            if refreshed:
                return self.access_token

        logger.error("No valid Twitch token available; reauthorization required.")
        return None

    async def wait_for_token(self):
        await self.token_event.wait()

def setup_twitch_routes(app, twitch_oauth):
    @app.route("/twitch_authorize")
    async def twitch_authorize():
        auth_url = twitch_oauth.get_authorization_url()
        logger.debug(f"Redirecting to Twitch OAuth URL: {auth_url}")
        return redirect(auth_url)

    @app.route("/twitch_callback")
    async def twitch_callback():
        code = request.args.get("code")
        if not code:
            logger.error("Missing Twitch authorization code in callback.")
            return "Missing authorization code", 400

        logger.debug(f"Twitch authorization code received: {code}")
        success = await twitch_oauth.exchange_token(code)
        if not success:
            return "Twitch token exchange failed.", 400

        return "Twitch authentication successful! You can close this window."

    @app.route("/twitch_validate_token")
    async def twitch_validate_token():
        token = await twitch_oauth.ensure_token()
        if not token:
            return "No valid Twitch token available. Reauthorization required.", 401

        return f"Twitch Access Token: {token}"
