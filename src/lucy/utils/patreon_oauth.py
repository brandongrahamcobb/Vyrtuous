''' patreon_oauth.py  The purpose of this program is to host the Quart app for Twitch OAuth 2.0.
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
from datetime import datetime, timedelta
from lucy.utils.setup_logging import logger
from quart import Quart, request, redirect, jsonify

import aiohttp
import asyncio
import logging

patreon_app = Quart(__name__)

class PatreonOAuth:
    def __init__(self, config):
        self.config = config
        self.access_token = None
        self.refresh_token = None
        self.expires_at = None
        self.token_event = asyncio.Event()
        self.client_id = self.config['api_keys']['Patreon']['client_id']
        self.client_secret = self.config['api_keys']['Patreon']['client_secret']
        self.redirect_uri = self.config['api_keys']['Patreon']['redirect_uri']

    def get_authorization_url(self):
        return (
            f"https://www.patreon.com/oauth2/authorize"
            f"?response_type=code"
            f"&client_id={self.client_id}"
            f"&redirect_uri={self.redirect_uri}"
            f"&scope={'%20'.join(self.config.get('scopes', ['identity', 'campaigns']))}"
        )

    async def exchange_token(self, code):
        data = {
            "client_id": self.client_id,
            "client_secret": self.client_secret,
            "code": code,
            "grant_type": "authorization_code",
            "redirect_uri": self.redirect_uri,
        }
        headers = {
            "Content-Type": "application/x-www-form-urlencoded"
        }
        async with aiohttp.ClientSession() as session:
            async with session.post("https://www.patreon.com/api/oauth2/token", data=data, headers=headers) as resp:
                if resp.status != 200:
                    error_text = await resp.text()
                    logger.error(f"Patreon token exchange failed with status {resp.status}: {error_text}")
                    return False
                token_data = await resp.json()
        if "access_token" not in token_data:
            logger.error(f"Patreon token exchange failed: {token_data}")
            return False
        self.access_token = token_data["access_token"]
        self.refresh_token = token_data.get("refresh_token")
        self.expires_at = datetime.utcnow() + timedelta(seconds=token_data.get("expires_in", 3600))
        logger.info("Successfully exchanged Patreon authorization code for an access token.")
        self.token_event.set()
        return True

    async def refresh_token_func(self):
        if not self.refresh_token:
            logger.error("No Patreon refresh token available.")
            return False
        data = {
            "client_id": self.client_id,
            "client_secret": self.client_secret,
            "refresh_token": self.refresh_token,
            "grant_type": "refresh_token",
        }
        headers = {
            "Content-Type": "application/x-www-form-urlencoded"
        }
        async with aiohttp.ClientSession() as session:
            async with session.post("https://www.patreon.com/api/oauth2/token", data=data, headers=headers) as resp:
                if resp.status != 200:
                    error_text = await resp.text()
                    logger.error(f"Patreon token refresh failed with status {resp.status}: {error_text}")
                    return False
                token_data = await resp.json()
        if "access_token" not in token_data:
            logger.error(f"Patreon token refresh failed: {token_data}")
            return False
        self.access_token = token_data["access_token"]
        self.refresh_token = token_data.get("refresh_token")
        self.expires_at = datetime.utcnow() + timedelta(seconds=token_data.get("expires_in", 3600))
        logger.info("Patreon token refreshed successfully.")
        return True

    async def ensure_token(self):
        if self.access_token and datetime.utcnow() < self.expires_at:
            logger.debug("Valid Patreon token available.")
            return self.access_token
        if self.refresh_token:
            logger.debug("Refreshing Patreon token...")
            refreshed = await self.refresh_token_func()
            if refreshed:
                return self.access_token
        logger.error("No valid Patreon token available; reauthorization required.")
        return None

    async def wait_for_token(self):
        await self.token_event.wait()

def setup_patreon_routes(app: Quart, patreon_oauth: PatreonOAuth):
    @app.route("/patreon_authorize")
    async def patreon_authorize():
        auth_url = patreon_oauth.get_authorization_url()
        logger.debug(f"Redirecting to Patreon OAuth URL: {auth_url}")
        return redirect(auth_url)

    @app.route("/patreon_callback")
    async def patreon_callback():
        code = request.args.get("code")
        if not code:
            logger.error("Missing authorization code in callback.")
            return "Missing authorization code", 400
        logger.debug(f"Patreon authorization code received: {code}")
        success = await patreon_oauth.exchange_token(code)
        if not success:
            return "Patreon token exchange failed.", 400
        return "Patreon authentication successful! You can close this window.", 200

    @app.route("/patreon_validate_token")
    async def patreon_validate_token():
        token = await patreon_oauth.ensure_token()
        if not token:
            return "No valid Patreon token available. Reauthorization required.", 401

        return jsonify({"Patreon Access Token": token})
