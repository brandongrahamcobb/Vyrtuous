''' linkedin_oauth.py  The purpose of this program is to host the Quart app for LinkedIn OAuth 2.0.
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
from quart import Quart, request, redirect

import aiohttp
import asyncio
import logging

linkedin_app = Quart(__name__)

TOKEN_URL = 'https://www.linkedin.com/oauth/v2/accessToken'
AUTH_URL_BASE = 'https://www.linkedin.com/oauth/v2/authorization'
SCOPES = ['profile', 'w_member_social']

class LinkedInOAuth:
    def __init__(self, config):
        self.config = config
        self.access_token = None
        self.refresh_token = None
        self.expires_at = None
        self.token_event = asyncio.Event()
        self.client_id = self.config['api_keys']['LinkedIn']['client_id']
        self.client_secret = self.config['api_keys']['LinkedIn']['client_secret']
        self.redirect_uri = self.config['api_keys']['LinkedIn']['redirect_uri']

    def get_authorization_url(self):
        scopes = [
            "openid",  # OpenID scope
            "profile",  # Use your name and photo
            "r_ads_reporting",  # Retrieve reporting for your advertising accounts
            "r_organization_social",  # Retrieve your organization's posts, comments, reactions, and other engagement data
            "rw_organization_admin",  # Manage your organization's pages and retrieve reporting data
            "w_member_social",  # Create, modify, and delete posts, comments, and reactions on your behalf
            "r_ads",  # Retrieve your advertising accounts
            "w_organization_social",  # Create, modify, and delete posts, comments, and reactions on your organization's behalf
            "rw_ads",  # Manage your advertising accounts
            "r_basicprofile",  # Use your basic profile including your name, photo, headline, and public profile URL
            "r_organization_admin",  # Retrieve your organization's pages and their reporting data
            "email",  # Use the primary email address associated with your LinkedIn account
            "r_1st_connections_size"  # Retrieve the number of 1st-degree connections within your network
        ]
        scope_str = "%20".join(scopes)
        return (
            f"https://www.linkedin.com/oauth/v2/authorization"
            f"?response_type=code"
            f"&client_id={self.client_id}"
            f"&redirect_uri={self.redirect_uri}"
            f"&scope={scope_str}"
        )

    async def exchange_token(self, code):
        data = {
            "grant_type": "authorization_code",
            "code": code,
            "redirect_uri": self.redirect_uri,
            "client_id": self.client_id,
            "client_secret": self.client_secret
        }
        headers = {
            "Content-Type": "application/x-www-form-urlencoded"
        }
        async with aiohttp.ClientSession() as session:
            async with session.post(TOKEN_URL, data=data, headers=headers) as resp:
                token_data = await resp.json()
        if "access_token" not in token_data:
            logger.error(f"LinkedIn token exchange failed: {token_data}")
            return False
        self.access_token = token_data["access_token"]
        self.expires_at = datetime.utcnow() + timedelta(seconds=token_data["expires_in"])
        logger.info("Successfully exchanged LinkedIn authorization code for an access token.")
        self.token_event.set()
        return True

    async def refresh_token_func(self):
        logger.error("LinkedIn token refresh not supported.")
        return False

    async def ensure_token(self):
        if self.access_token and datetime.utcnow() < self.expires_at:
            logger.debug("Valid LinkedIn token available.")
            return self.access_token
        if self.refresh_token:
            logger.debug("Refreshing LinkedIn token...")
            refreshed = await self.refresh_token_func()
            if refreshed:
                return self.access_token
        logger.error("No valid LinkedIn token available; reauthorization required.")
        return None

    async def wait_for_token(self):
        await self.token_event.wait()

def setup_linkedin_routes(app, linkedin_oauth):
    @app.route("/linkedin_authorize")
    async def linkedin_authorize():
        auth_url = linkedin_oauth.get_authorization_url()
        logger.debug(f"Redirecting to LinkedIn OAuth URL: {auth_url}")
        return redirect(auth_url)

    @app.route("/linkedin_callback")
    async def linkedin_callback():
        code = request.args.get("code")
        if not code:
            logger.error("Missing LinkedIn authorization code in callback.")
            return "Missing authorization code", 400
        logger.debug(f"LinkedIn authorization code received: {code}")
        success = await linkedin_oauth.exchange_token(code)
        if not success:
            return "LinkedIn token exchange failed.", 400
        return "LinkedIn authentication successful! You can close this window."

    @app.route("/linkedin_validate_token")
    async def linkedin_validate_token():
        token = await linkedin_oauth.ensure_token()
        if not token:
            return "No valid LinkedIn token available. Reauthorization required.", 401
        return f"LinkedIn Access Token: {token}"
