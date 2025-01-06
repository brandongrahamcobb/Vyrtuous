import aiohttp
from datetime import datetime, timedelta
import logging

logger = logging.getLogger(__name__)

LINKEDIN_CLIENT_ID = 'your_linkedin_client_id'
LINKEDIN_CLIENT_SECRET = 'your_linkedin_client_secret'
LINKEDIN_REDIRECT_URI = 'http://localhost:5000/callback'

TWITCH_CLIENT_ID = 'your_twitch_client_id'
TWITCH_CLIENT_SECRET = 'your_twitch_client_secret'
TWITCH_REDIRECT_URI = 'http://localhost:5000/callback'

def linkedin_authorization_url():
    return (
        f"https://www.linkedin.com/oauth/v2/authorization?response_type=code"
        f"&client_id={LINKEDIN_CLIENT_ID}&redirect_uri={LINKEDIN_REDIRECT_URI}"
        f"&scope=r_liteprofile%20r_emailaddress%20w_member_social"
    )

async def linkedin_exchange_token(code):
    url = "https://www.linkedin.com/oauth/v2/accessToken"
    data = {
        "grant_type": "authorization_code",
        "code": code,
        "redirect_uri": LINKEDIN_REDIRECT_URI,
        "client_id": LINKEDIN_CLIENT_ID,
        "client_secret": LINKEDIN_CLIENT_SECRET
    }
    async with aiohttp.ClientSession() as session:
        async with session.post(url, data=data) as response:
            token_data = await response.json()
            if 'access_token' in token_data:
                return token_data['access_token']
            logger.error(f"Error exchanging code for LinkedIn token: {token_data}")
            return None

# Twitch OAuth
def twitch_authorization_url():
    return (
        f"https://id.twitch.tv/oauth2/authorize?client_id={TWITCH_CLIENT_ID}"
        f"&redirect_uri={TWITCH_REDIRECT_URI}&response_type=code&scope=chat:read chat:edit"
    )

async def twitch_exchange_token(code):
    url = "https://id.twitch.tv/oauth2/token"
    data = {
        "client_id": TWITCH_CLIENT_ID,
        "client_secret": TWITCH_CLIENT_SECRET,
        "code": code,
        "grant_type": "authorization_code",
        "redirect_uri": TWITCH_REDIRECT_URI,
    }
    async with aiohttp.ClientSession() as session:
        async with session.post(url, data=data) as response:
            token_data = await response.json()
            if 'access_token' in token_data:
                return token_data['access_token']
            logger.error(f"Error exchanging code for Twitch token: {token_data}")
            return None

async def refresh_token(platform, refresh_token):
    if platform == 'linkedin':
        url = "https://www.linkedin.com/oauth/v2/accessToken"
        data = {
            "grant_type": "refresh_token",
            "refresh_token": refresh_token,
            "client_id": LINKEDIN_CLIENT_ID,
            "client_secret": LINKEDIN_CLIENT_SECRET
        }
    elif platform == 'twitch':
        url = "https://id.twitch.tv/oauth2/token"
        data = {
            "client_id": TWITCH_CLIENT_ID,
            "client_secret": TWITCH_CLIENT_SECRET,
            "refresh_token": refresh_token,
            "grant_type": "refresh_token"
        }
    else:
        return None

    async with aiohttp.ClientSession() as session:
        async with session.post(url, data=data) as response:
            token_data = await response.json()
            if 'access_token' in token_data:
                return token_data['access_token']
            logger.error(f"Error refreshing {platform} token: {token_data}")
            return None
