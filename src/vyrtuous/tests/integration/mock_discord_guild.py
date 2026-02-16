"""!/bin/python3
mock_discord_guild.py The purpose of this program is to support integration testing for Vyrtuous.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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
"""

import discord

from vyrtuous.tests.integration.mock_discord_bot import MockBot
from vyrtuous.tests.integration.mock_discord_state import MockState

GUILD_SNOWFLAKE = 10000000000000500
GUILD_NAME = "Guild Name"
GUILD_DATA = {
    "id": GUILD_SNOWFLAKE,
    "name": GUILD_NAME,
    "icon": None,
    "splash": None,
    "discovery_splash": None,
    "owner_id": 154749533429956608,
    "region": None,
    "afk_channel_id": None,
    "afk_timeout": 300,
    "widget_enabled": False,
    "widget_channel_id": None,
    "verification_level": 0,
    "default_message_notifications": 0,
    "explicit_content_filter": 0,
    "roles": [],
    "emojis": [],
    "features": [],
    "mfa_level": 0,
    "application_id": None,
    "system_channel_id": None,
    "system_channel_flags": 0,
    "rules_channel_id": None,
    "joined_at": "2025-01-01T00:00:00.000000+00:00",
    "large": False,
    "member_count": 0,
    "voice_states": {},
    "members": {},
    "channels": [],
    "presences": {},
    "max_presences": None,
    "max_members": None,
    "vanity_url_code": None,
    "description": None,
    "banner": None,
    "premium_tier": 0,
    "premium_subscription_count": 0,
    "preferred_locale": "en-US",
    "public_updates_channel_id": None,
    "max_video_channel_users": 25,
    "approximate_member_count": None,
    "approximate_presence_count": None,
    "nsfw_level": 0,
}


class MockGuild(discord.Guild):
    def __init__(
        self, bot: MockBot, channels, members, roles, state: MockState, **overrides
    ):
        data = GUILD_DATA.copy()
        data.update(overrides)
        super().__init__(data=data, state=state)
        self._channels = channels
        self._members = members
        self._voice_channels = {}
        self._roles = roles

    def get_channel(self, channel_snowflake):
        if channel_snowflake is None:
            return None
        for channel in self._channels:
            if str(channel.id) == str(channel_snowflake):
                return channel

    def get_member(self, member_snowflake):
        if member_snowflake is None:
            return None
        if member_snowflake in self._members.keys():
            return self._members[member_snowflake]
        else:
            return None
            # if str(member.id) == str(member_snowflake):
            #     return member

    def get_role(self, role_snowflake):
        if role_snowflake is None:
            return None
        if role_snowflake in self._roles.keys():
            return self._roles[role_snowflake]
        else:
            return None

    @property
    def channels(self):
        return self._channels

    @property
    def voice_channels(self):
        return self._voice_channels

    @property
    def members(self):
        return self._members.values()

    @property
    def roles(self):
        return self._roles
