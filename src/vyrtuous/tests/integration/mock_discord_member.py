"""!/bin/python3
mock_discord_member.py The purpose of this program is to support integration testing for Vyrtuous.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import discord

from vyrtuous.tests.integration.mock_discord_guild import MockGuild
from vyrtuous.tests.integration.mock_discord_role import ROLE_SNOWFLAKE
from vyrtuous.tests.integration.mock_discord_state import MockState

MEMBER_DATA = {
    "user": {
        "id": "",
        "username": "",
        "discriminator": "1234",
        "avatar": None,
        "bot": False,
    },
    "nick": "",
    "roles": [],
    "joined_at": "2025-01-01T00:00:00.000000+00:00",
    "premium_since": None,
    "pending": False,
    "deaf": False,
    "mute": False,
    "communication_disabled_until": None,
    "flags": 0,
    "avatar": None,
}


class MockMember(discord.Member):
    def __init__(
        self,
        guild: MockGuild,
        id: int,
        is_bot: bool,
        name: str,
        state: MockState,
        **overrides,
    ):
        self.data = MEMBER_DATA.copy()
        self.data["user"]["bot"] = is_bot
        self.data["user"]["id"] = id
        self.data["user"]["username"] = name
        self.data.update(overrides)
        super().__init__(data=self.data, guild=guild, state=state)
        self._mock_roles = []

    async def edit(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        return self

    @property
    def roles(self):
        return self._mock_roles

    async def add_roles(self, *roles, reason=None, **kwargs):
        for role in roles:
            if role not in self._mock_roles:
                self._mock_roles.append(role)
