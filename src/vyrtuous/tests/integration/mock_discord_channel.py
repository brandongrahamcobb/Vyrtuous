"""mock_discord_channel.py The purpose of this program is to support integration testing for Vyrtuous.

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

from types import SimpleNamespace

import discord

from vyrtuous.tests.integration.mock_discord_guild import MockGuild
from vyrtuous.tests.integration.mock_discord_message import MockMessage
from vyrtuous.tests.integration.mock_discord_state import MockState

GUILD_ID = 10000000000000500
TEXT_CHANNEL_ID = 10000000000000010
TEXT_CHANNEL_NAME = "text-channel"
CHANNEL_DATA = {
    "id": TEXT_CHANNEL_ID,
    "name": TEXT_CHANNEL_NAME,
    "type": discord.TextChannel,
    "position": 0,
    "nsfw": False,
    "topic": "",
    "rate_limit_per_user": 0,
    "parent_id": None,
    "guild_id": GUILD_ID,
}


class MockChannel(discord.TextChannel):

    def __init__(self, guild: MockGuild, state: MockState, **overrides):
        data = CHANNEL_DATA.copy()
        data.update(overrides)
        super().__init__(data=data, guild=guild, state=state)
        self._messages = []
        self._state = state

    def append_message(self, msg):
        self._messages.append(msg)

    async def fetch_message(self, message_snowflake):
        if message_snowflake is None:
            return None
        for msg in self._messages:
            if msg.id == message_snowflake:
                return msg
        raise ValueError(f"Message with snowflake {message_snowflake} not found")

    def permissions_for(self, member):
        return SimpleNamespace(send_messages=True)

    async def send(
        self, author=None, bot=None, content=None, embed=None, file=None, **kwargs
    ):
        msg = MockMessage(
            author=self.guild.me,
            channel=self,
            content=content,
            embed=embed,
            file=file,
            guild=self.guild,
            state=self._state,
        )
        self._messages.append(msg)
        return msg

    async def set_permissions(self, target, **overwrites):
        self.overwrites[target.id] = overwrites
        return True
