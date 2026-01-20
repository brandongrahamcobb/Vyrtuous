"""mock_discord_role.py The purpose of this program is to support integration testing for Vyrtuous.

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

from vyrtuous.tests.integration.mock_discord_guild import MockGuild
from vyrtuous.tests.integration.mock_discord_state import MockState

GUILD_SNOWFLAKE = 10000000000000500
ROLE_SNOWFLAKE = 10000000000000200
ROLE_NAME = "Role Name"
ROLE_DATA = {
    "id": ROLE_SNOWFLAKE,
    "name": ROLE_NAME,
    "color": 0,
    "hoist": False,
    "position": 0,
    "permissions": 0,
    "managed": False,
    "mentionable": False,
    "guild_id": GUILD_SNOWFLAKE,
}


class MockRole(discord.Role):

    def __init__(self, guild: MockGuild, state: MockState, **overrides):
        data = ROLE_DATA.copy()
        data.update(overrides)
        super().__init__(data=data, guild=guild, state=state)
