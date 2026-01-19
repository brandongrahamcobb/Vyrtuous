
import discord

from vyrtuous.tests.integration.mock_discord_guild import MockGuild
from vyrtuous.tests.integration.mock_discord_state import MockState

GUILD_ID = 10000000000000500
ROLE_ID = 10000000000000200
ROLE_NAME = "Role Name"
ROLE_DATA = {
    "id": ROLE_ID,
    "name": ROLE_NAME,
    "color": 0,
    "hoist": False,
    "position": 0,
    "permissions": 0,
    "managed": False,
    "mentionable": False,
    "guild_id": GUILD_ID,
}

class MockRole(discord.Role):
    
    def __init__(self, guild: MockGuild, state: MockState, **overrides):
        data = ROLE_DATA.copy()
        data.update(overrides)
        super().__init__(data=data, guild=guild, state=state)