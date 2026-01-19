import discord

from vyrtuous.tests.integration.mock_discord_guild import MockGuild
from vyrtuous.tests.integration.mock_discord_state import MockState

MEMBER_DATA =  {
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
    "avatar": None
}

class MockMember(discord.Member):

    def __init__(self, guild: MockGuild, id: int, is_bot: bool, name: str, state: MockState, **overrides):
        data = MEMBER_DATA.copy()
        data['user']['bot'] = is_bot
        data['user']['id'] = id
        data['user']['username'] = name
        data.update(overrides)
        super().__init__(data=data, guild=guild, state=state)

    async def edit(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        return self