import discord

from vyrtuous.tests.integration.mock_discord_member import MockMember
from vyrtuous.tests.integration.mock_discord_state import MockState

GUILD_ID = 10000000000000500
PRIVILEGED_AUTHOR_ID = 10000000000000001
PRIVILEGED_AUTHOR_NAME = "Privileged Author Name"
MESSAGE_ID = 10000000000000100
TEXT_CHANNEL_ID = 10000000000000010
MESSAGE_DATA = {
    "id": MESSAGE_ID,
    "channel_id": TEXT_CHANNEL_ID,
    "guild_id": GUILD_ID,
    "author": {
        "id": PRIVILEGED_AUTHOR_ID,
        "username": PRIVILEGED_AUTHOR_NAME,
        "discriminator": "1234",
        "avatar": None,
        "bot": False,
        "system": False,
        "public_flags": 0,
    },
    "member": {
        "nick": "Tester",
        "roles": [],
        "joined_at": "2025-01-01T00:00:00.000000+00:00",
        "premium_since": None,
        "deaf": False,
        "mute": False,
        "pending": False,
        "permissions": 104189505,
        "communication_disabled_until": None,
    },
    "content": "",
    "timestamp": "2025-01-04T12:00:00.000000+00:00",
    "edited_timestamp": None,
    "tts": False,
    "mention_everyone": False,
    "mentions": [],
    "mention_roles": [],
    "attachments": [],
    "embeds": [],
    "reactions": [],
    "pinned": False,
    "webhook_id": None,
    "type": 0,
    "activity": None,
    "application": {
        "id": 0,
        "flags": 0,
        "cover_image": None,
        "description": "Test application",
        "icon": None,
        "name": "TestApp",
        "team": None,
        "verify_key": None,
    },
    "message_reference": {
        "message_id": 123,
        "channel_id": TEXT_CHANNEL_ID,
        "guild_id": GUILD_ID,
        "type": 0,
    },
    "flags": 0,
}

class MockMessage(discord.Message):

    def __init__(self, author: MockMember, channel, state: MockState, **overrides):
        data = MESSAGE_DATA.copy()
        data['author']['id'] = author.id
        data.update(overrides)
        super().__init__(channel=channel, data=data, state=state)

    async def add_reaction(self, emoji):
        self.reactions.append(emoji)

    async def clear_reactions(self):
        self.reactions.clear()

    async def edit(self, *, content=None, embed=None, embeds=None, view=None, **kwargs):
        if content is not None:
            self.content = content
        if embed is not None:
            self.edited_embeds.append(embed)
        if embeds is not None:
            self.embeds = embeds
        if view is not None:
            self.view = view
        return self

    async def remove_reaction(self, emoji, user):
        if emoji in self.reactions:
            self.reactions.remove(emoji)