"""!/bin/python3
test_suite.py The purpose of this program is to provide the shared test variables for tests using Discord objects.

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

import asyncio
import threading
from contextlib import asynccontextmanager
from types import SimpleNamespace
from unittest.mock import AsyncMock

from vyrtuous.tests.integration.mock_discord_channel import (
    MockTextChannel,
    MockVoiceChannel,
)
from vyrtuous.tests.integration.mock_discord_guild import MockGuild
from vyrtuous.tests.integration.mock_discord_member import MockMember
from vyrtuous.tests.integration.mock_discord_message import MockMessage
from vyrtuous.tests.integration.mock_discord_role import MockRole
from vyrtuous.tests.integration.mock_discord_state import MockState
from vyrtuous.utils.state_service import StateService

PRIVILEGED_AUTHOR_SNOWFLAKE = 10000000000000001
PRIVILEGED_AUTHOR_NAME = "Privileged Author Name"
NOT_PRIVILEGED_AUTHOR_SNOWFLAKE_ONE = 10000000000000002
NOT_PRIVILEGED_AUTHOR_NAME_ONE = "Not Privileged Author Name One"
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
DUMMY_MEMBER_NAME = "Not Privileged Author Name Two"
RED = "\033[91m"
YELLOW = "\033[93m"
GREEN = "\033[92m"
RESET = "\033[0m"
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@asynccontextmanager
async def capture(channel):
    called = []
    before = list(channel._messages)
    channel._end_result = []
    original_end = StateService.end

    async def patched_end(self, *args, **kwargs):
        called.append((args, kwargs))
        self._add_reactions = AsyncMock()
        if "success" in kwargs:
            channel._end_result.append("success")
        elif "warning" in kwargs:
            channel._end_result.append("warning")
        elif "error" in kwargs:
            channel._end_result.append("error")
        else:
            print(f"   ⚪ No success/warning/error found in kwargs")
        return None

    StateService.end = patched_end

    try:
        yield
        await asyncio.sleep(0.5)
    finally:
        StateService.end = original_end
        after = channel._messages
        if not after:
            await asyncio.sleep(1)
        new_messages = after[len(before) :]
        channel._captured = new_messages


def build_guild(bot, state):
    guild = MockGuild(bot=bot, channels={}, members={}, roles={}, state=state)
    return guild


def build_role(guild, state):
    role = MockRole(guild=guild, state=state)
    return role


def build_voice_channel(bot, guild, state, id=None):
    channel = MockVoiceChannel(bot=bot, id=id, guild=guild, state=state)
    return channel


def build_text_channel(bot, guild, state, id=None):
    channel = MockTextChannel(bot=bot, id=id, guild=guild, state=state)
    return channel


def build_member(bot, guild, id, is_bot, name, state):
    member = MockMember(
        bot=bot,
        guild=guild,
        id=id,
        is_bot=is_bot,
        name=name,
        state=state,
    )
    return member


def build_message(author, channel, content, guild, state):
    msg = MockMessage(
        author=author, channel=channel, content=content, guild=guild, state=state
    )
    return msg


def setup(bot):
    state = MockState()
    guild = build_guild(bot, state)
    role = build_role(guild, state)
    guild._roles.update({role.id: role})
    bot._guilds.update({guild.id: guild})
    text_channel = build_text_channel(bot, guild, state, id=TEXT_CHANNEL_SNOWFLAKE)
    voice_channel = build_voice_channel(bot, guild, state, id=VOICE_CHANNEL_SNOWFLAKE)
    guild._channels.update({text_channel.id: text_channel})
    guild._channels.update({voice_channel.id: voice_channel})
    guild._voice_channels.update({voice_channel.id: voice_channel})
    author = build_member(
        bot=bot,
        guild=guild,
        id=NOT_PRIVILEGED_AUTHOR_SNOWFLAKE_ONE,
        is_bot=False,
        name=NOT_PRIVILEGED_AUTHOR_NAME_ONE,
        state=state,
    )
    dummy = build_member(
        bot=bot,
        guild=guild,
        id=DUMMY_MEMBER_SNOWFLAKE,
        is_bot=False,
        name=DUMMY_MEMBER_NAME,
        state=state,
    )
    bot_member = build_member(
        bot=bot,
        guild=guild,
        id=PRIVILEGED_AUTHOR_SNOWFLAKE,
        is_bot=True,
        name=PRIVILEGED_AUTHOR_NAME,
        state=state,
    )
    # voice_channel._members.append(dummy)
    guild._members.update({author.id: author})
    guild._members.update({dummy.id: dummy})
    guild._members.update({bot_member.id: bot_member})
    state.user = bot_member
    bot._connection = state
    bot._get_websocket = lambda *a, **k: SimpleNamespace(is_ratelimited=lambda: False)
    bot.me = bot_member
    bot._state = state
    objects = {
        "author": author,
        "bot": bot,
        "text_channel": text_channel,
        "voice_channel": voice_channel,
        "guild": guild,
        "state": state,
    }
    return objects


async def send_message(bot, content: str = None):
    objects = setup(bot)

    msg = build_message(
        author=objects.get("author", None),
        channel=objects.get("text_channel", None),
        content=content,
        guild=objects.get("guild", None),
        state=objects.get("state", None),
    )

    async with capture(objects.get("text_channel", None)):
        bot.loop = asyncio.get_running_loop()
        bot.dispatch("message", msg)
    return objects.get("text_channel", None)._end_result


@asynccontextmanager
async def capture_command():
    results = []
    original_end = StateService.end

    async def patched_end(self, *args, **kwargs):
        if getattr(self, "_ended", False):
            return
        self._ended = True
        if "success" in kwargs:
            results.append(("success", kwargs["success"]))
        elif "warning" in kwargs:
            results.append(("warning", kwargs["warning"]))
        elif "error" in kwargs:
            results.append(("error", kwargs["error"]))
        self._add_reactions = AsyncMock()

    StateService.end = patched_end
    try:
        yield results
    finally:
        StateService.end = original_end
