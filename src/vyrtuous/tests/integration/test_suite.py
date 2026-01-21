"""test_suite.py The purpose of this program is to provide the shared test variables for tests using Discord objects.

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

from contextlib import asynccontextmanager
import asyncio

from vyrtuous.tests.integration.mock_discord_channel import MockChannel
from vyrtuous.tests.integration.mock_discord_guild import MockGuild
from vyrtuous.tests.integration.mock_discord_member import MockMember
from vyrtuous.tests.integration.mock_discord_message import MockMessage
from vyrtuous.tests.integration.mock_discord_role import MockRole
from vyrtuous.tests.integration.mock_discord_state import MockState

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


@asynccontextmanager
async def capture(channel):
    before = list(channel._messages)
    yield
    after = channel._messages
    new_messages = after[len(before) :]
    channel._captured = new_messages


async def send_message(bot, content: str = None):
    state = MockState()
    guild = MockGuild(bot=bot, channels={}, members=[], roles=[], state=state)
    role = MockRole(guild=guild, state=state)
    guild._roles.append(role)
    bot._guilds.append(guild)
    channel = MockChannel(bot=bot, guild=guild, state=state)
    guild._channels.update({channel.id: channel})
    author = MockMember(
        bot=bot,
        guild=guild,
        id=NOT_PRIVILEGED_AUTHOR_SNOWFLAKE_ONE,
        is_bot=False,
        name=NOT_PRIVILEGED_AUTHOR_NAME_ONE,
        state=state,
    )
    dummy = MockMember(
        bot=bot,
        guild=guild,
        id=DUMMY_MEMBER_SNOWFLAKE,
        is_bot=False,
        name=DUMMY_MEMBER_NAME,
        state=state,
    )
    bot_member = MockMember(
        bot=bot,
        guild=guild,
        id=PRIVILEGED_AUTHOR_SNOWFLAKE,
        is_bot=True,
        name=PRIVILEGED_AUTHOR_NAME,
        state=state,
    )
    channel._members.append(dummy)
    guild._members.append(author)
    guild._members.append(dummy)
    guild._members.append(bot_member)
    state.user = bot_member
    bot._connection = state
    bot.me = bot_member
    bot._state = state
    msg = MockMessage(
        author=author, channel=channel, content=content, guild=guild, state=state
    )

    async with capture(channel):
        bot.loop = asyncio.get_running_loop()
        bot.dispatch("message", msg)
        await asyncio.sleep(0.1)

    return channel._captured[-1]
