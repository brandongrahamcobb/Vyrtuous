"""test_rmute_xrmute.py The purpose of this program is to be the integration test for the rmute and xrmute commands for Vyrtuous.

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

from typing import Optional

import pytest

from vyrtuous.tests.integration.conftest import context
from vyrtuous.tests.integration.test_suite import build_message, send_message, setup

ROLE_SNOWFLAKE = 10000000000000200
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, channel",
    [
        ("!ow", "{channel_snowflake}"),
        ("!ow", "<#{channel_snowflake}>"),
    ],
)
async def test_overwrites(bot, command: str, channel):
    """
    Voice-mute a whole room and undo it by adding and removing
    entries in the PostgreSQL database 'vyrtuous' in the table
    'active_voice_mutes'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !rmute 10000000000000010
    [{emoji} Room Muted\n Member1\n Member2]

    >>> !xrmute 10000000000000010
    [{emoji} Room Unmuted\n Member1\n Member2]

    >>> !rmute <#10000000000000010>
    [{emoji} Room Muted\n Member1\n Member2]

    >>> !xrmute <#10000000000000010>
    [{emoji} Room Unmuted\n Member1\n Member2]
    """
    c = channel.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
    )
    full = f"{command} {c}"
    captured = await send_message(bot=bot, content=full)
    assert captured
    objects = setup(bot)
    msg = build_message(
        author=objects.get("author", None),
        channel=objects.get("text_channel", None),
        content=full,
        guild=objects.get("guild", None),
        state=objects.get("state", None),
    )
    ctx = context(
        bot=bot,
        channel=objects.get("text_channel", None),
        guild=objects.get("guild", None),
        message=msg,
        prefix="!",
    )
    dev_commands = bot.get_cog("DevTextCommands")
    command = await dev_commands.list_overwrites_text_command(ctx, channel=c)
