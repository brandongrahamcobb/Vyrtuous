"""test_ping.py The purpose of this program is to be the integration test for the ping command for Vyrtuous.

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


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, lines",
    [
        ("!debug", "10"),
    ],
)
async def test_ping(bot, command: str, lines):
    """
    Fetch last lines from the logs.

    Parameters
    ----------
    lines
        The number of lines

    Examples
    --------
    >>> !debug 10
    [{emoji} Logger Information]\n Scheduled blah blah blah.
    """
    full = f"{command} {lines}"
    captured = await send_message(bot=bot, content=command)
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
    admin_commands = bot.get_cog("AdminTextCommands")
    command = await admin_commands.debug_text_command(ctx, lines=int(lines))
