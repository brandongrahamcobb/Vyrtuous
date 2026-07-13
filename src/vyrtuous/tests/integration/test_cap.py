"""!/bin/python3
test_cap.py The purpose of this program is to be the integration test for the cap command for Vyrtuous.

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

import os

import pytest

from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)

VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, channel, category, limit",
    [
        ("Administrator", "!cap", "{channel_snowflake}", "ban", "8h"),
    ],
)
async def test_cap(bot, command: str, channel, category, limit, permission_role):
    """
    Set a expires in limit for a channel by populating the PostgresSQL database
    'vyrtuous' in the table 'active_caps'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with cap
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !cap 10000000000000010 ban 8
    [{emoji} Ban cap created\n Guild1\n Channel1]
    """
    c = channel.format(
        channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
    )
    full = f"{command} {c} {category} {limit}"
    if os.environ["TEST_MODE"].lower() == "text" or os.environ["TEST_MODE"].lower() == "all":
        captured = await send_message(bot=bot, content=full)
        assert captured == ["success"]
    elif os.environ["TEST_MODE"].lower() == "app" or os.environ["TEST_MODE"].lower() == "all":
        objects = setup(bot)
        msg = build_message(
            author=objects.get("author", None),
            channel=objects.get("text_channel", None),
            content=full,
            guild=objects.get("guild", None),
            state=objects.get("state", None),
        )
        inx = interaction(
            bot=bot,
            channel=objects.get("text_channel", None),
            guild=objects.get("guild", None),
            message=msg,
        )
        async with capture_command() as end_results:
            cog = bot.get_cog("HiddenAdministratorAppCommands")
            command = cog.cap_app_command
            await command.callback(
                cog, interaction=inx, channel=c, category=category, limit=limit
            )
        for kind, content in end_results:
            assert kind == "success"
