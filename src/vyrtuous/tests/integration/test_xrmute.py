"""!/bin/python3
test_rmute_xrmute.py The purpose of this program is to be the integration test for the rmute and xrmute commands for Vyrtuous.

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

from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)

ROLE_SNOWFLAKE = 10000000000000200
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, channel",
    [
        ("Administrator", "xrmute", "{channel_snowflake}"),
        ("Administrator", "xrmute", "<#{channel_snowflake}>"),
    ],
)
async def test_rmute_xrmute(bot, command: str, prefix: str, channel, permission_role):
    """
    Voice-mute a whole channel and undo it by adding and removing
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
        channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
    )
    full = f"{prefix}{command} {c}"
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        captured = await send_message(bot=bot, content=full)
        assert captured == ["success"]
    if (
        os.environ["TEST_MODE"].lower() == "app"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
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
            cog = bot.get_cog("AdministratorAppCommands")
            command = cog.channel_mute_app_command
            transformer = AppTarget()
            if c:
                resolved = await transformer.transform(inx, c)
            else:
                resolved = None
            command = cog.channel_unmute_app_command
            await command.callback(cog, interaction=inx, channel=resolved)
            for kind, content in end_results:
                assert kind == "success"
