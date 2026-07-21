"""!/bin/python3
test_survey.py The purpose of this program is to be the integration test for the survey command for Vyrtuous.

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

VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, channel",
    [
        ("Moderator", "survey", None),
        ("Moderator", "survey", "{channel_snowflake}"),
        ("Moderator", "survey", "<#{channel_snowflake}>"),
    ],
)
async def test_survey(bot, command: str, prefix: str, channel, permission_role):
    """
    Server mute a member localized to the guild

    Parameters
    ----------
    None
        No paramers neccesary

    Examples
    --------

    >>> !survey <#10000000000000010>
    [{emoji} Survey results for Channel1]

    >>> !survey 10000000000000010
    [{emoji} Survey results for Channel1]
    """
    full = f"{prefix}{command}"
    c = None
    if channel:
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
            channel=objects.get("voice_channel", None),
            content=full,
            guild=objects.get("guild", None),
            state=objects.get("state", None),
        )
        inx = interaction(
            bot=bot,
            channel=objects.get("voice_channel", None),
            guild=objects.get("guild", None),
            message=msg,
        )
        async with capture_command() as end_results:
            cog = bot.get_cog("HiddenModeratorAppCommands")
            command = cog.survey_app_command
            transformer = AppTarget()
            if c:
                resolved_channel = await transformer.transform(inx, c)
            else:
                resolved_channel = None
            await command.callback(cog, interaction=inx, channel=resolved_channel)
        for kind, content in end_results:
            assert kind == "success"
