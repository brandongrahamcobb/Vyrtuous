"""!/bin/python3
test_stage.py The purpose of this program is to be the integration test for the stage command for Vyrtuous.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import os

import pytest

from vyrtuous.models.duration import AppDuration
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
    "permission_role, command, channel, duration",
    [
        ("Coordinator", "automute", None, None),
        ("Coordinator", "automute", "{channel_snowflake}", None),
        ("Coordinator", "automute", "{channel_snowflake}", "1h"),
    ],
)
async def test_automute(
    bot, command: str, prefix: str, channel, duration, permission_role
):
    """
    Create or teardown a stage by accessing
    the PostgresSQL database 'vyrtuous' in the table 'video_channels'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with stage
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !stage <@10000000000000010>
    [{emoji} Stage Room has been created]

    >>> !stage 10000000000000010
    [{emoji} Stage has been ended]
    """
    c = None
    d = None
    full = f"{prefix}{command}"
    if channel and duration is None:
        c = channel.format(channel_snowflake=VOICE_CHANNEL_SNOWFLAKE)
        full = f"{prefix}{command} {c}"
    elif channel and duration:
        c = channel.format(channel_snowflake=VOICE_CHANNEL_SNOWFLAKE)
        d = duration
        full = f"{prefix}{command} {c} {d}"
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
            channel=objects.get("voice_channel", None),
            guild=objects.get("guild", None),
            message=msg,
        )
        async with capture_command() as end_results:
            cog = bot.get_cog("CoordinatorAppCommands")
            command = cog.toggle_automute_app_command
            channel_transformer = AppTarget()
            if c:
                resolved_channel = await channel_transformer.transform(inx, c)
            else:
                resolved_channel = None
            duration_transformer = AppDuration()
            if d:
                resolved_duration = await duration_transformer.transform(inx, d)
            else:
                resolved_duration = None
            await command.callback(
                cog,
                interaction=inx,
                channel=resolved_channel,
                duration=resolved_duration,
            )
        for kind, content in end_results:
            assert kind == "success"
