"""!/bin/python3
test_vr.py The purpose of this program is to be the integration test for the vr command for Vyrtuous.

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
        ("Administrator", "v", "{channel_snowflake}"),
        ("Administrator", "v", "{channel_snowflake}"),
    ],
)
async def test_v(bot, command: str, prefix: str, channel, permission_role):
    """
    Create or teardown a video channel by accessing
    the PostgresSQL database 'vyrtuous' in the table 'video_channels'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with vr
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !vr <@10000000000000010>
    [{emoji} Video Room has been created]

    >>> !vr 10000000000000010
    [{emoji} Video Rooms has been deleted]
    """
    c = None
    # d = None
    full = f"{prefix}{command}"
    if channel:  # and duration is None:
        c = channel.format(channel_snowflake=VOICE_CHANNEL_SNOWFLAKE)
        full = f"{prefix}{command} {c}"
    # elif channel: and duration:
    #     c = channel.format(channel_snowflake=VOICE_CHANNEL_SNOWFLAKE)
    #     d = duration
    #     full = f"{command} {c} {d}
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
            cog = bot.get_cog("ChannelManagementAppCommands")
            command = cog.toggle_video_channel_app_command
            transformer = AppTarget()
            if c:
                resolved_channel = await transformer.transform(inx, c)
            else:
                resolved_channel = None
            await command.callback(cog, interaction=inx, channel=resolved_channel)
        for kind, content in end_results:
            assert kind == "success"
