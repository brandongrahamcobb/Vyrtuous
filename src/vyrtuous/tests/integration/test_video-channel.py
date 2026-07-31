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
from contextlib import ExitStack
from unittest.mock import patch

import pytest

from vyrtuous.db.video_channel import VideoChannel
from vyrtuous.models.duration import AppDuration
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    send_message,
    setup,
)

VOICE_CHANNEL_SNOWFLAKE = 10000000000000011

COMMAND = "video-only"
BASE_PERMISSIONS = ["command.channel.video-channel"]
TABLE_NAME = VideoChannel.__tablename__


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "channel, duration",
    [
        ("{channel_snowflake}", None),
        ("<#{channel_snowflake}>", "8h"),
    ],
)
async def test_video_channel_text_command(
    bot, prefix: str, channel: str | None, duration: str | None
):
    docstring = """
    Toggle video-only for a channel in the PostgreSQL database table "active_video_only_channels".

    Parameters
    ----------
    channel : str | int
        Resolves to: discord.VoiceChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>
    duration : str
        Resolves to: DurationObject
        Examples: 0 | 1 | 30m | 1h | 1d | 1w 

    Example
    --------
    >>> !video-channel
    Embed
    """
    assert TABLE_NAME in docstring
    assert COMMAND in docstring
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions",
                    side_effect=check_permissions(BASE_PERMISSIONS),
                )
            )
            full = f"{prefix}{COMMAND}"
            if channel is None:
                c = channel
            else:
                c = channel.format(channel_snowflake=VOICE_CHANNEL_SNOWFLAKE)
                full += f" {c}"
            if duration is None:
                d = duration
            else:
                d = duration
                full += f" {d}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "channel, duration",
    [
        ("{channel_snowflake}", None),
        ("<#{channel_snowflake}>", "8h"),
    ],
)
async def test_video_channel_app_command(
    bot, channel: str | None, duration: str | None
):
    docstring = """
    Toggle video-only for a channel in the PostgreSQL database table "active_video_only_channels".

    Parameters
    ----------
    channel : str | int
        Resolves to: discord.VoiceChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>
    duration : str
        Resolves to: DurationObject
        Examples: 0 | 1 | 30m | 1h | 1d | 1w 

    Example
    --------
    >>> /video-channel
    Embed
    """
    assert TABLE_NAME in docstring
    assert COMMAND in docstring
    if (
        os.environ["TEST_MODE"].lower() == "app"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions",
                    side_effect=check_permissions(BASE_PERMISSIONS),
                )
            )
            cog = bot.get_cog("ChannelManagementAppCommands")
            command = cog.toggle_video_only_channel_app_command
            if channel is None:
                c = channel
            else:
                c = channel.format(channel_snowflake=VOICE_CHANNEL_SNOWFLAKE)
            if duration is None:
                d = duration
            else:
                d = duration
            objects = setup(bot)
            msg = build_message(
                author=objects.get("author", None),
                channel=objects.get("voice_channel", None),
                content="",
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
                target_transformer = AppTarget()
                if c:
                    resolved_channel = await target_transformer.transform(inx, c)
                else:
                    resolved_channel = c
                duration_transformer = AppDuration()
                if d:
                    resolved_duration = await duration_transformer.transform(inx, d)
                else:
                    resolved_duration = d
                await command.callback(
                    cog,
                    interaction=inx,
                    channel=resolved_channel,
                    duration=resolved_duration,
                )
            for kind, content in end_results:
                assert kind == "success"
