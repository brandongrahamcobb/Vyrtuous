"""!/bin/python3
test_survey.py The purpose of this program is to be the integration test for the survey command for Vyrtuous.

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

COMMAND = "survey"
BASE_PERMISSIONS = ["command.info.survey", "command.info.scope.channel"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "channel",
    [
        (None),
        ("{channel_snowflake}"),
        ("<#{channel_snowflake}>"),
    ],
)
async def test_survey_text_command(
    bot,
    prefix: str,
    channel: str | None,
):
    docstring = """
    List member-group association from the PermissionState cache.

    Parameters
    ----------
    channel (Optional) : str | int
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> !survey
    Embed
    """
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
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions_at_all",
                    side_effect=check_permissions(BASE_PERMISSIONS),
                )
            )

            full = f"{prefix}{COMMAND}"
            if channel is None:
                c = channel
            else:
                c = channel.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                )
                full += f" {c}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "channel",
    [
        (None),
        ("{channel_snowflake}"),
        ("<#{channel_snowflake}>"),
    ],
)
async def test_survey_app_command(
    bot,
    channel: str | None,
):
    docstring = """
    List member-group association from the PermissionState cache.

    Parameters
    ----------
    channel (Optional) : str | int
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> /survey
    Embed
    """
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
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions_at_all",
                    side_effect=check_permissions(BASE_PERMISSIONS),
                )
            )

            cog = bot.get_cog("InfoAppCommands")
            command = cog.survey_app_command
            if channel is None:
                c = channel
            else:
                c = channel.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                )
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
                transformer = AppTarget()
                if c:
                    resolved_channel = await transformer.transform(inx, c)
                else:
                    resolved_channel = None
                await command.callback(
                    cog,
                    interaction=inx,
                    channel=resolved_channel,
                )
            for kind, content in end_results:
                assert kind == "success"
