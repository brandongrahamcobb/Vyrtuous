"""!/bin/python3
test_roleid.py The purpose of this program is to be the integration test for the roleid list command for Vyrtuous.

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
VOICE_CHANNEL_TWO_SNOWFLAKE = 10000000000000012

COMMAND = "rmv"
BASE_PERMISSIONS = ["command.utility.move"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target_channel, source_channel",
    [
        (
            "{target_channel_snowflake}",
            "{source_channel_snowflake}",
        ),
        (
            "{target_channel_snowflake}",
            None,
        ),
    ],
)
async def test_rmv_text_command(
    bot, prefix: str, target_channel: str, source_channel: str | None
):
    docstring = """
    Move members from one voice-channel to another voice-channel.

    Parameters
    ----------
    target_channel : str | int
        Resolves to: discord.VoiceChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> !rmv 10000000000000011
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
            tc = target_channel.format(target_channel_snowflake=VOICE_CHANNEL_SNOWFLAKE)
            full = f"{prefix}{COMMAND} {tc}"
            if source_channel is None:
                sc = source_channel
            else:
                sc = source_channel.format(
                    source_channel_snowflake=VOICE_CHANNEL_TWO_SNOWFLAKE,
                )
                full += f" {sc}"
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target_channel, source_channel",
    [
        (
            "{target_channel_snowflake}",
            "{source_channel_snowflake}",
        ),
        (
            "{target_channel_snowflake}",
            None,
        ),
    ],
)
async def test_rmv_app_command(bot, target_channel: str, source_channel: str | None):
    docstring = """
    Move members from one voice-channel to another voice-channel.

    Parameters
    ----------
    target_channel : str | int
        Resolves to: discord.VoiceChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> /rmv 10000000000000011
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
            cog = bot.get_cog("UtilityAppCommands")
            command = cog.channel_move_all_app_command
            tc = target_channel.format(target_channel_snowflake=VOICE_CHANNEL_SNOWFLAKE)
            if source_channel is None:
                sc = source_channel
            else:
                sc = source_channel.format(
                    source_channel_snowflake=VOICE_CHANNEL_TWO_SNOWFLAKE,
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
                resolved_target_channel = await transformer.transform(inx, tc)
                if sc:
                    resolved_source_channel = await transformer.transform(inx, sc)
                else:
                    resolved_source_channel = None
                await command.callback(
                    cog,
                    interaction=inx,
                    target_channel=resolved_target_channel,
                    source_channel=resolved_source_channel,
                )
            for kind, content in end_results:
                assert kind == "success"
