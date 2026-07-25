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
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, target_channel, source_channel",
    [
        (
            "Administrator",
            "rmv",
            "{target_channel_snowflake}",
            "{source_channel_snowflake}",
        ),
        (
            "Administrator",
            "rmv",
            "{target_channel_snowflake}",
            None,
        ),
    ],
)
async def test_rmv(
    bot, command: str, prefix: str, source_channel, target_channel, permission_role
):
    """
    Move all members from one VC to another

    Parameters
    ----------
    source_channel_snowflake
        The snowflake or mention of a channel

    target_channel_snowflake
        The snowflake or mention of a channel

    Examples
    --------
    >>> !!rmv 1000000000000010 1000000000000011
    [{emoji} Members moved succesfully to Voice Channel One\n Member1\b Member2]

    """
    tc = target_channel.format(
        target_channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
    )
    full = f"{prefix}{command} {tc}"
    sc = None
    if source_channel:
        sc = source_channel.format(
            source_channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
        )
        full = f"{prefix}{command} {sc} {tc}"
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
            cog = bot.get_cog("AdministratorAppCommands")
            command = cog.channel_move_all_app_command
            transformer = AppTarget()
            resolved_target = await transformer.transform(inx, tc)
            if sc:
                resolved_source = await transformer.transform(inx, sc)
            else:
                resolved_source = None
            await command.callback(
                cog,
                interaction=inx,
                target_channel=resolved_target,
                source_channel=resolved_source,
            )
        for kind, content in end_results:
            assert kind == "success"
