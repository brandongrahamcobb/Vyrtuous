"""!/bin/python3
test_streams.py The purpose of this program is to be the integration test for the streams list command for Vyrtuous.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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

GUILD_SNOWFLAKE = 10000000000000500
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, target, source",
    [
        (
            "Administrator",
            "stream",
            "{target_channel_snowflake}",
            None,
        ),
        (
            "Administrator",
            "stream",
            "{target_channel_snowflake}",
            "{source_channel_snowflake}",
        ),
        (
            "Administrator",
            "stream",
            "{target_channel_snowflake}",
            "{source_guild_snowflake}",
        ),
        (
            "Administrator",
            "stream",
            "<#{target_channel_snowflake}>",
            None,
        ),
        (
            "Administrator",
            "stream",
            "<#{target_channel_snowflake}>",
            None,
        ),
        (
            "Administrator",
            "stream",
            "<#{target_channel_snowflake}>",
            "{source_guild_snowflake}",
        ),
        (
            "Administrator",
            "stream",
            "{target_channel_snowflake}",
            "{source_channel_snowflake}",
        ),
    ],
)
async def test_stream(bot, command: str, prefix: str, source, target, permission_role):
    """
    Setup, modify or teardown a streaming route, modifying the
    the PostgresSQL database 'vyrtuous' in the table 'streaming'.

    Parameters
    ----------
    source_channel_snowflake : int | str
        Mention or snowflake of a channel with streams
        in any of the guilds Vyrtuous has access inside.
    target_channel_snowflake : int | str, optional
        Mention or snowflake of a channel with streams
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !stream 1000000000000010 create all
    [{emoji} Streaming Route created for Channel11]

    >>> !stream <@10000000000000010> create all
    [{emoji} Streaming Route created for Channel1]

    >>> !stream 10000000000000010 modify channel {channel_snowflake}
    [{emoji} Streaming Route modified for Channel1]
    """
    tc = target.format(
        target_channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
    )
    full = f"{prefix}{command} {tc}"
    s = None
    if source:
        s = source.format(
            source_guild_snowflake=GUILD_SNOWFLAKE,
            source_channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
        )
        full = f"{prefix}{command} {tc} {s}"
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
            cog = bot.get_cog("HiddenAdministratorAppCommands")
            command = cog.modify_streaming_app_command
            transformer = AppTarget()
            resolved_target = await transformer.transform(inx, tc)
            if s:
                resolved = await transformer.transform(inx, s)
            else:
                resolved = None
            await command.callback(
                cog,
                interaction=inx,
                target_channel=resolved_target,
                source=resolved,
            )
        for kind, content in end_results:
            assert kind == "success"
