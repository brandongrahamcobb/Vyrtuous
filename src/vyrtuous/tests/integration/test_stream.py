"""!/bin/python3
test_streams.py The purpose of this program is to be the integration test for the streams list command for Vyrtuous.

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

from contextlib import ExitStack
from typing import Optional
from unittest.mock import patch

import pytest

from vyrtuous.tests.conftest import context
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)

TEXT_CHANNEL_SNOWFLAKE = 10000000000000010
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, source_channel, action, type, target",
    [
        (
            "Administrator",
            "!stream",
            "{source_channel_snowflake}",
            "create",
            "all",
            None,
        ),
        (
            "Administrator",
            "!stream",
            "{source_channel_snowflake}",
            "modify",
            "channel",
            "{target_channel_snowflake}",
        ),
        (
            "Administrator",
            "!stream",
            "{source_channel_snowflake}",
            "delete",
            None,
            None,
        ),
        (
            "Administrator",
            "!stream",
            "<#{source_channel_snowflake}>",
            "create",
            "all",
            None,
        ),
        (
            "Administrator",
            "!stream",
            "<#{source_channel_snowflake}>",
            "modify",
            None,
            None,
        ),
        (
            "Administrator",
            "!stream",
            "<#{source_channel_snowflake}>",
            "delete",
            None,
            None,
        ),
        (
            "Administrator",
            "!stream",
            "{source_channel_snowflake}",
            "create",
            "channel",
            "{target_channel_snowflake}",
        ),
        (
            "Administrator",
            "!stream",
            "{source_channel_snowflake}",
            "create",
            "channel",
            "Fail",
        ),
    ],
)
async def test_stream(
    bot, command: str, source_channel, action, target, type, permission_role
):
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
    snowflakes = None
    sc = source_channel.format(source_channel_snowflake=TEXT_CHANNEL_SNOWFLAKE)
    kwargs = {"channel": sc, "action": action}
    full = f"{command} {sc} {action}"
    if type:
        kwargs.update({"entry_type": type})
        full = f"{command} {sc} {action} {type}"
    if target:
        tc = target.format(
            target_channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
        )
        snowflakes = [tc]
        full = f"{command} {sc} {action} {type} {tc}"
    # captured = await send_message(bot=bot, content=full)
    # assert captured
    objects = setup(bot)
    msg = build_message(
        author=objects.get("author", None),
        channel=objects.get("text_channel", None),
        content=full,
        guild=objects.get("guild", None),
        state=objects.get("state", None),
    )
    ctx = context(
        bot=bot,
        channel=objects.get("text_channel", None),
        guild=objects.get("guild", None),
        message=msg,
        prefix="!",
    )
    admin_commands = bot.get_cog("AdminTextCommands")
    with ExitStack() as stack:
        stack.enter_context(
            patch(
                "vyrtuous.db.roles.admin.administrator_service.administrator_predicator",
                return_value=True,
            )
        )
        stack.enter_context(
            patch(
                "vyrtuous.commands.permissions.permission_service.PermissionService.has_equal_or_lower_role",
                return_value=permission_role,
            )
        )
        stack.enter_context(
            patch(
                "vyrtuous.commands.permissions.permission_service.PermissionService.resolve_highest_role",
                return_value=permission_role,
            )
        )
        async with capture_command() as end_results:
            command = await admin_commands.modify_streaming_text_command(ctx, **kwargs)
        for kind, content in end_results:
            assert kind == "success"
