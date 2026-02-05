"""test_roleid.py The purpose of this program is to be the integration test for the roleid list command for Vyrtuous.

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

from vyrtuous.tests.integration.conftest import context
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
    "permission_role, command, source_channel, target_channel",
    [
        (
            "Administrator",
            "!rmv",
            "{source_channel_snowflake}",
            "{target_channel_snowflake}",
        ),
    ],
)
async def test_rmv(bot, command: str, source_channel, target_channel, permission_role):
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
    sc = source_channel.format(
        source_channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
    )
    tc = target_channel.format(
        target_channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
    )
    full = f"{command} {sc} {tc}"
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
            command = await admin_commands.room_move_all_text_command(
                ctx, source_channel=sc, target_channel=tc
            )
        for kind, content in end_results:
            assert kind == "success"
