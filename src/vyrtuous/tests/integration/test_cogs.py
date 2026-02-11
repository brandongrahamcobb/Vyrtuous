"""!/bin/python3
test_cogs.py The purpose of this program is to be the integration test for the cogs list command for Vyrtuous.

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


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command",
    [
        ("Developer", "!cogs"),
    ],
)
async def test_cogs(bot, command: str, permission_role):
    """
    List cogs loaded by 'Vyrtuous'.

    Parameters
    ----------
    None
        Generic showing all cogs

    Examples
    --------
    >>> !cogs
    [{emoji} Cogs\n Cog1\n Cog2]
    """
    # captured = await send_message(bot=bot, content=command)
    # assert captured
    objects = setup(bot)
    msg = build_message(
        author=objects.get("author", None),
        channel=objects.get("text_channel", None),
        content=command,
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
    dev_commands = bot.get_cog("DevTextCommands")
    with ExitStack() as stack:
        stack.enter_context(
            patch(
                "vyrtuous.developer.developer_service.developer_predicator",
                return_value=True,
            )
        )
        stack.enter_context(
            patch(
                "vyrtuous.utils.permission_service.PermissionService.has_equal_or_lower_role",
                return_value=permission_role,
            )
        )
        stack.enter_context(
            patch(
                "vyrtuous.utils.permission_service.PermissionService.resolve_highest_role",
                return_value=permission_role,
            )
        )
        async with capture_command() as end_results:
            command = await dev_commands.list_cogs_text_command(ctx)
        for kind, content in end_results:
            assert kind == "success"
