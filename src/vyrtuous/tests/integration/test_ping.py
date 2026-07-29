"""!/bin/python3
test_ping.py The purpose of this program is to be the integration test for the ping command for Vyrtuous.

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
from unittest.mock import AsyncMock, patch

import pytest

from vyrtuous.cache.registry import PermissionState
from vyrtuous.tests.conftest import interaction
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
        ("Developer", "ping"),
    ],
)
async def test_ping(bot, command: str, prefix: str, permission_role):
    """
    Backup the database 'vyrtuous'.

    Parameters
    ----------
    None
        No parameter required

    Examples
    --------
    >>> !ping
    [{emoji} Pong!]
    """
    permission_state = bot.registry.get(PermissionState)
    full = f"{prefix}{command}"
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.utils.permissions.permission_service.resolve_effective_group",
                    new=AsyncMock(
                        return_value=permission_state.groups.get(
                            permission_role.lower()
                        )
                    ),
                )
            )
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
            content="",
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
            cog = bot.get_cog("UtilityAppCommands")
            command = cog.ping_app_command
            with ExitStack() as stack:
                stack.enter_context(
                    patch(
                        "vyrtuous.utils.permissions.permission_service.resolve_effective_group",
                        new=AsyncMock(
                            return_value=permission_state.groups.get(
                                permission_role.lower()
                            )
                        ),
                    )
                )
                await command.callback(cog, interaction=inx)
        for kind, content in end_results:
            assert kind == "success"
