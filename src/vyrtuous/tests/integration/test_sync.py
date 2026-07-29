"""!/bin/python3
test_sync.py The purpose of this program is to be the integration test for the sync list command for Vyrtuous.

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
    "permission_role, command, spec",
    [
        ("Guild_Owner", "sync", None),
        ("Guild_Owner", "sync", "*"),
        ("Guild_Owner", "sync", "^"),
        ("Guild_Owner", "sync", "~"),
    ],
)
async def test_sync(bot, command: str, prefix: str, spec, permission_role):
    """
    Syncs app commands.

    Parameters
    ----------
    spec
        Syncs app commands globally (None), syncs to the current guild (~),
        syncs to from global to the current guild (*), cleans and syncs to the current guild (^)

    Examples
    --------
    >>> !sync
    [{emoji} Synced # commands globally]

    >>> !sync *
    [{emoji} Synced # commands to the current guild]

    >>> !sync ~
    [{emoji} Synced # commands to the current guild]

    >>> !sync ^
    [{emoji} Synced 0 commands to the current guild]
    """
    permission_state = bot.registry.get(PermissionState)
    if spec:
        full = f"{prefix}{command} {spec}"
    else:
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
    elif (
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
            cog = bot.get_cog("GuildOwnerAppCommands")
            command = cog.sync_app_command
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
                await command.callback(cog, interaction=inx, spec=spec)
            for kind, content in end_results:
                assert kind == "success"
