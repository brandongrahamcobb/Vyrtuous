"""!/bin/python3
test_load_reload_unload.py The purpose of this program is to be the integration test for the load, reload, and unload cog commands for Vyrtuous.

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
from vyrtuous.models.module import AppModule
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, cog",
    [
        ("Developer", "unload", "{cog}"),
        ("Developer", "load", "{cog}"),
        ("Developer", "reload", "{cog}"),
    ],
)
async def test_load_reload_unload(bot, command: str, prefix: str, cog, permission_role):
    """
    Load, reload or unload cogs.

    Parameters
    ----------
    cog
        The cog file path starting with vyrtuous.cog.*

    Examples
    --------
    >>> !load vyrtuous.cog.scheduled_tasks
    [{emoji} Loaded ScheduledTasks]

    >>> !reload vyrtuous.cog.scheduled_tasks
    [{emoji} Reloaded ScheduledTasks]

    >>> !unload vyrtuous.cog.scheduled_tasks
    [{emoji} Unloaded ScheduledTasks]

    """

    permission_state = bot.registry.get(PermissionState)
    c = cog.format(cog="vyrtuous.listeners.scheduled_tasks")
    full = f"{prefix}{command} {c}"
    if command == "load":
        await bot.unload_extension(c)
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
            cog = bot.get_cog("DevelopmentAppCommands")
            transformer = AppModule()
            resolved = await transformer.transform(inx, c)
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
                if command == "reload":
                    command = cog.reload_app_command
                    await command.callback(cog, interaction=inx, module=resolved)
                    for kind, content in end_results:
                        assert kind == "success"
                elif command == "load":
                    await bot.unload_extension(c)
                    command = cog.load_app_command
                    await command.callback(cog, interaction=inx, module=resolved)
                    for kind, content in end_results:
                        assert kind == "success"
                elif command == "unload":
                    await bot.load_extension(c)
                    command = cog.unload_app_command
                    await command.callback(cog, interaction=inx, module=resolved)
            for kind, content in end_results:
                assert kind == "success"
