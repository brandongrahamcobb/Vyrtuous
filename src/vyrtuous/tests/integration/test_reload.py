"""!/bin/python3
test_reload.py The purpose of this program is to be the test for the reload cog command.

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

from vyrtuous.models.module import AppModule
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    send_message,
    setup,
)

COG = "vyrtuous.listeners.scheduled_tasks"
COMMAND = "reload"
BASE_PERMISSIONS = ["command.dev.reload"]


@pytest.mark.asyncio
async def test_reload_text_command(bot, prefix: str):
    docstring = """
    Reload a cog.

    Parameters
    ----------
    module : str
        Resolves to: ModuleObject
        Examples: vyrtuous.listeners.scheduled_tasks

    Examples
    --------
    >>> !reload
    {emoji} Reloaded `vyrtuous.listeners.scheduled_tasks`
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
            full = f"{prefix}{COMMAND} {COG}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
async def test_reload_app_command(bot):
    docstring = """
    Reload a cog.

    Parameters
    ----------
    module : str
        Resolves to: ModuleObject
        Examples: vyrtuous.listeners.scheduled_tasks

    Examples
    --------
    >>> !reload
    {emoji} Reloaded `vyrtuous.listeners.scheduled_tasks`
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
            cog = bot.get_cog("DevelopmentAppCommands")
            command = cog.reload_app_command
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
            module_transformer = AppModule()
            resolved_module = await module_transformer.transform(inx, COG)
            async with capture_command() as end_results:
                await command.callback(cog, interaction=inx, module=resolved_module)
            for kind, content in end_results:
                assert kind == "success"
