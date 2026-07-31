"""!/bin/python3
test_cogs.py The purpose of this program is to be the integration test for the cogs list command for Vyrtuous.

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

import importlib
import inspect
import os
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

import pytest
from discord.ext import commands

from vyrtuous.inc.helpers import DISCORD_COGS, DISCORD_COGS_CLASSES
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    send_message,
    setup,
)

COMMAND = "cogs"
BASE_PERMISSIONS = ["command.info.cogs"]


@pytest.mark.asyncio
async def test_cogs_text_command(bot, prefix: str):
    docstring = """
    List cogs.

    Parameters
    ----------
    None

    Examples
    --------
    >>> !cogs
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
            full = f"{prefix}{COMMAND}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
async def test_cogs_app_command(bot):
    docstring = """
    List cogs.

    Parameters
    ----------
    None

    Examples
    --------
    >>> !cogs
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
            cog = bot.get_cog("InfoAppCommands")
            command = cog.list_cogs_app_command
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
                await command.callback(cog, interaction=inx)
            for kind, content in end_results:
                assert kind == "success"


def test_all_cogs():
    cogs = []
    base = Path.cwd().parent
    for path in base.rglob("*.py"):
        module = importlib.import_module(
            str(path.with_suffix("").relative_to(base)).replace("/", ".")
        )
        for _, obj in inspect.getmembers(module, inspect.isclass):
            if issubclass(obj, commands.Cog) and obj is not commands.Cog:
                cogs.append(obj)
    assert cogs
    for cog in cogs:
        assert cog.__module__ in DISCORD_COGS
        assert cog.__name__ in DISCORD_COGS_CLASSES
