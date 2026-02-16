"""!/bin/python3
test_aroles.py The purpose of this program is to be the integration test for the aroles list command for Vyrtuous.

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

import os
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

GUILD_SNOWFLAKE = 10000000000000500
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, specifier",
    [
        ("Sysadmin", "!aroles", "all"),
        ("Administrator", "!aroles", "{guild_snowflake}"),
    ],
)
async def test_aroles(bot, command: str, specifier, permission_role):
    """
    List roles registered in the PostgresSQL database
    'vyrtuous' in the table 'administrator roles'.

    Parameters
    ----------
    all : str, optional
        Generic showing all administrator roles in all guilds
    guild_snowflake : int | str, optional
        Snowflake of a guild where administrator roles are present.

    Examples
    --------
    >>> !aroles "all"
    [{emoji} Administrator Roles\n Guild1\n Guild2]

    >>> !aroles 10000000000000500
    [{emoji} Administrator Roles\n Guild1]
    """
    formatted = specifier.format(guild_snowflake=GUILD_SNOWFLAKE)
    full = f"{command} {formatted}"
    if os.environ["TEST_MODE"].lower() == "integration":
        captured = await send_message(bot=bot, content=full)
        assert captured
    elif os.environ["TEST_MODE"].lower() == "unit":
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
                    "vyrtuous.administrator.administrator_service.administrator_predicator",
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
                command = await admin_commands.list_administrator_roles_text_command(
                    ctx, target=formatted
                )
            for kind, content in end_results:
                assert kind == "success"
