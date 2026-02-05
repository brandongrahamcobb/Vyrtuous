"""test_assign.py The purpose of this program is to be the integration test for the assign command for Vyrtuous.

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

GUILD_SNOWFLAKE = 10000000000000500
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
UUID = "7c772534-9528-4c3d-a065-ad3e29f754f8"


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, reference, action, notes",
    [
        ("Developer", "!bug", "{uuid}", "resolve", None),
        ("Developer", "!bug", "{uuid}", "append", "Test notes"),
        ("Developer", "!bug", "{uuid}", "overwrite", "Test notes"),
    ],
)
async def test_assign(bot, command: str, reference, action, notes, permission_role):
    """
    Assign a developer to a developer issues present in the PostgresSQL database
    'vyrtuous' in the table 'developer_logs'.

    Parameters
    ----------
    uuid : str, optional
        UUID reference of a developer issue in
        any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !assign 10000000000000500
    [{emoji} Member1 assign to Developer Issues\n Reference: {UUID}\n Message: {message.jump_url}]
    """
    ref = reference.format(
        uuid=UUID,
    )
    if notes:
        full = f"{command} {ref} {action} {notes}"
    else:
        full = f"{command} {ref} {action}"
    # captured = await send_message(bot=bot, content=full)
    # assert captured.content
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
    dev_commands = bot.get_cog("DevTextCommands")

    with ExitStack() as stack:
        stack.enter_context(
            patch(
                "vyrtuous.db.roles.dev.developer_service.developer_predicator",
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
            if notes:
                command = await dev_commands.update_bug_tracking_text_command(
                    ctx, reference=ref, action=action, notes=notes
                )
            else:
                command = await dev_commands.update_bug_tracking_text_command(
                    ctx, reference=ref, action=action
                )
        for kind, content in end_results:
            assert kind == "success"
