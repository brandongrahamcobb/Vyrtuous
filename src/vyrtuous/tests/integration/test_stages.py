"""!/bin/python3
test_stages.py The purpose of this program is to be the integration test for the stages list command for Vyrtuous.

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

GUILD_SNOWFLAKE = 10000000000000500
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, target",
    [
        ("Sysadmin", "!stages", "all"),
        ("Administrator", "!stages", "{channel_snowflake}"),
        ("Administrator", "!stages", "<#{channel_snowflake}>"),
        ("Administrator", "!stages", "{guild_snowflake}"),
    ],
)
async def test_stages(bot, command: str, target, permission_role):
    """
    List channels which are registered in the PostgresSQL database
    'vyrtuous' in the table 'stages'.

    Parameters
    ----------
    all : str, optional
        Generic showing all temporary rooms in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with stages
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where temporary rooms are present.

    Examples
    --------
    >>> !stages "all"
    [{emoji} Stages\n Guild1\n Guild2]

    >>> !stages 10000000000000500
    [{emoji} Stages\n Guild1]

    >>> !stages <@10000000000000010>
    [{emoji} Stages for Channel1]

    >>> !stages 10000000000000010
    [{emoji} Stages for Channel1]
    """
    t = target.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE, guild_snowflake=GUILD_SNOWFLAKE
    )
    full = f"{command} {t}"
    captured = await send_message(bot=bot, content=full)
    assert captured.content
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
            command = await admin_commands.list_stages_text_command(ctx, target=t)
        for kind, content in end_results:
            assert kind == "success"
