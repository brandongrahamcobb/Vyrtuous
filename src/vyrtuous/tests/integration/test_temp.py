"""test_temp.py The purpose of this program is to be the integration test for the temp command for Vyrtuous.

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

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, channel",
    [
        ("Administrator", "!temp", "{channel_snowflake}"),
        ("Administrator", "!temp", "{channel_snowflake}"),
    ],
)
async def test_temp(bot, command: str, channel, permission_role):
    """
    Create or teardown a temporary room by accessing
    the PostgresSQL database 'vyrtuous' in the table 'temporary_rooms'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with temp
        in any of the guilds Vyrtuous has access inside.
    mebmer_snowflake : int | str, optional
        Mention or snowflake of a member to own the channel.

    Examples
    --------
    >>> !temp <@10000000000000010> 10000000000000003
    [{emoji} Temporary Room has been created]

    >>> !temp 10000000000000010
    [{emoji} Temporary Rooms has been deleted]
    """
    c = channel.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
    )
    full = f"{command} {c}"
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
            command = await admin_commands.toggle_temp_room_text_command(ctx, channel=c)
        for kind, content in end_results:
            assert kind == "success"
