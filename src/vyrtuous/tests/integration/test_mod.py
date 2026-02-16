"""!/bin/python3
test_mod.py The purpose of this program is to be the integration test for the mod demotion/promotion command for Vyrtuous.

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
from vyrtuous.tests.integration.test_suite import (build_message,
                                                   capture_command,
                                                   send_message, setup)

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, member, channel",
    [
        ("Coordinator", "!mod", "{member_snowflake}", "{channel_snowflake}"),
        ("Coordinator", "!mod", "<@{member_snowflake}>", "<#{channel_snowflake}>"),
    ],
)
async def test_mod(bot, command: str, member, channel, permission_role):
    """
    Promote or demote a member with 'Moderator' by registering them in the PostgresSQL database
    'vyrtuous' in the table 'moderators'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with modorators
        in any of the guilds Vyrtuous has access inside.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is an moderator
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------

    >>> !mod <@10000000000000003> <@10000000000000010>
    [{emoji} Moderator granted for Member1]

    >>> !mod 10000000000000003 10000000000000010
    [{emoji} Moderator granted for Member1]
    """
    m = member.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
    )
    c = channel.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
    )
    full = f"{command} {m} {c}"
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
        coord_commands = bot.get_cog("CoordinatorTextCommands")
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.owner.guild_owner_service.guild_owner_predicator",
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
                command = await coord_commands.toggle_moderator_text_command(
                    ctx, member=m, channel=c
                )
            for kind, content in end_results:
                assert kind == "success"
