"""!/bin/python3
test_coords.py The purpose of this program is to be the integration test for the coords list command for Vyrtuous.

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
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, target",
    [
        ("Sysadmin", "!coords", "all"),
        ("Moderator", "!coords", "{channel_snowflake}"),
        ("Moderator", "!coords", "<#{channel_snowflake}>"),
        ("Administrator", "!coords", "{guild_snowflake}"),
        ("Moderator", "!coords", "{member_snowflake}"),
        ("Moderator", "!coords", "<@{member_snowflake}>"),
    ],
)
async def test_coords(bot, command: str, target, permission_role):
    """
    List members who are registered in the PostgresSQL database
    'vyrtuous' in the table 'coordinators'.

    Parameters
    ----------
    all : str, optional
        Generic showing all coordinators in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with coordinators
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where coordinators are present.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is an coordinator
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !coords "all"
    [{emoji} Coordinators\n Guild1\n Guild2]

    >>> !coords <#10000000000000010>
    [{emoji} Coordinators for Channel1\n Member1\n Member2]

    >>> !coords 10000000000000010
    [{emoji} Coordinators for Channel1\n Member1\n Member2]

    >>> !coords 10000000000000500
    [{emoji} Coordinators\n Guild1]

    >>> !coords <@10000000000000003>
    [{emoji} Coordinators for Member1\n Guild1\n Guild2]

    >>> !coords 10000000000000003
    [{emoji} Coordinators for Member1\n Guild1\n Guild2]
    """
    t = target.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
        guild_snowflake=GUILD_SNOWFLAKE,
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
    )
    full = f"{command} {t}"
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
    mod_commands = bot.get_cog("ModeratorTextCommands")
    with ExitStack() as stack:
        stack.enter_context(
            patch(
                "vyrtuous.moderator.moderator_service.moderator_predicator",
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
            command = await mod_commands.list_coordinators_text_command(ctx, target=t)
        for kind, content in end_results:
            assert kind == "success"
