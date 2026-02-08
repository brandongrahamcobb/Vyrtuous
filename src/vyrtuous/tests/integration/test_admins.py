"""!/bin/python3
test_admins.py The purpose of this program is to be the integration test for the admins list command for Vyrtuous.

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

import inspect
from contextlib import ExitStack
from typing import Optional
from unittest.mock import patch

import discord
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
NOT_PRIVILEGED_AUTHOR_SNOWFLAKE_ONE = 10000000000000002
NOT_PRIVILEGED_AUTHOR_NAME_ONE = "Not Privileged Author Name One"


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, target",
    [
        ("Sysadmin", "!admins", "all"),
        ("Administrator", "!admins", "{guild_snowflake}"),
        ("Moderator", "!admins", "{member_snowflake}"),
        ("Moderator", "!admins", "<@{member_snowflake}>"),
    ],
)
async def test_admins(bot, command: str, target, permission_role):
    """
    List members who are registered in the PostgresSQL database
    'vyrtuous' in the table 'administrators'.

    Parameters
    ----------
    all : str, optional
        Generic showing all administrators in all guilds
    guild_snowflake : int | str, optional
        Snowflake of a guild where administrators are present.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is an administrator
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !admins "all"
    [{emoji} Administrators\n Guild1\n Guild2]

    >>> !admins 10000000000000500
    [{emoji} Administrators\n Guild1]

    >>> !admins <@10000000000000003>
    [{emoji} Administrators for Member1\n Guild1\n Guild2]

    >>> !admins 10000000000000003
    [{emoji} Administrators for Member1\n Guild1\n Guild2]
    """
    t = target.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE, guild_snowflake=GUILD_SNOWFLAKE
    )
    full = f"{command} {t}"
    # captured = await send_message(bot=bot, content=full)
    # assert captured
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
                "vyrtuous.db.roles.mod.moderator_service.moderator_predicator",
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
            command = await mod_commands.list_administrators_text_command(ctx, target=t)
        for kind, content in end_results:
            print(kind, content)
            assert kind == "success"
