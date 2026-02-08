"""!/bin/python3
test_hero.py The purpose of this program is to be the integration test for the hero promotion command for Vyrtuous.

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

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, target",
    [
        ("Guild Owner", "!hero", "{member_snowflake}"),
        ("Guild Owner", "!hero", "<@{member_snowflake}>"),
    ],
)
async def test_hero(bot, command: str, target, permission_role):
    """
    Promote or demote member to 'Hero' in memory (lost on reload).

    Parameters
    ----------
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is will be or is a hero
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------

    >>> !hero <@10000000000000003>
    [{emoji} Invincibility granted for Member1]

    >>> !hero 10000000000000003
    [{emoji} Invincibility granted for Member1]
    """
    t = target.format(
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
    go_commands = bot.get_cog("GuildOwnerTextCommands")
    with ExitStack() as stack:
        stack.enter_context(
            patch(
                "vyrtuous.db.roles.owner.guild_owner_service.guild_owner_predicator",
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
            command = await go_commands.invincibility_text_command(ctx, member=t)
        for kind, content in end_results:
            assert kind == "success"
