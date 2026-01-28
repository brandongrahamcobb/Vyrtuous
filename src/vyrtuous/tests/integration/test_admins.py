"""test_admins.py The purpose of this program is to be the integration test for the admins list command for Vyrtuous.

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
from typing import Optional

import discord
import pytest

from vyrtuous.tests.integration.conftest import context
from vyrtuous.tests.integration.test_suite import (
    build_message,
    send_message,
    setup
)

GUILD_SNOWFLAKE = 10000000000000500
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
NOT_PRIVILEGED_AUTHOR_SNOWFLAKE_ONE = 10000000000000002
NOT_PRIVILEGED_AUTHOR_NAME_ONE = "Not Privileged Author Name One"


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, target",
    [
        ("!admins", "all"),
        ("!admins", "{guild_snowflake}"),
        ("!admins", "{member_snowflake}"),
        ("!admins", "<@{member_snowflake}>"),
    ],
)
async def test_admins(bot, command: Optional[str], target):
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
    formatted = target.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE, guild_snowflake=GUILD_SNOWFLAKE
    )
    full = f"{command} {formatted}"
    captured = await send_message(bot=bot, content=full)
    assert captured.content
    objects = setup(bot)
    msg = build_message(
        author=objects.get("author", None), channel=objects.get("channel", None), content=full, guild=objects.get("guild", None), state=objects.get("state", None)
    )
    ctx = context(bot=bot, message=msg, prefix="!")
    mod_commands = bot.get_cog("ModeratorCommands")
    command = await mod_commands.list_administrators_text_command(ctx, target=formatted)
