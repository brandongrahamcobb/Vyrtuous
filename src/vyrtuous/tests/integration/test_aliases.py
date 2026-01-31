"""test_cogs.py The purpose of this program is to be the integration test for the cogs list command for Vyrtuous.

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

import asyncio
from typing import Optional

import pytest

from vyrtuous.tests.integration.conftest import context
from vyrtuous.tests.integration.test_suite import build_message, send_message, setup

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
ROLE_SNOWFLAKE = 10000000000000200
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, member_snowflake",
    [
        ("!testban", "{member_snowflake}"),
        ("!testmute", "{member_snowflake}"),
        ("!testflag", "{member_snowflake}"),
        ("!testvegan", "{member_snowflake}"),
        ("!testtmute", "{member_snowflake}"),
        ("!testrole", "{member_snowflake}"),
    ],
)
async def test_aliases(bot, command: str, member_snowflake):
    """
    Create and delete command aliases in the PostgreSQL
    database 'vyrtuous' in the table 'command_aliases'.

    Parameters
    ----------
    alias_type
        The type of alias. Can be one of ban, unban, vmute, unvmute
        flag, unflag, vegan, carnist, tmute, untmute, role and unrole.
    alias_name
        The name of the alias.
    role_snowflake
        The snowflake or mention of a role

    Examples
    --------
    >>> !alias ban testban
    [{emoji} Alias `testban` created]

    >>> !testban 10000000000000003
    [{emoji} Member Name was Banned]

    >>> !xalias testban
    [{emoji} Alias `testban` deleted]
    """
    member = member_snowflake.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
    )
    full = f"{command} {member}"
    captured = await send_message(bot=bot, content=full)
    assert captured
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
    generic_event_listeners = bot.get_cog("GenericEventListeners")
    command = await generic_event_listeners.on_message(message=msg)
