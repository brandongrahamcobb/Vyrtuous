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
    "command, category, alias_name, channel_snowflake, role_snowflake",
    [
        ("!alias", "role", "testrole", "fail_channel", "fail_role"),
        ("!alias", "role", "testrole", "fail_channel", None),
        ("!alias", "vmute", "testmute", "{channel_snowflake}", None),
        ("!alias", "flag", "testflag", "{channel_snowflake}", None),
        ("!alias", "vegan", "testvegan", "{channel_snowflake}", None),
        ("!alias", "tmute", "testtmute", "{channel_snowflake}", "{role_snowflake}"),
        ("!alias", "ban", "testban", "{channel_snowflake}", "{role_snowflake}"),
        ("!alias", "role", "testrole", "{channel_snowflake}", "{role_snowflake}"),
    ],
)
async def test_alias(
    bot, command: str, category, alias_name, channel_snowflake, role_snowflake
):
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
    channel = channel_snowflake.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
    )
    kwargs = {"category": category, "alias_name": alias_name, "channel": channel}
    if role_snowflake:
        role = role_snowflake.format(
            role_snowflake=ROLE_SNOWFLAKE,
        )
        kwargs.update({"role": role})
        full = f"{command} {category} {alias_name} {channel} {role}"
    else:
        full = f"{command} {category} {alias_name} {channel}"
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
    admin_commands = bot.get_cog("AdminTextCommands")
    command = await admin_commands.create_alias_text_command(ctx, **kwargs)
