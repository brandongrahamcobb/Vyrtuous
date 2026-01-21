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

from typing import Optional

import pytest

from vyrtuous.tests.integration.test_suite import send_message

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
ROLE_SNOWFLAKE = 10000000000000200
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        # ("!alias ban testban {channel_snowflake}"),
        # ("!testban {member_snowflake}"),
        # ("!xalias testban"),
        # ("!alias unban testunban {channel_snowflake}"),
        # ("!testunban {member_snowflake}"),
        # ("!xalias testunban"),
        # ("!alias vmute testmute {channel_snowflake}"),
        # ("!testmute {member_snowflake}"),
        # ("!xalias testmute"),
        # ("!alias unvmute testunmute {channel_snowflake}"),
        # ("!testunmute {member_snowflake}"),
        # ("!xalias testunmute"),
        # ("!alias flag testflag {channel_snowflake}"),
        # ("!testflag {member_snowflake}"),
        # ("!xalias testflag"),
        # ("!alias unflag testunflag {channel_snowflake}"),
        # ("!testunflag {member_snowflake}"),
        # ("!xalias testunflag"),
        # ("!alias vegan testvegan {channel_snowflake}"),
        # ("!testvegan {member_snowflake}"),
        # ("!xalias testvegan"),
        # ("!alias carnist testcarnist {channel_snowflake}"),
        # ("!testcarnist {member_snowflake}"),
        # ("!xalias testcarnist"),
        # ("!alias tmute testtmute {channel_snowflake}"),
        # ("!testtmute {member_snowflake}"),
        # ("!xalias testtmute"),
        # ("!alias untmute testuntmute {channel_snowflake}"),
        # ("!testuntmute {member_snowflake}"),
        # ("!xalias testuntmute"),
        ("!alias role testrole {channel_snowflake} {role_snowflake}"),
        ("!testrole {member_snowflake}"),
        ("!xalias role testrole"),
        ("!alias unrole testunrole {channel_snowflake} {role_snowflake}"),
        ("!testunrole {member_snowflake}"),
        ("!xalias unrole testunrole"),
    ],
)
async def test_aliases(bot, command: Optional[str]):
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
    formatted = command.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
        role_snowflake=ROLE_SNOWFLAKE,
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured
