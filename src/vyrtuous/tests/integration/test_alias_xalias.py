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

ROLE_SNOWFLAKE = 10000000000000200
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!alias ban testban {channel_snowflake}"),
        ("!xalias testban"),
        ("!alias unban testunban {channel_snowflake}"),
        ("!xalias testunban"),
        ("!alias vmute testmute {channel_snowflake}"),
        ("!xalias testmute"),
        ("!alias unvmute testunmute {channel_snowflake}"),
        ("!xalias testunmute"),
        ("!alias flag testflag {channel_snowflake}"),
        ("!xalias testflag"),
        ("!alias unflag testunflag {channel_snowflake}"),
        ("!xalias testunflag"),
        ("!alias vegan testvegan {channel_snowflake}"),
        ("!xalias testvegan"),
        ("!alias carnist testcarnist {channel_snowflake}"),
        ("!xalias testcarnist"),
        ("!alias tmute testtmute {channel_snowflake}"),
        ("!xalias testtmute"),
        ("!alias untmute testuntmute {channel_snowflake}"),
        ("!xalias testuntmute"),
        ("!alias role testrole {channel_snowflake} {role_snowflake}"),
        ("!xalias role testrole"),
        ("!alias unrole testunrole {channel_snowflake} {role_snowflake}"),
        ("!xalias unrole testunrole"),
    ],
)
async def test_cogs(bot, command: Optional[str]):
    """
    Create and delete command aliases in the PostgreSQL
    database 'vyrtuous' in the table 'command_aliases'.

    Parameters
    ----------
    alias_type
        The type of alias. Can be one of ban, unban, voice_mute, unvoice_mute
        flag, unflag, vegan, carnist, tmute, untmute, role and unrole.
    alias_name
        The name of the alias.

    Examples
    --------
    >>> !alias ban testban
    [{emoji} Alias `testban` created]

    >>> !xalias testban
    [{emoji} Alias `testban` deleted]
    """
    formatted = command.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE, role_snowflake=ROLE_SNOWFLAKE
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
