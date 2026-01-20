"""test_aroles.py The purpose of this program is to be the integration test for the aroles list command for Vyrtuous.

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

GUILD_SNOWFLAKE = 10000000000000500
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!aroles all"),
        ("!aroles {guild_snowflake}"),
    ],
)
async def test_aroles(bot, command: Optional[str]):
    """
    List roles registered in the PostgresSQL database
    'vyrtuous' in the table 'administrator roles'.

    Parameters
    ----------
    all : str, optional
        Generic showing all administrator roles in all guilds
    guild_snowflake : int | str, optional
        Snowflake of a guild where administrator roles are present.

    Examples
    --------
    >>> !aroles "all"
    [{emoji} Administrator Roles\n Guild1\n Guild2]

    >>> !aroles 10000000000000500
    [{emoji} Administrator Roles\n Guild1]
    """
    formatted = command.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
        guild_snowflake=GUILD_SNOWFLAKE
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
