"""test_devs.py The purpose of this program is to be the integration test for the devs list command for Vyrtuous.

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
NOT_PRIVILEGED_AUTHOR_SNOWFLAKE = 10000000000000002


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!devs all"),
        ("!devs {member_snowflake}"),
        ("!devs {guild_snowflake}"),
    ],
)
async def test_devs(bot, command: Optional[str]):
    """
    List members who are registered in the PostgresSQL database
    'vyrtuous' in the table 'developers'.

    Parameters
    ----------
    all : str, optional
        Generic showing all developers in all guilds
    guild_snowflake : int, optional
        Snowflake of a guild where developers are present.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is an developer
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !devs "all"
    [{emoji} Developers\n Guild1\n Guild2]

    >>> !devs 10000000000000500
    [{emoji} Developers\n Guild1]

    >>> !devs <@10000000000000002>
    [{emoji} Developers for Member1\n Guild1\n Guild2]

    >>> !devs 10000000000000002
    [{emoji} Developers for Member1\n Guild1\n Guild2]
    """
    formatted = command.format(
        member_snowflake=NOT_PRIVILEGED_AUTHOR_SNOWFLAKE,
        guild_snowflake=GUILD_SNOWFLAKE
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
