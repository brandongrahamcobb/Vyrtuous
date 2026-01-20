"""test_coords.py The purpose of this program is to be the integration test for the coords list command for Vyrtuous.

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
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!coords all"),
        ("!coords {guild_snowflake}"),
        ("!coords {member_snowflake}"),
        ("!coords <@{member_snowflake}>"),
    ],
)
async def test_coords(bot, command: Optional[str]):
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

    >>> !coords <@10000000000000002>
    [{emoji} Coordinators for Member1\n Guild1\n Guild2]

    >>> !coords 10000000000000002
    [{emoji} Coordinators for Member1\n Guild1\n Guild2]
    """
    formatted = command.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
        guild_snowflake=GUILD_SNOWFLAKE,
        member_snowflake=NOT_PRIVILEGED_AUTHOR_SNOWFLAKE,
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
