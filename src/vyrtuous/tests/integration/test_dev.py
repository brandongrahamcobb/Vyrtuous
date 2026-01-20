"""test_dev.py The purpose of this program is to be the integration test for the dev promotion command for Vyrtuous.

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


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!dev {member_snowflake}"),
        ("!dev <@{member_snowflake}>"),
    ],
)
async def test_dev(bot, command: Optional[str]):
    """
    Promote a member to 'Developer' by registering them in the PostgresSQL database
    'vyrtuous' in the table 'developers'.

    Parameters
    ----------
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is an developer
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------

    >>> !dev <@10000000000000002>
    [{emoji} Developer granted for Member1]

    >>> !dev 10000000000000002
    [{emoji} Developer granted for Member1]
    """
    formatted = command.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
