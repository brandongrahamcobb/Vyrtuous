"""test_cap.py The purpose of this program is to be the integration test for the cap command for Vyrtuous.

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

TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!cap {channel_snowflake} ban 8"),
    ],
)
async def test_cap(bot, command: Optional[str]):
    """
    Set a expires in limit for a channel by populating the PostgresSQL database
    'vyrtuous' in the table 'active_caps'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with cap
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !cap 10000000000000010 ban 8
    [{emoji} Ban cap created\n Guild1\n Channel1]
    """
    formatted = command.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
