"""test_caps.py The purpose of this program is to be the integration test for the caps list command for Vyrtuous.

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
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!caps all"),
        ("!caps {channel_snowflake}"),
        ("!caps <#{channel_snowflake}>"),
        ("!caps {guild_snowflake}"),
    ],
)
async def test_caps(bot, command: Optional[str]):
    """
    List caps in the PostgresSQL database
    'vyrtuous' in the table 'active_caps'.

    Parameters
    ----------
    all : str, optional
        Generic showing all caps in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with caps
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where caps are present.

    Examples
    --------
    >>> !caps "all"
    [{emoji} Caps\n Guild1\n Guild2]

    >>> !caps 10000000000000500
    [{emoji} Caps\n Guild1]

    >>> !caps <@10000000000000002>
    [{emoji} Caps for Channel1]

    >>> !caps 10000000000000002
    [{emoji} Caps for Channel1]
    """
    formatted = command.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
        guild_snowflake=GUILD_SNOWFLAKE
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
