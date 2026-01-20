"""test_rmute_xrmute.py The purpose of this program is to be the integration test for the rmute and xrmute commands for Vyrtuous.

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
        ("!rmute {channel_snowflake}"),
        ("!xrmute {channel_snowflake}"),
        ("!rmute <#{channel_snowflake}>"),
        ("!xrmute <#{channel_snowflake}>"),
    ],
)
async def test_rmute_xrmute(bot, command: Optional[str]):
    """
    Voice-mute a whole room and undo it by adding and removing
    entries in the PostgreSQL database 'vyrtuous' in the table
    'active_voice_mutes'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !rmute 10000000000000010
    [{emoji} Room Muted\n Member1\n Member2]

    >>> !xrmute 10000000000000010
    [{emoji} Room Unmuted\n Member1\n Member2]

    >>> !rmute <#10000000000000010>
    [{emoji} Room Muted\n Member1\n Member2]

    >>> !xrmute <#10000000000000010>
    [{emoji} Room Unmuted\n Member1\n Member2]
    """
    formatted = command.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
