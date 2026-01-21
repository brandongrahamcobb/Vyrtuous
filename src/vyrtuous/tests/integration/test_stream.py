"""test_streams.py The purpose of this program is to be the integration test for the streams list command for Vyrtuous.

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
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!stream {source_channel_snowflake} create all"),
        (
            "!stream {source_channel_snowflake} modify channel {target_channel_snowflake}"
        ),
        ("!stream {source_channel_snowflake} delete"),
        ("!stream <#{source_channel_snowflake}> create all"),
        ("!stream <#{source_channel_snowflake}> modify"),
        ("!stream <#{source_channel_snowflake}> delete"),
        (
            "!stream {source_channel_snowflake} create channel {target_channel_snowflake}"
        ),
    ],
)
async def test_stream(bot, command: Optional[str]):
    """
    Setup, modify or teardown a streaming route, modifying the
    the PostgresSQL database 'vyrtuous' in the table 'streaming'.

    Parameters
    ----------
    source_channel_snowflake : int | str
        Mention or snowflake of a channel with streams
        in any of the guilds Vyrtuous has access inside.
    target_channel_snowflake : int | str, optional
        Mention or snowflake of a channel with streams
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !stream 1000000000000010 create all
    [{emoji} Streaming Route created for Channel11]

    >>> !stream <@10000000000000010> create all
    [{emoji} Streaming Route created for Channel1]

    >>> !stream 10000000000000010 modify channel {channel_snowflake}
    [{emoji} Streaming Route modified for Channel1]
    """
    formatted = command.format(
        source_channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
        target_channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured
