"""test_vr.py The purpose of this program is to be the integration test for the vr command for Vyrtuous.

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
        ("!vr {channel_snowflake}"),
        ("!vr {channel_snowflake}"),
    ],
)
async def test_vr(bot, command: Optional[str]):
    """
    Create or teardown a video room by accessing
    the PostgresSQL database 'vyrtuous' in the table 'video_rooms'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with vr
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !vr <@10000000000000010>
    [{emoji} Video Room has been created]

    >>> !vr 10000000000000010
    [{emoji} Video Rooms has been deleted]
    """
    formatted = command.format(channel_snowflake=TEXT_CHANNEL_SNOWFLAKE)
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
