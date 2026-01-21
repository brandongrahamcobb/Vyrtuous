"""test_survey.py The purpose of this program is to be the integration test for the survey command for Vyrtuous.

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
        ("!survey {channel_snowflake} "),
        ("!survey <#{channel_snowflake}>"),
    ],
)
async def test_survey(bot, command: Optional[str]):
    """
    Server mute a member localized to the guild

    Parameters
    ----------
    None
        No paramers neccesary

    Examples
    --------

    >>> !survey <#10000000000000010>
    [{emoji} Survey results for Channel1]

    >>> !survey 10000000000000010
    [{emoji} Survey results for Channel1]
    """
    formatted = command.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
