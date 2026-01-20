"""test_dlogs.py The purpose of this program is to be the integration test for the dlogs list command for Vyrtuous.

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
UUID = "7c772534-9528-4c3d-a065-ad3e29f754f8"

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!dlogs all"),
        ("!dlogs all resolved"),
        ("!dlogs all unresolved"),
        ("!dlogs {uuid}"),
        ("!dlogs {uuid} resolved"),
        ("!dlogs {uuid} unresolved"),
        ("!dlogs {guild_snowflake}"),
        ("!dlogs {guild_snowflake} resolved"),
        ("!dlogs {guild_snowflake} unresolved"),
    ],
)
async def test_cogs(bot, command: Optional[str]):
    """
    List developer issues which are registered in the PostgresSQL database
    'vyrtuous' in the table 'developer_logs'.

    Parameters
    ----------
    all : str, optional
        Generic showing all developer issues in all guilds
    guild_snowflake : int | str, optional
        Snowflake of a guild where developer issues are present.
    uuid : str, optional
        UUID reference of a developer issue in
        any of the guilds Vyrtuous has access inside.
    resolved : str, optional
        Filter which specifies whether the listed developer
        issues are resolved
    unresolved : str, optional
        Filter which specifies whether the listed developer
        issues are unresolved

    Examples
    --------
    >>> !dlogs "all"
    [{emoji} Developer Issues\n Guild1\n Guild2]

    >>> !dlogs "all" "resolved"
    [{emoji} Developer Issues\n Guild1\n Guild2]

    >>> !dlogs "all" "unresolved
    [{emoji} Developer Issues\n Guild1\n Guild2]

    >>> !dlogs 10000000000000500
    [{emoji} Developer Issues\n Guild1]

    >>> !dlogs 10000000000000500 resolved
    [{emoji} Developer Issues\n Guild1]

    >>> !dlogs 10000000000000500 unresolved
    [{emoji} Developer Issues\n Guild1]

    >>> !dlogs "7c772534-9528-4c3d-a065-ad3e29f754f8"
    [{emoji} Developer Issues\n Guild1]

    >>> !dlogs "7c772534-9528-4c3d-a065-ad3e29f754f8" resolved
    [{emoji} Developer Issues\n Guild1]

    >>> !dlogs "7c772534-9528-4c3d-a065-ad3e29f754f8" unresolved
    [{emoji} Developer Issues\n Guild1]
    """
    formatted = command.format(
        guild_snowflake=GUILD_SNOWFLAKE,
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
        uuid=UUID
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured
