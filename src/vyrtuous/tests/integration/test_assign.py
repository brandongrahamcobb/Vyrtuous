"""test_assign.py The purpose of this program is to be the integration test for the assign command for Vyrtuous.

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
        ("!assign {uuid}"),
    ],
)
async def test_assign(bot, command: Optional[str]):
    """
    Assign a developer to a developer issues present in the PostgresSQL database
    'vyrtuous' in the table 'developer_logs'.

    Parameters
    ----------
    uuid : str, optional
        UUID reference of a developer issue in
        any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !assign 10000000000000500
    [{emoji} Member1 assign to Developer Issues\n Reference: {UUID}\n Message: {message.jump_url}]
    """
    formatted = command.format(
        uuid=UUID,
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured
