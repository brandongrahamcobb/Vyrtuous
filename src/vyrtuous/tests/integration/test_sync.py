"""test_sync.py The purpose of this program is to be the integration test for the sync list command for Vyrtuous.

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


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!sync"),
        ("!sync *"),
        ("!sync ^"),
        ("!sync ~"),
    ],
)
async def test_sync(bot, command: Optional[str]):
    """
    List sync loaded by 'Vyrtuous'.

    Parameters
    ----------
    spec
        Syncs app commands globally (None), syncs to the current guild (~),
        syncs to from global to the current guild (*), cleans and syncs to the current guild (^)

    Examples
    --------
    >>> !sync
    [{emoji} Synced # commands globally]

    >>> !sync *
    [{emoji} Synced # commands to the current guild]

    >>> !sync ~
    [{emoji} Synced # commands to the current guild]

    >>> !sync ^
    [{emoji} Synced 0 commands to the current guild]
    """
    captured = await send_message(bot=bot, content=command)
    assert captured
