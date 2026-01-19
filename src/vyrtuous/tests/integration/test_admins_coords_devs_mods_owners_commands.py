"""test_command.py The purpose of this program is to be the introductory integration test for Vyrtuous.

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

NOT_PRIVILEGED_AUTHOR_ID = 10000000000000002

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!admins 10000000000000001"),
    ],
)
async def test_admins_coords_devs_mods_owners_commands(bot, command: Optional[str]):
    captured = await send_message(bot=bot, content="!admins")
    print(captured.content)
