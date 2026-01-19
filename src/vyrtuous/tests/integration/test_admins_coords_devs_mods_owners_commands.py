"""test_chown_temp_temps_xtemp_commands.py The purpose of this program is to black box test the temporary room commands.
Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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

from vyrtuous.tests.integration.conftest import bot
from vyrtuous.tests.integration.test_suite import send_message

NOT_PRIVILEGED_AUTHOR_ID = 10000000000000002

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!admins"),
    ],
)
async def test_admins_coords_devs_mods_owners_commands(bot, command: Optional[str]):
    captured = await send_message(bot=bot, content="!admins")
    print(captured.content)
