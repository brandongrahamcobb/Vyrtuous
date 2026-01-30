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

from vyrtuous.tests.integration.conftest import context
from vyrtuous.tests.integration.test_suite import build_message, send_message, setup


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, spec",
    [
        ("!sync", None),
        ("!sync", "*"),
        ("!sync", "^"),
        ("!sync", "~"),
    ],
)
async def test_sync(bot, command: Optional[str], spec):
    """
    Syncs app commands.

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
    if spec:
        full = f"{command} {spec}"
    else:
        full = f"{command}"
    captured = await send_message(bot=bot, content=full)
    assert captured
    objects = setup(bot)
    msg = build_message(
        author=objects.get("author", None),
        channel=objects.get("channel", None),
        content=full,
        guild=objects.get("guild", None),
        state=objects.get("state", None),
    )
    ctx = context(bot=bot, message=msg, prefix="!")
    dev_commands = bot.get_cog("DevCommands")
    command = await dev_commands.sync_text_command(ctx, spec=spec)
