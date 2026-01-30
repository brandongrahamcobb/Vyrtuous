"""test_load_reload_unload.py The purpose of this program is to be the integration test for the load, reload, and unload cog commands for Vyrtuous.

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
    "command, cog",
    [
        ("!unload", "{cog}"),
        ("!load", "{cog}"),
        ("!reload", "{cog}"),
    ],
)
async def test_load_reload_unload(bot, command: Optional[str], cog):
    """
    Load, reload or unload cogs.

    Parameters
    ----------
    cog
        The cog file path starting with vyrtuous.cogs.*

    Examples
    --------
    >>> !load vyrtuous.cogs.scheduled_tasks
    [{emoji} Loaded ScheduledTasks]

    >>> !reload vyrtuous.cogs.scheduled_tasks
    [{emoji} Reloaded ScheduledTasks]

    >>> !unload vyrtuous.cogs.scheduled_tasks
    [{emoji} Unloaded ScheduledTasks]

    """
    c = cog.format(cog="vyrtuous.cogs.scheduled_tasks")
    full = f"{command} {c}"
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
    command = await dev_commands.load_text_command(ctx, module=c)
    command = await dev_commands.reload_text_command(ctx, module=c)
    command = await dev_commands.unload_text_command(ctx, module=c)
