"""test_aroles.py The purpose of this program is to be the integration test for the aroles list command for Vyrtuous.

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

ROLE_SNOWFLAKE = 10000000000000200


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, role",
    [
        ("!admin", "{role_snowflake}"),
        ("!admin", "<@&{role_snowflake}>"),
    ],
)
async def test_admin(bot, command: Optional[str], role):
    """
    List roles registered in the PostgresSQL database
    'vyrtuous' in the table 'administrator roles'.

    Parameters
    ----------
    all : str, optional
        Generic showing all administrator roles in all guilds
    guild_snowflake : int | str, optional
        Snowflake of a guild where administrator roles are present.

    Examples
    --------
    >>> !aroles "all"
    [{emoji} Administrator Roles\n Guild1\n Guild2]

    >>> !aroles 10000000000000500
    [{emoji} Administrator Roles\n Guild1]
    """
    r = role.format(
        role_snowflake=ROLE_SNOWFLAKE,
    )
    full = f"{command} {r}"
    captured = await send_message(bot=bot, content=full)
    assert captured.content
    objects = setup(bot)
    msg = build_message(
        author=objects.get("author", None),
        channel=objects.get("text_channel", None),
        content=full,
        guild=objects.get("guild", None),
        state=objects.get("state", None),
    )
    ctx = context(
        bot=bot,
        channel=objects.get("text_channel", None),
        guild=objects.get("guild", None),
        message=msg,
        prefix="!",
    )
    go_commands = bot.get_cog("GuildOwnerTextCommands")
    command = await go_commands.toggle_administrator_by_role_text_command(ctx, role=r)
