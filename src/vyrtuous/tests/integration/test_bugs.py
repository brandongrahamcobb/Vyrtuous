"""test_bugs.py The purpose of this program is to be the integration test for the bugs list command for Vyrtuous.

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

GUILD_SNOWFLAKE = 10000000000000500
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
UUID = "7c772534-9528-4c3d-a065-ad3e29f754f8"


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, target, filter",
    [
        ("!bugs", "all", None),
        ("!bugs", "all", "resolved"),
        ("!bugs", "all", "unresolved"),
        ("!bugs", "{uuid}", None),
        ("!bugs", "{uuid}", "resolved"),
        ("!bugs", "{uuid}", "unresolved"),
        ("!bugs", "{guild_snowflake}", None),
        ("!bugs", "{guild_snowflake}", "resolved"),
        ("!bugs", "{guild_snowflake}", "unresolved"),
    ],
)
async def test_bugs(bot, command: Optional[str], target, filter):
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
    >>> !bugs "all"
    [{emoji} Developer Issues\n Guild1\n Guild2]

    >>> !bugs "all" "resolved"
    [{emoji} Developer Issues\n Guild1\n Guild2]

    >>> !bugs "all" "unresolved
    [{emoji} Developer Issues\n Guild1\n Guild2]

    >>> !bugs 10000000000000500
    [{emoji} Developer Issues\n Guild1]

    >>> !bugs 10000000000000500 resolved
    [{emoji} Developer Issues\n Guild1]

    >>> !bugs 10000000000000500 unresolved
    [{emoji} Developer Issues\n Guild1]

    >>> !bugs "7c772534-9528-4c3d-a065-ad3e29f754f8"
    [{emoji} Developer Issues\n Guild1]

    >>> !bugs "7c772534-9528-4c3d-a065-ad3e29f754f8" resolved
    [{emoji} Developer Issues\n Guild1]

    >>> !bugs "7c772534-9528-4c3d-a065-ad3e29f754f8" unresolved
    [{emoji} Developer Issues\n Guild1]
    """
    formatted = target.format(
        guild_snowflake=GUILD_SNOWFLAKE,
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
        uuid=UUID,
    )
    if filter:
        full = f"{command} {formatted} {filter}"
    else:
        full = f"{command} {formatted}"
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
    command = await dev_commands.list_bugs_text_command(
        ctx, target=formatted, filter=filter
    )
