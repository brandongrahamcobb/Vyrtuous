"""!/bin/python3
test_cogs.py The purpose of this program is to be the integration test for the cogs list command for Vyrtuous.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import os
from contextlib import ExitStack
from unittest.mock import patch

import pytest

from vyrtuous.tests.integration.test_suite import check_permissions, send_message

GUILD_SNOWFLAKE = 10000000000000500
COMMAND = "xalias"
BASE_PERMISSIONS = ["command.alias.delete", "command.alias.scope.channel"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "alias_name, guild, extra_permissions",
    [
        ("testban", None, []),
        ("testflag", None, []),
        ("testtmute", None, []),
        ("testmute", None, []),
        ("testban", "{guild_snowflake}", []),
        ("testflag", "{guild_snowflake}", []),
        ("testtmute", "{guild_snowflake}", []),
        ("testmute", "{guild_snowflake}", []),
        (
            "testrole",
            "{guild_snowflake}",
            ["command.alias.scope.role"],
        ),
        (
            "testrole",
            None,
            ["command.alias.scope.role"],
        ),
    ],
)
async def test_xalias_text_command(
    bot,
    prefix: str,
    alias_name: str,
    guild: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Delete command aliases in the PostgreSQL
    database 'vyrtuous' in the table 'command_aliases'.

    Parameters
    ----------
    alias_name : str
        Examples: testban, testmute

    guild (Optional) : str | int
        Resolves to: discord.Guild
        Examples: 10000000000000010

    Examples
    --------
    >>> !xalias testban
    Embed
    """
    assert COMMAND in docstring
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        extra_permissions.extend(BASE_PERMISSIONS)
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions",
                    side_effect=check_permissions(extra_permissions),
                )
            )
            full = f"{prefix}{COMMAND} {alias_name}"
            if guild is None:
                g = guild
            else:
                g = guild.format(
                    guild_snowflake=GUILD_SNOWFLAKE,
                )
                full += f" {g}"
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]
