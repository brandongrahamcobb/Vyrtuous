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
from unittest.mock import AsyncMock, patch

import pytest

from vyrtuous.cache.permissions import PermissionGroup, PermissionScope
from vyrtuous.cache.registry import PermissionState
from vyrtuous.tests.integration.test_suite import send_message
from vyrtuous.utils.permissions import permission_service

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
ROLE_SNOWFLAKE = 10000000000000200
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, category, alias_name, channel_snowflake, role_snowflake",
    [
        ("Administrator", "alias", "vmute", "testmute", "{channel_snowflake}", None),
        ("Administrator", "alias", "flag", "testflag", "{channel_snowflake}", None),
        (
            "Administrator",
            "alias",
            "tmute",
            "testtmute",
            "{channel_snowflake}",
            "{role_snowflake}",
        ),
        (
            "Administrator",
            "alias",
            "ban",
            "testban",
            "{channel_snowflake}",
            "{role_snowflake}",
        ),
        (
            "Administrator",
            "alias",
            "role",
            "testrole",
            "{channel_snowflake}",
            "{role_snowflake}",
        ),
    ],
)
async def test_alias(
    bot,
    command: str,
    prefix: str,
    category,
    alias_name,
    channel_snowflake,
    role_snowflake,
    permission_role,
):
    """
    Create and delete command aliases in the PostgreSQL
    database 'vyrtuous' in the table 'command_aliases'.

    Parameters
    ----------
    alias_type
        The type of alias. Can be one of ban, unban, vmute, unvmute
        flag, unflag, vegan, carnist, tmute, untmute, role and unrole.
    alias_name
        The name of the alias.
    role_snowflake
        The snowflake or mention of a role

    Examples
    --------
    >>> !alias ban testban
    [{emoji} Alias `testban` created]

    >>> !testban 10000000000000003
    [{emoji} Member Name was Banned]

    >>> !xalias testban
    [{emoji} Alias `testban` deleted]
    """
    permission_state = bot.registry.get(PermissionState)
    channel = channel_snowflake.format(
        channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
    )
    kwargs = {"category": category, "alias_name": alias_name, "channel": channel}
    if role_snowflake:
        role = role_snowflake.format(
            role_snowflake=ROLE_SNOWFLAKE,
        )
        kwargs.update({"role": role})
        full = f"{prefix}{command} {category} {alias_name} {channel} {role}"
    else:
        full = f"{prefix}{command} {category} {alias_name} {channel}"
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.utils.permissions.permission_service.resolve_effective_group",
                    new=AsyncMock(
                        return_value=permission_state.groups.get(
                            permission_role.lower()
                        )
                    ),
                )
            )
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]
