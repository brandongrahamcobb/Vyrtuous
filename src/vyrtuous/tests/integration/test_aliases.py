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

from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
ROLE_SNOWFLAKE = 10000000000000200
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, member_snowflake",
    [
        ("Moderator", "testban", "{member_snowflake}"),
        ("Moderator", "testban", "{member_snowflake}"),
        ("Moderator", "testmute", "{member_snowflake}"),
        ("Moderator", "testmute", "{member_snowflake}"),
        ("Moderator", "testflag", "{member_snowflake}"),
        ("Moderator", "testflag", "{member_snowflake}"),
        ("Moderator", "testtmute", "{member_snowflake}"),
        ("Moderator", "testtmute", "{member_snowflake}"),
        ("Coordinator", "testrole", "{member_snowflake}"),
    ],
)
async def test_aliases(bot, prefix, command: str, member_snowflake, permission_role):
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
    member = member_snowflake.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
    )
    full = f"{prefix}{command} {member}"
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        objects = setup(bot)
        msg = build_message(
            author=objects.get("author", None),
            channel=objects.get("text_channel", None),
            content=full,
            guild=objects.get("guild", None),
            state=objects.get("state", None),
        )
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.tests.integration.mock_discord_channel.MockVoiceChannel.fetch_message",
                    new=AsyncMock(return_value=msg),
                )
            )
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]
