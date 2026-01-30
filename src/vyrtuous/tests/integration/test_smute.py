"""test_smute.py The purpose of this program is to be the integration test for the smute command for Vyrtuous.

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

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, member",
    [
        ("!smute", "{member_snowflake}"),
        ("!smute", "<@{member_snowflake}>"),
    ],
)
async def test_smute(bot, command: Optional[str], member):
    """
    Server mute a member localized to the guild

    Parameters
    ----------
    None
        No paramers neccesary

    Examples
    --------

    >>> !smute <@10000000000000003>
    [{emoji} Member1 was Server Muted]

    >>> !smute 10000000000000003
    [{emoji} Member1 was Server Muted]
    """
    m = member.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
    )
    full = f"{command} {m} test_reason"
    captured = await send_message(bot=bot, content=full)
    assert captured.content
    objects = setup(bot)
    msg = build_message(
        author=objects.get("author", None),
        channel=objects.get("channel", None),
        content=full,
        guild=objects.get("guild", None),
        state=objects.get("state", None),
    )
    ctx = context(bot=bot, message=msg, prefix="!")
    admin_commands = bot.get_cog("AdminCommands")
    command = await admin_commands.toggle_server_mute_text_command(
        ctx, member=m, reason="test_reason"
    )
