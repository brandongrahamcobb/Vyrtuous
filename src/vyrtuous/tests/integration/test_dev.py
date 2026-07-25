"""!/bin/python3
test_dev.py The purpose of this program is to be the integration test for the dev demotion/promotion command for Vyrtuous.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

import os
from datetime import datetime, timezone

import pytest

from vyrtuous.cache.registry import MemberState
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, member",
    [
        ("Sysadmin", "dev", "{member_snowflake}"),
        ("Sysadmin", "dev", "<@{member_snowflake}>"),
        ("Sysadmin", "dev", "{simplified_member_snowflake}"),
        ("Sysadmin", "dev", "<@{simplified_member_snowflake}>"),
    ],
)
async def test_dev(bot, command: str, prefix: str, member, permission_role):
    """
    Promote or demote a member with 'Developer' by registering them in the PostgresSQL database
    'vyrtuous' in the table 'developers'.

    Parameters
    ----------
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is an developer
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------

    >>> !dev <@10000000000000003>
    [{emoji} Developer granted for Member1]

    >>> !dev 10000000000000003
    [{emoji} Developer granted for Member1]
    """
    bot.registry.get(MemberState).active.update(
        {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
    )
    m = member.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
        simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
    )
    full = f"{prefix}{command} {m}"
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        captured = await send_message(bot=bot, content=full)
        assert captured == ["success"]
    if (
        os.environ["TEST_MODE"].lower() == "app"
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
        inx = interaction(
            bot=bot,
            channel=objects.get("text_channel", None),
            guild=objects.get("guild", None),
            message=msg,
        )
        async with capture_command() as end_results:
            cog = bot.get_cog("SysadminAppCommands")
            command = cog.toggle_developer_app_command
            transformer = AppTarget()
            resolved_member = await transformer.transform(inx, argument=m)
            await command.callback(cog, interaction=inx, member=resolved_member)
        for kind, content in end_results:
            assert kind == "success"
