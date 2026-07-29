"""!/bin/python3
test_smute.py The purpose of this program is to be the integration test for the smute command for Vyrtuous.

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
GUILD_SNOWFLAKE = 10000000000000500
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, member, guild, reason",
    [
        ("Administrator", "smute", "{member_snowflake}", None, None),
        (
            "Administrator",
            "smute",
            "<@{member_snowflake}>",
            "{guild_snowflake}",
            "test_reason",
        ),
        ("Administrator", "smute", "{simplified_member_snowflake}", None, None),
        (
            "Administrator",
            "smute",
            "<@{simplified_member_snowflake}>",
            "{guild_snowflake}",
            "test_reason",
        ),
        (
            "Administrator",
            "smute",
            "<@{member_snowflake}>",
            "{guild_snowflake}",
            None,
        ),
    ],
)
async def test_smute(
    bot, command: str, prefix: str, member, reason, guild, permission_role
):
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
    bot.registry.get(MemberState).active.update(
        {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
    )

    m = member.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
        simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
    )
    g = None
    r = None
    full = f"{prefix}{command} {m}"
    if reason and not guild:
        r = reason
        full = f"{prefix}{command} {m} {r}"
    if reason and guild:
        r = reason
        g = guild.format(guild_snowflake=GUILD_SNOWFLAKE)
        full = f"{prefix}{command} {m} {g} {r}"
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
            cog = bot.get_cog("ModerationAppCommands")
            command = cog.toggle_server_mute_app_command
            transformer = AppTarget()
            resolved_member = await transformer.transform(inx, m)
            if g:
                resolved_guild = await transformer.transform(inx, g)
            else:
                resolved_guild = None
            await command.callback(
                cog,
                interaction=inx,
                member=resolved_member,
                reason=r,
                guild=resolved_guild,
            )
        for kind, content in end_results:
            assert kind == "success"
