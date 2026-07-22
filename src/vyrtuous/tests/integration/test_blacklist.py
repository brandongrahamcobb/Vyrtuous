"""!/bin/python3
test_mod.py The purpose of this program is to be the integration test for the mod demotion/promotion command for Vyrtuous.

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

GUILD_SNOWFLAKE = 10000000000000500
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, member, channel",
    [
        (
            "Coordinator",
            "blacklist",
            "{member_snowflake}",
            "{channel_snowflake}",
        ),
        (
            "Coordinator",
            "blacklist",
            "<@{member_snowflake}>",
            None,
        ),
        (
            "Coordinator",
            "blacklist",
            "{simplified_member_snowflake}",
            "{channel_snowflake}",
        ),
        (
            "Coordinator",
            "blacklist",
            "<@{simplified_member_snowflake}>",
            "<#{channel_snowflake}>",
        ),
        (
            "Coordinator",
            "blacklist",
            "<@{member_snowflake}>",
            "<#{channel_snowflake}>",
        ),
    ],
)
async def test_blacklist(
    bot, command: str, prefix: str, member, channel, permission_role
):
    """
    Blacklist or unlisted a member's ban in the PostgresSQL database
    'vyrtuous' in the table 'active_bans'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with ban
        in any of the guilds Vyrtuous has access inside.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is banned
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------

    >>> !blacklist <@10000000000000003> <@10000000000000010>
    [{emoji} Member blacklisted in Channel]

    >>> !mod 10000000000000003 10000000000000010
    [{emoji} Member unlisted in Channel]
    """
    bot.registry.get(MemberState).active.update(
        {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
    )
    m = member.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
        simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
    )
    full = f"{prefix}{command} {m}"
    c = None
    if channel:
        c = channel.format(
            channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
        )
        full = f"{prefix}{command} {m} {c}"
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
            cog = bot.get_cog("HiddenCoordinatorAppCommands")
            command = cog.toggle_blacklist_app_command
            transformer = AppTarget()
            resolved_member = await transformer.transform(inx, m)
            if c:
                resolved_channel = await transformer.transform(inx, c)
            else:
                resolved_channel = None
            await command.callback(
                cog,
                interaction=inx,
                member=resolved_member,
                channel=resolved_channel,
            )
        for kind, content in end_results:
            assert kind == "success"
