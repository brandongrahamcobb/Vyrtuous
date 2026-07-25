"""!/bin/python3
test_cap.py The purpose of this program is to be the integration test for the cap command for Vyrtuous.

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

import pytest

from vyrtuous.models.category import AppCategory
from vyrtuous.models.duration import AppDuration
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)

VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, category, channel, limit",
    [
        ("Administrator", "cap", "ban", None, None),
        ("Administrator", "cap", "ban", "{channel_snowflake}", None),
        ("Administrator", "cap", "ban", "{channel_snowflake}", "8h"),
    ],
)
async def test_cap(bot, command: str, prefix: str, channel, category, limit, permission_role):
    """
    Set a expires in limit for a channel by populating the PostgresSQL database
    'vyrtuous' in the table 'active_caps'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with cap
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !cap 10000000000000010 ban 8
    [{emoji} Ban cap created\n Guild1\n Channel1]
    """
    c = None
    if channel:
        c = channel.format(
            channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
        )
    l = None
    if limit:
        l = limit
    full = f"{prefix}{command} {category}"
    if c and not limit:
        full = f"{prefix}{command} {category} {c}"
    elif c and limit:
        full = f"{prefix}{command} {category} {c} {limit}"
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
            cog = bot.get_cog("HiddenAdministratorAppCommands")
            command = cog.cap_app_command
            category_transformer = AppCategory()
            resolved_category = await category_transformer.transform(inx, category)
            target_transformer = AppTarget()
            if c:
                resolved_channel = await target_transformer.transform(inx, c)
            else:
                resolved_channel = None
            duration_transformer = AppDuration()
            if l:
                resolved_limit = await duration_transformer.transform(inx, l)
            else:
                resolved_limit = None
            await command.callback(
                cog,
                interaction=inx,
                category=resolved_category,
                channel=resolved_channel,
                limit=resolved_limit,
            )
        for kind, content in end_results:
            assert kind == "success"
