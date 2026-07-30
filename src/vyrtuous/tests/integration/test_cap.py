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
from contextlib import ExitStack
from unittest.mock import AsyncMock, patch

import pytest

from vyrtuous.db.cap import Cap
from vyrtuous.models.category import AppCategory
from vyrtuous.models.duration import AppDuration
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    send_message,
    setup,
)

VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
COMMAND = "cap"
BASE_PERMISSIONS = ["command.channel.cap"]
TABLE_NAME = Cap.__tablename__


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "channel, category, duration",
    [
        ("{channel_snowflake}", "ban", None),
        ("{channel_snowflake}", "tmute", None),
        ("{channel_snowflake}", "vmute", None),
        ("<#{channel_snowflake}>", "ban", None),
        ("<#{channel_snowflake}>", "tmute", None),
        ("<#{channel_snowflake}>", "vmute", None),
        ("{channel_snowflake}", "ban", "8h"),
        ("{channel_snowflake}", "tmute", "8h"),
        ("{channel_snowflake}", "vmute", "8h"),
        ("<#{channel_snowflake}>", "ban", "8h"),
        ("<#{channel_snowflake}>", "tmute", "8h"),
        ("<#{channel_snowflake}>", "vmute", "8h"),
    ],
)
async def test_cap_text_command(
    bot, prefix: str, category: str, channel: str, duration: str | None
):
    docstring = """
    Set a moderation duration limit for capped members in a channel by
    populating the PostgresSQL database 'vyrtuous' in the table 'active_caps'.

    Parameters
    ----------
    channel : int | str
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    category : str
        Resolves to: CategoryObject
        Examples: ban, tmute or vmute

    duration : str
        Resolves to: DurationObject
        Examples: 0 | 1 | 30m | 1h | 1d | 1w 

    Example
    --------
    >>> !cap 10000000000000010 ban
    "Cap `ban` created for Ask a Vegan successfully."
    """
    assert COMMAND in docstring
    assert TABLE_NAME in docstring
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions",
                    side_effect=check_permissions(BASE_PERMISSIONS),
                )
            )
            c = channel.format(
                channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
            )
            full = f"{prefix}{COMMAND} {category} {c}"
            if duration is None:
                d = duration
            else:
                d = duration
                full += f" {d}"
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "channel, category, duration",
    [
        ("{channel_snowflake}", "ban", None),
        ("{channel_snowflake}", "tmute", None),
        ("{channel_snowflake}", "vmute", None),
        ("<#{channel_snowflake}>", "ban", None),
        ("<#{channel_snowflake}>", "tmute", None),
        ("<#{channel_snowflake}>", "vmute", None),
        ("{channel_snowflake}", "ban", "8h"),
        ("{channel_snowflake}", "tmute", "8h"),
        ("{channel_snowflake}", "vmute", "8h"),
        ("<#{channel_snowflake}>", "ban", "8h"),
        ("<#{channel_snowflake}>", "tmute", "8h"),
        ("<#{channel_snowflake}>", "vmute", "8h"),
    ],
)
async def test_cap_app_command(bot, category: str, channel: str, duration: str | None):
    docstring = """
    Set a moderation duration limit for capped members in a channel by
    populating the PostgresSQL database 'vyrtuous' in the table 'active_caps'.

    Parameters
    ----------
    channel : int | str
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    category : str
        Resolves to: CategoryObject
        Examples: ban, tmute or vmute

    duration : str
        Resolves to: DurationObject
        Examples: 0 | 1 | 30m | 1h | 1d | 1w 

    Example
    --------
    >>> /cap 10000000000000010 ban
    "Cap `ban` created for Ask a Vegan successfully."
    """
    assert COMMAND in docstring
    assert TABLE_NAME in docstring
    if (
        os.environ["TEST_MODE"].lower() == "app"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions",
                    side_effect=check_permissions(BASE_PERMISSIONS),
                )
            )
            cog = bot.get_cog("ChannelManagementAppCommands")
            command = cog.cap_app_command
            c = channel.format(
                channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
            )
            if duration is None:
                d = duration
            else:
                d = duration
            objects = setup(bot)
            msg = build_message(
                author=objects.get("author", None),
                channel=objects.get("text_channel", None),
                content="",
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
                category_transformer = AppCategory()
                duration_transformer = AppDuration()
                target_transformer = AppTarget()
                resolved_channel = await target_transformer.transform(inx, c)
                resolved_category = await category_transformer.transform(inx, category)
                if d:
                    resolved_duration = await duration_transformer.transform(inx, d)
                else:
                    resolved_duration = d
                await command.callback(
                    cog,
                    interaction=inx,
                    channel=resolved_channel,
                    category=resolved_category,
                    limit=resolved_duration,
                )
            for kind, content in end_results:
                assert kind == "success"
