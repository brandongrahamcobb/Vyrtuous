"""!/bin/python3
test_streams.py The purpose of this program is to be the integration test for the streams list command for Vyrtuous.

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

from vyrtuous.db.stream import Stream
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    send_message,
    setup,
)

GUILD_SNOWFLAKE = 10000000000000500
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
OTHER_GUILD_CHANNEL_SNOWFLAKE = 10000000000000013
COMMAND = "stream"
BASE_PERMISSIONS = ["command.channel.stream"]
TABLE_NAME = Stream.__tablename__


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target_channel, source, extra_permissions",
    [
        ("{target_channel_snowflake}", None, []),
        ("{target_channel_snowflake}", "{source_channel_snowflake}", []),
        (
            "{target_channel_snowflake}",
            "{other_guild_channel_snowflake}",
            ["other_guilds"],
        ),
        ("{target_channel_snowflake}", "{other_guild_snowflake}", ["other_guilds"]),
        ("<#{target_channel_snowflake}>", None, []),
        ("<#{target_channel_snowflake}>", "<#{source_channel_snowflake}>", []),
        (
            "<#{target_channel_snowflake}>",
            "{other_guild_channel_snowflake}",
            ["other_guilds"],
        ),
        ("<#{target_channel_snowflake}>", "{other_guild_snowflake}", ["other_guilds"]),
    ],
)
async def test_stream_text_command(
    bot,
    prefix: str,
    target_channel: str,
    source: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Setup, modify or teardown a streaming route, modifying the
    the PostgresSQL database 'vyrtuous' in the table 'streaming'.

    Parameters
    ----------
    target_channel : int | str
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    source (Optional) : int | str
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel |  discord.Guild
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> !stream 10000000000000010
    Embed
    """
    assert COMMAND in docstring
    assert TABLE_NAME in docstring
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
            tc = target_channel.format(
                target_channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
            )
            full = f"{prefix}{COMMAND} {tc}"
            if source is None:
                s = source
            else:
                s = source.format(
                    other_guild_snowflake=GUILD_SNOWFLAKE,
                    other_guild_channel_snowflake=OTHER_GUILD_CHANNEL_SNOWFLAKE,
                    source_channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                )
                full += f" {s}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target_channel, source, extra_permissions",
    [
        ("{target_channel_snowflake}", None, []),
        ("{target_channel_snowflake}", "{source_channel_snowflake}", []),
        (
            "{target_channel_snowflake}",
            "{other_guild_channel_snowflake}",
            ["other_guilds"],
        ),
        ("{target_channel_snowflake}", "{other_guild_snowflake}", ["other_guilds"]),
        ("<#{target_channel_snowflake}>", None, []),
        ("<#{target_channel_snowflake}>", "<#{source_channel_snowflake}>", []),
        (
            "<#{target_channel_snowflake}>",
            "{other_guild_channel_snowflake}",
            ["other_guilds"],
        ),
        ("<#{target_channel_snowflake}>", "{other_guild_snowflake}", ["other_guilds"]),
    ],
)
async def test_stream_app_command(
    bot, target_channel: str, source: str | None, extra_permissions: list[str]
):
    docstring = """
    Setup, modify or teardown a streaming route, modifying the
    the PostgresSQL database 'vyrtuous' in the table 'streaming'.

    Parameters
    ----------
    target_channel : int | str
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    source (Optional) : int | str
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel |  discord.Guild
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> /stream 10000000000000010
    Embed
    """
    assert COMMAND in docstring
    assert TABLE_NAME in docstring
    if (
        os.environ["TEST_MODE"].lower() == "app"
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
            cog = bot.get_cog("ChannelManagementAppCommands")
            command = cog.modify_streaming_app_command
            tc = target_channel.format(
                target_channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
            )
            if source is None:
                s = source
            else:
                s = source.format(
                    other_guild_snowflake=GUILD_SNOWFLAKE,
                    other_guild_channel_snowflake=OTHER_GUILD_CHANNEL_SNOWFLAKE,
                    source_channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                )
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
                transformer = AppTarget()
                resolved_target_channel = await transformer.transform(inx, tc)
                if s:
                    resolved_source = await transformer.transform(inx, s)
                else:
                    resolved_source = s
                await command.callback(
                    cog,
                    interaction=inx,
                    target_channel=resolved_target_channel,
                    source=resolved_source,
                )
            for kind, content in end_results:
                assert kind == "success"
