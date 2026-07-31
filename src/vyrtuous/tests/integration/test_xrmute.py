"""!/bin/python3
test_rmute_xrmute.py The purpose of this program is to be the integration test for the rmute and xrmute commands for Vyrtuous.

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
OTHER_GUILD_CHANNEL_SNOWFLAKE = 10000000000000013

BASE_PERMISSIONS = [
    "command.moderation.unvoice-mute.channel_unmute",
]
COMMAND = "xrmute"


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "channel, extra_permissions",
    [
        (None, []),
        (None, []),
        ("{channel_snowflake}", []),
        ("<#{channel_snowflake}>", []),
        ("{channel_snowflake}", []),
        ("<#{channel_snowflake}>", []),
        ("{channel_snowflake}", []),
        ("<#{channel_snowflake}>", []),
        ("{other_guild_channel_snowflake}", ["other_guilds"]),
        ("<#{other_guild_channel_snowflake}>", ["other_guilds"]),
        ("{other_guild_channel_snowflake}", ["other_guilds"]),
        ("<#{other_guild_channel_snowflake}>", ["other_guilds"]),
        ("{other_guild_channel_snowflake}", ["other_guilds"]),
        ("<#{other_guild_channel_snowflake}>", ["other_guilds"]),
    ],
)
async def test_xrmute_text_command(
    bot,
    prefix: str,
    channel: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Simulate a right-click unmute for all member in the voice channel except for the executor
    and members with higher group privileges.

    Parameters
    ----------
    channel (Optional) : str | int | None
        Resolves to: discord.VoiceChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> !xrmute
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
            full = f"{prefix}{COMMAND}"
            if channel is None:
                c = channel
            else:
                c = channel.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    other_guild_channel_snowflake=OTHER_GUILD_CHANNEL_SNOWFLAKE,
                )
                full += f" {c}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "channel, extra_permissions",
    [
        (None, []),
        (None, []),
        ("{channel_snowflake}", []),
        ("<#{channel_snowflake}>", []),
        ("{channel_snowflake}", []),
        ("<#{channel_snowflake}>", []),
        ("{channel_snowflake}", []),
        ("<#{channel_snowflake}>", []),
        ("{other_guild_channel_snowflake}", ["other_guilds"]),
        ("<#{other_guild_channel_snowflake}>", ["other_guilds"]),
        ("{other_guild_channel_snowflake}", ["other_guilds"]),
        ("<#{other_guild_channel_snowflake}>", ["other_guilds"]),
        ("{other_guild_channel_snowflake}", ["other_guilds"]),
        ("<#{other_guild_channel_snowflake}>", ["other_guilds"]),
    ],
)
async def test_xrmute_app_command(
    bot,
    channel: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Simulate a right-click unmute for all member in the voice channel except for the executor
    and members with higher group privileges.

    Parameters
    ----------
    channel (Optional) : str | int | None
        Resolves to: discord.VoiceChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> /xrmute
    Embed
    """
    assert COMMAND in docstring
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
            cog = bot.get_cog("ModerationAppCommands")
            command = cog.channel_unmute_app_command
            if channel is None:
                c = channel
            else:
                c = channel.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    other_guild_channel_snowflake=OTHER_GUILD_CHANNEL_SNOWFLAKE,
                )
            objects = setup(bot)
            msg = build_message(
                author=objects.get("author", None),
                channel=objects.get("voice_channel", None),
                content="",
                guild=objects.get("guild", None),
                state=objects.get("state", None),
            )
            inx = interaction(
                bot=bot,
                channel=objects.get("voice_channel", None),
                guild=objects.get("guild", None),
                message=msg,
            )
            async with capture_command() as end_results:
                target_transformer = AppTarget()
                if c:
                    resolved_channel = await target_transformer.transform(inx, c)
                else:
                    resolved_channel = c
                await command.callback(
                    cog,
                    interaction=inx,
                    channel=resolved_channel,
                )
            for kind, content in end_results:
                assert kind == "success"
