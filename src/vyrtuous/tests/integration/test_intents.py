"""!/bin/python3
test_pc.py The purpose of this program is to be the integration test for the pc list command for Vyrtuous.

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

GUILD_SNOWFLAKE = 10000000000000500
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
COMMAND = "intents"
BASE_PERMISSIONS = ["command.info.intents"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target, extra_permissions",
    [
        ("{channel_snowflake}", ["command.info.scope.channel"]),
        ("<#{channel_snowflake}>", ["command.info.scope.channel"]),
        ("{guild_snowflake}", ["command.info.scope.guild"]),
    ],
)
async def test_intents_text_command(
    bot, prefix: str, target: str | None, extra_permissions
):
    docstring = """
    List member and role intents for a channel or guild.

    Parameters
    ----------
    target (Optional) : str | int
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel |  discord.Guild
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> !intents
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
            if target is None:
                t = target
            else:
                t = target.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    guild_snowflake=GUILD_SNOWFLAKE,
                )
                full += f" {t}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target, extra_permissions",
    [
        ("{channel_snowflake}", ["command.info.scope.channel"]),
        ("<#{channel_snowflake}>", ["command.info.scope.channel"]),
        ("{guild_snowflake}", ["command.info.scope.guild"]),
    ],
)
async def test_intents_app_command(bot, target: str | None, extra_permissions):
    docstring = """
    List member and role intents for a channel or guild.

    Parameters
    ----------
    target (Optional) : str | int
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel |  discord.Guild
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> /intents
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
            cog = bot.get_cog("InfoAppCommands")
            command = cog.list_intents_app_command
            if target is None:
                t = target
            else:
                t = target.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    guild_snowflake=GUILD_SNOWFLAKE,
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
                if t:
                    resolved_target = await transformer.transform(inx, t)
                else:
                    resolved_target = None
                await command.callback(
                    cog,
                    interaction=inx,
                    target=resolved_target,
                )
            for kind, content in end_results:
                assert kind == "success"
