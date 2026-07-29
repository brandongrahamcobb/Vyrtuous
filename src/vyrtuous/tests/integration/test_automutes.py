"""!/bin/python3
test_stages.py The purpose of this program is to be the integration test for the stages list command for Vyrtuous.

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

from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)

GUILD_SNOWFLAKE = 10000000000000500
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, target",
    [
        ("Administrator", "automutes", "{channel_snowflake}"),
        ("Administrator", "automutes", "<#{channel_snowflake}>"),
        ("Administrator", "automutes", "{guild_snowflake}"),
    ],
)
async def test_ams(bot, command: str, prefix: str, target, permission_role):
    """
    List channels which are registered in the PostgresSQL database
    'vyrtuous' in the table 'stages'.

    Parameters
    ----------
    all : str, optional
        Generic showing all temporary channels in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with stages
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where temporary channels are present.

    Examples
    --------
    >>> !stages "all"
    [{emoji} Stages\n Guild1\n Guild2]

    >>> !stages 10000000000000500
    [{emoji} Stages\n Guild1]

    >>> !stages <@10000000000000010>
    [{emoji} Stages for Channel1]

    >>> !stages 10000000000000010
    [{emoji} Stages for Channel1]
    """
    t = target.format(
        channel_snowflake=VOICE_CHANNEL_SNOWFLAKE, guild_snowflake=GUILD_SNOWFLAKE
    )
    full = f"{prefix}{command} {t}"
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
            cog = bot.get_cog("InfoAppCommands")
            command = cog.list_automute_channels_app_command
            transformer = AppTarget()
            resolved_target = await transformer.transform(inx, t)
            await command.callback(
                cog,
                interaction=inx,
                target=resolved_target,
            )
        for kind, content in end_results:
            assert kind == "success"
