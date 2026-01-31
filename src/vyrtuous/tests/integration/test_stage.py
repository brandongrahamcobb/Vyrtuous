"""test_stage.py The purpose of this program is to be the integration test for the stage command for Vyrtuous.

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

TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, channel",
    [
        ("!stage", "{channel_snowflake}"),
        ("!stage", "{channel_snowflake}"),
    ],
)
async def test_stage(bot, command: str, channel):
    """
    Create or teardown a stage by accessing
    the PostgresSQL database 'vyrtuous' in the table 'video_rooms'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with stage
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !stage <@10000000000000010>
    [{emoji} Stage Room has been created]

    >>> !stage 10000000000000010
    [{emoji} Stage has been ended]
    """
    c = channel.format(channel_snowflake=TEXT_CHANNEL_SNOWFLAKE)
    full = f"{command} {c}"
    captured = await send_message(bot=bot, content=full)
    assert captured.content
    objects = setup(bot)
    msg = build_message(
        author=objects.get("author", None),
        channel=objects.get("text_channel", None),
        content=full,
        guild=objects.get("guild", None),
        state=objects.get("state", None),
    )
    ctx = context(
        bot=bot,
        channel=objects.get("text_channel", None),
        guild=objects.get("guild", None),
        message=msg,
        prefix="!",
    )
    admin_commands = bot.get_cog("AdminTextCommands")
    command = await admin_commands.toggle_stage_text_command(ctx, channel=c)
