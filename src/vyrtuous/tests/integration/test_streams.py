"""test_streams.py The purpose of this program is to be the integration test for the streams list command for Vyrtuous.

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

GUILD_SNOWFLAKE = 10000000000000500
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, target",
    [
        ("!streams", "all"),
        ("!streams", "{channel_snowflake}"),
        ("!streams", "<#{channel_snowflake}>"),
        ("!streams", "{guild_snowflake}"),
    ],
)
async def test_streams(bot, command: Optional[str], target):
    """
    List members who are registered in the PostgresSQL database
    'vyrtuous' in the table 'streaming'.

    Parameters
    ----------
    all : str, optional
        Generic showing all streams in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with streams
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where streams are present.

    Examples
    --------
    >>> !streams "all"
    [{emoji} Streaming Channels\n Guild1\n Guild2]

    >>> !streams 10000000000000500
    [{emoji} Streaming Channels\n Guild1]

    >>> !streams <@10000000000000010>
    [{emoji} Streaming Channels for Channel1]

    >>> !streams 10000000000000010
    [{emoji} Streaming Channels for Channel1]
    """
    t = target.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE, guild_snowflake=GUILD_SNOWFLAKE
    )
    full = f"{command} {t}"
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
    command = await admin_commands.list_streaming_text_command(ctx, target=t)
