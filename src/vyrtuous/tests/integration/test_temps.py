"""test_temps.py The purpose of this program is to be the integration test for the temps list command for Vyrtuous.

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
        ("!temps", "all"),
        ("!temps", "{channel_snowflake}"),
        ("!temps", "<#{channel_snowflake}>"),
        ("!temps", "{guild_snowflake}"),
    ],
)
async def test_temps(bot, command: Optional[str], target):
    """
    List channels which are registered in the PostgresSQL database
    'vyrtuous' in the table 'temporary_rooms'.

    Parameters
    ----------
    all : str, optional
        Generic showing all temporary rooms in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with temps
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where temporary rooms are present.

    Examples
    --------
    >>> !temps "all"
    [{emoji} Temporary Rooms\n Guild1\n Guild2]

    >>> !temps 10000000000000500
    [{emoji} Temporary Rooms\n Guild1]

    >>> !temps <@10000000000000010>
    [{emoji} Temporary Rooms for Channel1]

    >>> !temps 10000000000000010
    [{emoji} Temporary Rooms for Channel1]
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
    command = await admin_commands.list_temp_rooms_text_command(ctx, target=t)
