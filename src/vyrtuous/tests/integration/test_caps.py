"""!/bin/python3
test_caps.py The purpose of this program is to be the integration test for the caps list command for Vyrtuous.

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
    "permission_role, command, target, guild",
    [
        ("Administrator", "caps", "{channel_snowflake}", None),
        ("Administrator", "caps", "<#{channel_snowflake}>", "{guild_snowflake}"),
        ("Administrator", "caps", "{guild_snowflake}", None),
    ],
)
async def test_caps(bot, command: str, prefix: str, target, guild, permission_role):
    """
    List caps in the PostgresSQL database
    'vyrtuous' in the table 'active_caps'.

    Parameters
    ----------
    all : str, optional
        Generic showing all caps in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with caps
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where caps are present.

    Examples
    --------
    >>> !caps "all"
    [{emoji} Caps\n Guild1\n Guild2]

    >>> !caps 10000000000000500
    [{emoji} Caps\n Guild1]

    >>> !caps <@10000000000000010>
    [{emoji} Caps for Channel1]

    >>> !caps 10000000000000010
    [{emoji} Caps for Channel1]
    """
    t = target.format(
        channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
        guild_snowflake=GUILD_SNOWFLAKE,
    )
    g = None
    if guild is None:
        full = f"{prefix}{command} {t}"
    else:
        g = guild.format(
            channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
            guild_snowflake=GUILD_SNOWFLAKE,
        )
        full = f"{prefix}{command} {t} {g}"
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
            command = cog.list_caps_app_command
            transformer = AppTarget()
            if t:
                resolved = await transformer.transform(inx, t)
            else:
                resolved = None
            if g:
                resolved_guild = await transformer.transform(inx, g)
            else:
                resolved_guild = None
            await command.callback(
                cog, interaction=inx, target=resolved, guild=resolved_guild
            )
        for kind, content in end_results:
            assert kind == "success"
