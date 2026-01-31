"""test_mutes.py The purpose of this program is to be the integration test for the mutes list command for Vyrtuous.

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
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, target",
    [
        ("!smutes", "all"),
        ("!smutes", "{guild_snowflake}"),
        ("!smutes", "{member_snowflake}"),
        ("!smutes", "<@{member_snowflake}>"),
    ],
)
async def test_server_mutes(bot, command: str, target):
    """
    List voice-mutes on members which are registered in the PostgresSQL database
    'vyrtuous' in the table 'active_voice_mutes'.

    Parameters
    ----------
    all : str, optional
        Generic showing all voice mutes in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with voice-mutes on members
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where mutes are present.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who has been voice-muted
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !smutes "all"
    [{emoji} Voice Mutes\n Guild1\n Guild2]

    >>> !smutes <#10000000000000010>
    [{emoji} Voice Mutes for Channel1\n Member1\n Member2]

    >>> !smutes 10000000000000010
    [{emoji} Voice Mutes for Channel1\n Member1\n Member2]

    >>> !smutes 10000000000000500
    [{emoji} Voice Mutes\n Guild1]

    >>> !smutes <@10000000000000003>
    [{emoji} Voice Mutes for Member1\n Guild1\n Guild2]

    >>> !smutes 10000000000000003
    [{emoji} Voice Mutes for Member1\n Guild1\n Guild2]
    """
    t = target.format(
        guild_snowflake=GUILD_SNOWFLAKE,
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
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
    command = await admin_commands.list_server_mutes_text_command(ctx, target=t)
