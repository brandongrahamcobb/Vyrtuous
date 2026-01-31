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

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, member",
    [
        ("!summary", "{member_snowflake}"),
        ("!summary", "<@{member_snowflake}>"),
    ],
)
async def test_summary(bot, command: str, member):
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
    >>> !summary <@10000000000000003>
    [{emoji} Infractions for Member1\n Guild1\n Guild2]

    >>> !summary 10000000000000003
    [{emoji} Infractions for Member1\n Guild1\n Guild2]
    """
    m = member.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
    )
    full = f"{command} {m}"
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
    mod_commands = bot.get_cog("ModeratorTextCommands")
    command = await mod_commands.list_moderation_summary_text_command(ctx, member=m)
