"""test_tmutes.py The purpose of this program is to be the integration test for the tmutes list command for Vyrtuous.

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

from vyrtuous.tests.integration.test_suite import send_message

GUILD_SNOWFLAKE = 10000000000000500
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
TEXT_CHANNEL_SNOWFLAKE = 10000000000000010


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("!tmutes all"),
        ("!tmutes {channel_snowflake}"),
        ("!tmutes <#{channel_snowflake}>"),
        ("!tmutes {guild_snowflake}"),
        ("!tmutes {member_snowflake}"),
        ("!tmutes <@{member_snowflake}>"),
    ],
)
async def test_text_mutes(bot, command: Optional[str]):
    """
    List text-mutes on members which are registered in the PostgresSQL database
    'vyrtuous' in the table 'active_text_mutes'.

    Parameters
    ----------
    all : str, optional
        Generic showing all text mutes in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with text-mutes on members
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where text mutes are present.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who has been text-muted
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !tmutes "all"
    [{emoji} Text-Mutes\n Guild1\n Guild2]

    >>> !tmutes <#10000000000000010>
    [{emoji} Text-Mutes for Channel1\n Member1\n Member2]

    >>> !tmutes 10000000000000010
    [{emoji} Text-Mutes for Channel1\n Member1\n Member2]

    >>> !tmutes 10000000000000500
    [{emoji} Text-Mutes\n Guild1]

    >>> !tmutes <@10000000000000002>
    [{emoji} Text-Mutes for Member1\n Guild1\n Guild2]

    >>> !tmutes 10000000000000002
    [{emoji} Text-Mutes for Member1\n Guild1\n Guild2]
    """
    formatted = command.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE,
        guild_snowflake=GUILD_SNOWFLAKE,
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE
    )
    captured = await send_message(bot=bot, content=formatted)
    assert captured.content
