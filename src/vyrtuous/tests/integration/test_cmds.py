"""test_cmds.py The purpose of this program is to be the integration test for the cmds list command for Vyrtuous.

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
        ("!cmds", "all"),
        ("!cmds", "{channel_snowflake}"),
        ("!cmds", "<#{channel_snowflake}>"),
        ("!cmds", "{guild_snowflake}"),
    ],
)
async def test_cmds(bot, command: Optional[str], target):
    """
    List channels which are registered in the PostgresSQL database
    'vyrtuous' in the table 'command_aliases'.

    Parameters
    ----------
    all : str, optional
        Generic showing all command aliases in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with command aliases
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where command aliases are present.

    Examples
    --------
    >>> !cmds "all"
    [{emoji} Aliases\n Guild1\n Guild2]

    >>> !cmds 10000000000000500
    [{emoji} Aliases\n Guild1]

    >>> !cmds <@10000000000000010>
    [{emoji} Aliases for Channel1]

    >>> !cmds 10000000000000010
    [{emoji} Aliases for Channel1]
    """
    t = target.format(
        channel_snowflake=TEXT_CHANNEL_SNOWFLAKE, guild_snowflake=GUILD_SNOWFLAKE
    )
    full = f"{command} {t}"
    captured = await send_message(bot=bot, content=full)
    assert captured.content
    objects = setup(bot)
    msg = build_message(
        author=objects.get("author", None), channel=objects.get("channel", None), content=full, guild=objects.get("guild", None), state=objects.get("state", None)
    )
    ctx = context(bot=bot, message=msg, prefix="!")
    mod_commands = bot.get_cog("ModeratorCommands")
    command = await mod_commands.list_commands_text_command(ctx, target=t)
