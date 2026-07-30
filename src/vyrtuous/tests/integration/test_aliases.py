"""!/bin/python3
test_cmds.py The purpose of this program is to be the integration test for the cmds list command for Vyrtuous.

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
from unittest.mock import AsyncMock, patch

import pytest

from vyrtuous.cache.registry import PermissionState
from vyrtuous.tests.integration.test_suite import send_message

GUILD_SNOWFLAKE = 10000000000000500
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011

VOICE_CHANNEL_SNOWFLAKE = 10000000000000011


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, target",
    [
        ("Moderator", "aliases", "{channel_snowflake}"),
        ("Moderator", "aliases", "<#{channel_snowflake}>"),
        ("Administrator", "aliases", "{guild_snowflake}"),
    ],
)
async def test_cmds(bot, command: str, prefix: str, target, permission_role):
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
    permission_state = bot.registry.get(PermissionState)
    t = target.format(
        channel_snowflake=VOICE_CHANNEL_SNOWFLAKE, guild_snowflake=GUILD_SNOWFLAKE
    )
    full = f"{prefix}{command} {t}"
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.utils.permissions.permission_service.resolve_effective_group",
                    new=AsyncMock(
                        return_value=permission_state.groups.get(
                            permission_role.lower()
                        )
                    ),
                )
            )
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]
