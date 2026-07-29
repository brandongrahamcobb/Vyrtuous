"""!/bin/python3
test_roleid.py The purpose of this program is to be the integration test for the roleid list command for Vyrtuous.

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

ROLE_NAME = "Vegan"
GUILD_SNOWFLAKE = 10000000000000500


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, role, guild",
    [
        ("Administrator", "roleid", "{role_name}", None),
        ("Administrator", "roleid", "{role_name}", "{guild_snowflake}"),
    ],
)
async def test_roleid(bot, command: str, prefix: str, role, guild, permission_role):
    """
    Fetch a role snowflake in a guild

    Parameters
    ----------
    role_snowflake
        The snowflake or mention of a role

    Examples
    --------
    >>> !roleid Vegan
    [{emoji} Role `Vegan` has the id `10000000000000200`]

    """
    r = role.format(role_name=ROLE_NAME)
    g = None
    full = f"{prefix}{command} {r}"
    if guild:
        g = guild.format(guild_snowflake=GUILD_SNOWFLAKE)
        full = f"{prefix}{command} {r} {g}"
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        captured = await send_message(bot=bot, content=full)
        assert captured == ["success"]
    elif (
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
            command = cog.get_role_id_app_command
            transformer = AppTarget()
            if g:
                resolved_guild = await transformer.transform(inx, g)
            else:
                resolved_guild = None
            await command.callback(
                cog, interaction=inx, role_name=r, guild=resolved_guild
            )
        for kind, content in end_results:
            assert kind == "success"
