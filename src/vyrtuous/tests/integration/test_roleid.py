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
from contextlib import ExitStack
from unittest.mock import AsyncMock, patch

import pytest

from vyrtuous.cache.registry import PermissionState
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    send_message,
    setup,
)

ROLE_NAME = "Vegan"
ROLE_ID = 10000000000000200
OTHER_GUILD_SNOWFLAKE = 10000000000000501


COMMAND = "roleid"
BASE_PERMISSIONS = ["command.info.roleid", "command.info.scope.role"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "other_guild, extra_permissions",
    [
        (None, []),
        ("{other_guild_snowflake}", ["other_guilds"]),
    ],
)
@pytest.mark.asyncio
async def test_roleid_text_command(
    bot,
    prefix: str,
    other_guild: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Get the role ID by name.

    Parameters
    ----------
    None

    Examples
    --------
    >>> !roleid
    {emoji} Role `Vegan` has ID `10000000000000200`
    """
    assert COMMAND in docstring
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        extra_permissions.extend(BASE_PERMISSIONS)
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions",
                    side_effect=check_permissions(BASE_PERMISSIONS),
                )
            )
            full = f"{prefix}{COMMAND} {ROLE_NAME}"
            if other_guild is None:
                g = other_guild
            else:
                g = other_guild.format(other_guild_snowflake=OTHER_GUILD_SNOWFLAKE)
                full += f" {g}"
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "other_guild, extra_permissions",
    [
        (None, []),
        ("{other_guild_snowflake}", ["other_guilds"]),
    ],
)
async def test_roleid_app_command(
    bot, other_guild: str | None, extra_permissions: list[str]
):
    docstring = """
    Get the role ID by name.

    Parameters
    ----------
    None

    Examples
    --------
    >>> /roleid
    {emoji} Role `Vegan` has ID `10000000000000200`
    """
    assert COMMAND in docstring
    if (
        os.environ["TEST_MODE"].lower() == "app"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        extra_permissions.extend(BASE_PERMISSIONS)
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions",
                    side_effect=check_permissions(BASE_PERMISSIONS),
                )
            )
            cog = bot.get_cog("InfoAppCommands")
            command = cog.get_role_snowflake_app_command
            if other_guild is None:
                g = other_guild
            else:
                g = other_guild.format(other_guild_snowflake=OTHER_GUILD_SNOWFLAKE)
            objects = setup(bot)
            msg = build_message(
                author=objects.get("author", None),
                channel=objects.get("text_channel", None),
                content="",
                guild=objects.get("guild", None),
                state=objects.get("state", None),
            )
            inx = interaction(
                bot=bot,
                channel=objects.get("text_channel", None),
                guild=objects.get("guild", None),
                message=msg,
            )
            transformer = AppTarget()
            if g:
                resolved_guild = await transformer.transform(inx, g)
            else:
                resolved_guild = None

            async with capture_command() as end_results:
                await command.callback(
                    cog, interaction=inx, role_name=ROLE_NAME, guild=resolved_guild
                )
            for kind, content in end_results:
                assert kind == "success"
