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
from unittest.mock import patch

import pytest

from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    send_message,
    setup,
)

ROLE_ID = 10000000000000200
OTHER_GUILD_SNOWFLAKE = 10000000000000501


COMMAND = "members"
BASE_PERMISSIONS = ["command.info.members", "command.info.scope.role"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "role, other_guild, extra_permissions",
    [
        ("{role_snowflake}", None, []),
        ("{role_snowflake}", "{other_guild_snowflake}", ["other_guilds"]),
    ],
)
@pytest.mark.asyncio
async def test_members_text_command(
    bot,
    prefix: str,
    role: str,
    other_guild: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Get the role members by ID.

    Parameters
    ----------
    role : str | int
        Resolves to: TargetObject
        Examples: 10000000000000010 | <@&10000000000000010>

    Examples
    --------
    >>> !members 10000000000000010
    Embed
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
            r = role.format(role_snowflake=ROLE_ID)
            full = f"{prefix}{COMMAND} {r}"
            if other_guild is None:
                g = other_guild
            else:
                g = other_guild.format(other_guild_snowflake=OTHER_GUILD_SNOWFLAKE)
                full += f" {g}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "role, other_guild, extra_permissions",
    [
        ("{role_snowflake}", None, []),
        ("{role_snowflake}", "{other_guild_snowflake}", ["other_guilds"]),
    ],
)
async def test_members_app_command(
    bot, role: str, other_guild: str | None, extra_permissions: list[str]
):
    docstring = """
    Get the role members by ID.

    Parameters
    ----------
    role : str | int
        Resolves to: TargetObject
        Examples: 10000000000000010 | <@&10000000000000010>

    Examples
    --------
    >>> /members 10000000000000010
    Embed
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
            command = cog.list_role_members_app_command
            r = role.format(role_snowflake=ROLE_ID)
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
            resolved_role = await transformer.transform(inx, r)
            if g:
                resolved_guild = await transformer.transform(inx, g)
            else:
                resolved_guild = None

            async with capture_command() as end_results:
                await command.callback(
                    cog, interaction=inx, role=resolved_role, guild=resolved_guild
                )
            for kind, content in end_results:
                assert kind == "success"
