"""!/bin/python3
test_aroles.py The purpose of this program is to be the integration test for the aroles list command for Vyrtuous.

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

from vyrtuous.db.autoassign import AutoAssignRole
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    send_message,
    setup,
)

GUILD_SNOWFLAKE = 10000000000000500
OTHER_GUILD_SNOWFLAKE = 10000000000000501


COMMAND = "autoassigns"
BASE_PERMISSIONS = [
    "command.info.autoassigns",
    "command.info.scope.guild",
    "other_channels",
]
TABLE_NAME = AutoAssignRole.__tablename__


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "guild, extra_permissions",
    [
        (None, []),
        (None, []),
        ("{guild_snowflake}", []),
        ("{other_guild_snowflake}", ["other_guilds"]),
    ],
)
async def test_autoassigns_text_command(
    bot, prefix: str, guild: str | None, extra_permissions: list[str]
):
    docstring = """
    List autoassigns which are registered in the PostgreSQL database
    'vyrtuous' in the table 'autoassign_roles'.

    Parameters
    ----------
    guild (Optional) : int
        Resolves to: discord.Guild
        Examples: 10000000000000010 

    Example
    --------
    >>> !autoassigns
    Embed
    """
    assert TABLE_NAME in docstring
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
                    side_effect=check_permissions(extra_permissions),
                )
            )
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions_at_all",
                    side_effect=check_permissions(extra_permissions),
                )
            )

            full = f"{prefix}{COMMAND}"
            if guild is None:
                g = guild
            else:
                g = guild.format(
                    guild_snowflake=GUILD_SNOWFLAKE,
                    other_guild_snowflake=OTHER_GUILD_SNOWFLAKE,
                )
                full += f" {g}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "guild, extra_permissions",
    [
        (None, []),
        (None, []),
        ("{guild_snowflake}", []),
        ("{other_guild_snowflake}", ["other_guilds"]),
    ],
)
async def test_autoassigns_app_command(
    bot, prefix: str, guild: str | None, extra_permissions: list[str]
):
    docstring = """
    List autoassigns which are registered in the PostgreSQL database
    'vyrtuous' in the table 'autoassign_roles'.

    Parameters
    ----------
    guild (Optional) : int
        Resolves to: discord.Guild
        Examples: 10000000000000010 

    Example
    --------
    >>> !autoassigns
    Embed
    """
    assert TABLE_NAME in docstring
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
                    side_effect=check_permissions(extra_permissions),
                )
            )
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions_at_all",
                    side_effect=check_permissions(extra_permissions),
                )
            )

            cog = bot.get_cog("InfoAppCommands")
            command = cog.list_autoassignment_roles_app_command
            if guild is None:
                g = guild
            else:
                g = guild.format(
                    guild_snowflake=GUILD_SNOWFLAKE,
                    other_guild_snowflake=OTHER_GUILD_SNOWFLAKE,
                )

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
            async with capture_command() as end_results:
                transformer = AppTarget()
                if g:
                    resolved_guild = await transformer.transform(inx, g)
                else:
                    resolved_guild = None
                await command.callback(cog, interaction=inx, guild=resolved_guild)
            for kind, content in end_results:
                assert kind == "success"


COLUMNS = [
    ("group_alias", "text", False),
    ("channel_snowflake", "bigint", True),
    ("guild_snowflake", "bigint", False),
    ("role_snowflake", "bigint", False),
    ("created_at", "timestamp with time zone", True),
    ("updated_at", "timestamp with time zone", True),
    ("id", "integer", False)
]


@pytest.mark.asyncio
@pytest.mark.parametrize("field, datatype, nullable", COLUMNS)
async def test_autoassign_channels_database_table(
    bot, field: str, datatype: str, nullable: bool
):
    async with bot.db_pool.acquire() as conn:
        statement = await conn.prepare(f"SELECT * FROM {TABLE_NAME}")
        columns = statement.get_attributes()
        assert len(columns) == len(COLUMNS)
        row = await conn.fetchrow(
            f"""
            SELECT
                column_name,
                data_type,
                is_nullable
            FROM information_schema.columns
            WHERE table_schema = 'public'
              AND table_name = $1
              AND column_name = $2
            ORDER BY ordinal_position
        """,
            TABLE_NAME,
            field,
        )
    assert row is not None
    assert row["column_name"] == field
    assert row["data_type"] == datatype
    assert row["is_nullable"] == ("YES" if nullable else "NO")
