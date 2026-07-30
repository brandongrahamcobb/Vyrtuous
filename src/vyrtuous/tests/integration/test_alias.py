"""!/bin/python3
test_cogs.py The purpose of this program is to be the integration test for the cogs list command for Vyrtuous.

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

from vyrtuous.db.alias import Alias
from vyrtuous.tests.integration.test_suite import check_permissions, send_message

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
ROLE_SNOWFLAKE = 10000000000000200
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
COMMAND = "alias"
BASE_PERMISSIONS = ["command.alias.create", "command.alias.scope.channel"]
TABLE_NAME = Alias.__tablename__


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "category, alias_name, channel, role, extra_permissions",
    [
        ("ban", "testban", "{channel_snowflake}", None, []),
        ("flag", "testflag", "{channel_snowflake}", None, []),
        ("tmute", "testtmute", "{channel_snowflake}", None, []),
        ("vmute", "testmute", "{channel_snowflake}", None, []),
        ("ban", "testban", None, None, []),
        ("flag", "testflag", None, None, []),
        (
            "role",
            "testrole",
            "{channel_snowflake}",
            "{role_snowflake}",
            ["command.alias.scope.role"],
        ),
        ("tmute", "testtmute", None, None, []),
        ("vmute", "testmute", None, None, []),
    ],
)
async def test_alias_text_command(
    bot,
    prefix: str,
    category: str,
    alias_name: str,
    channel: str | None,
    role: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Create command aliases in the PostgreSQL
    database 'vyrtuous' in the table 'command_aliases'.

    Parameters
    ----------
    category : str
        Resolves to: CategoryObject
        Examples: ban, flag, role, tmute or vmute

    alias_name : str
        Examples: testban, testmute

    channel (Optional) : str | int
        Resolves to: discord.VoiceChannel | discord.StageChannel | discord.TextChannel
        Examples: 10000000000000010 | <#10000000000000010>

    role (Optional) : str | int
        Resolves to: discord.Role
        Examples: 10000000000000010 | <@&10000000000000010>

    Examples
    --------
    >>> !alias ban testban
    Embed
    """
    assert COMMAND in docstring
    assert TABLE_NAME in docstring
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
            full = f"{prefix}{COMMAND} {category} {alias_name}"
            if channel is None:
                c = channel
            else:
                c = channel.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                )
                full += f" {c}"
            if role is None:
                r = role
            else:
                r = role.format(
                    role_snowflake=ROLE_SNOWFLAKE,
                )
                full += f" {r}"
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]


COLUMNS = [
    ("channel_snowflake", "bigint", False),
    ("guild_snowflake", "bigint", False),
    ("role_snowflake", "bigint", True),
    ("category", "text", False),
    ("alias_name", "text", False),
    ("created_at", "timestamp with time zone", True),
    ("updated_at", "timestamp with time zone", True),
]


@pytest.mark.asyncio
@pytest.mark.parametrize("field, datatype, nullable", COLUMNS)
async def test_alias_database_table(bot, field: str, datatype: str, nullable: bool):
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
