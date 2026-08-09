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
from unittest.mock import patch

import pytest

from vyrtuous.db.alias import Alias
from vyrtuous.tests.integration.test_suite import check_permissions, send_message

GUILD_SNOWFLAKE = 10000000000000500
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011

VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
COMMAND = "aliases"
BASE_PERMISSIONS = ["command.info.aliases"]
TABLE_NAME = Alias.__tablename__


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target, extra_permissions",
    [
        ("{channel_snowflake}", ["command.info.scope.channel"]),
        ("<#{channel_snowflake}>", ["command.info.scope.channel"]),
        ("{guild_snowflake}", ["command.info.scope.guild", "other_channels"]),
    ],
)
async def test_aliases_text_command(
    bot, prefix: str, target: str | None, extra_permissions: list[str]
):
    docstring = """
    List aliases which are registered in the PostgreSQL database
    'vyrtuous' in the table 'command_aliases'.

    Parameters
    ----------
    target (Optional) : str | int
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel |  discord.Guild
        Examples: 10000000000000010 | <#10000000000000010>

    Example
    --------
    >>> !aliases
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
            if target is None:
                t = target
            else:
                t = target.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    guild_snowflake=GUILD_SNOWFLAKE,
                )
                full += f" {t}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


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
