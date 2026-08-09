"""!/bin/python3
test_tmutes.py The purpose of this program is to be a shared test for the tmutes app and text commands.

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
from datetime import datetime, timezone
from unittest.mock import patch

import pytest

from vyrtuous.cache.registry import MemberState
from vyrtuous.db.text_mute import TextMute
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
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005
COMMAND = "tmutes"
BASE_PERMISSIONS = [
    "command.info.text-mutes",
]
TABLE_NAME = TextMute.__tablename__


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target, other_guild, extra_permissions",
    [
        ("{channel_snowflake}", None, ["command.info.scope.channel"]),
        ("<#{channel_snowflake}>", None, ["command.info.scope.channel"]),
        (
            "{channel_snowflake}",
            "{other_guild_snowflake}",
            ["command.info.scope.channel", "other_guilds"],
        ),
        (
            "<#{channel_snowflake}>",
            "{other_guild_snowflake}",
            ["command.info.scope.channel", "other_guilds"],
        ),
        ("{guild_snowflake}", None, ["command.info.scope.guild", "other_channels"]),
        ("{member_snowflake}", None, ["command.info.scope.member", "other_channels"]),
        (
            "<@{member_snowflake}>",
            None,
            ["command.info.scope.member", "other_channels"],
        ),
        (
            "{member_snowflake}",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_channels", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_channels", "other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            None,
            ["command.info.scope.member", "other_channels"],
        ),
        (
            "{simplified_member_snowflake}",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_channels", "other_guilds"],
        ),
    ],
)
async def test_tmutes_text_command(
    bot,
    prefix: str,
    target: str | None,
    other_guild: str | None,
    extra_permissions: list[str],
):
    docstring = """
    List text-mutes on members which are registered in the PostgresSQL database
    'vyrtuous' in the table 'active_text_mutes'.

    Parameters
    ----------
    target (Optional) : str | int
        Defaults to context channel
        Resolves to: int | discord.VoiceChannel | discord.TextChannel | discord.StageChannel | discord.Member | discord.Guild
        Examples: 10000000000000010 | <#10000000000000010> | <@10000000000000010>

    guild (Optional) : str | int
        Resolves to: discord.Guild
        Examples: 10000000000000010

    Example
    --------
    >>> !tmutes
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
            if target is None:
                t = target
            else:
                bot.registry.get(MemberState).active.update(
                    {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
                )
                t = target.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    guild_snowflake=GUILD_SNOWFLAKE,
                    member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                    simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
                )
            if other_guild is None:
                g = other_guild
                full = f"{prefix}{COMMAND} {t}"
            else:
                g = other_guild.format(other_guild_snowflake=OTHER_GUILD_SNOWFLAKE)
                full = f"{prefix}{COMMAND} {t} {g}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target, other_guild, extra_permissions",
    [
        ("{channel_snowflake}", None, ["command.info.scope.channel"]),
        ("<#{channel_snowflake}>", None, ["command.info.scope.channel"]),
        (
            "{channel_snowflake}",
            "{other_guild_snowflake}",
            ["command.info.scope.channel", "other_guilds"],
        ),
        (
            "<#{channel_snowflake}>",
            "{other_guild_snowflake}",
            ["command.info.scope.channel", "other_guilds"],
        ),
        ("{guild_snowflake}", None, ["command.info.scope.guild", "other_channels"]),
        ("{member_snowflake}", None, ["command.info.scope.member", "other_channels"]),
        (
            "<@{member_snowflake}>",
            None,
            ["command.info.scope.member", "other_channels"],
        ),
        (
            "{member_snowflake}",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_channels", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_channels", "other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            None,
            ["command.info.scope.member", "other_channels"],
        ),
        (
            "{simplified_member_snowflake}",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_channels", "other_guilds"],
        ),
    ],
)
async def test_tmutes_app_command(
    bot,
    target: str | None,
    other_guild: str | None,
    extra_permissions: list[str],
):
    docstring = """
    List text-mutes on members which are registered in the PostgresSQL database
    'vyrtuous' in the table 'active_text_mutes'.

    Parameters
    ----------
    target (Optional) : str | int
        Resolves to: int | discord.VoiceChannel | discord.TextChannel | discord.StageChannel | discord.Member | discord.Guild
        Examples: 10000000000000010 | <#10000000000000010> | <@10000000000000010>

    guild (Optional) : str | int
        Resolves to: discord.Guild
        Examples: 10000000000000010

    Example
    --------
    >>> /tmutes
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
            command = cog.list_text_mutes_app_command
            if target is None:
                t = target
            else:
                bot.registry.get(MemberState).active.update(
                    {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
                )
                t = target.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    guild_snowflake=GUILD_SNOWFLAKE,
                    member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                    simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
                )
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
            if t:
                resolved_target = await transformer.transform(inx, t)
            else:
                resolved_target = None
            if g:
                resolved_guild = await transformer.transform(inx, g)
            else:
                resolved_guild = None
            async with capture_command() as end_results:
                await command.callback(
                    cog, interaction=inx, target=resolved_target, guild=resolved_guild
                )
        for kind, content in end_results:
            assert kind == "success"


COLUMNS = [
    ("channel_snowflake", "bigint", False),
    ("guild_snowflake", "bigint", False),
    ("member_snowflake", "bigint", False),
    ("expires_in", "timestamp with time zone", True),
    ("created_at", "timestamp with time zone", True),
    ("updated_at", "timestamp with time zone", True),
    ("last_muted", "timestamp with time zone", True),
    ("reset", "boolean", True),
    ("reason", "text", False),
]


@pytest.mark.asyncio
@pytest.mark.parametrize("field, datatype, nullable", COLUMNS)
async def test_active_text_mutes_database_table(
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
