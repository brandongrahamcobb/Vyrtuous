"""!/bin/python3
test_bans.py The purpose of this program is to be the integration test for the bans list command for Vyrtuous.

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
from vyrtuous.db.ban import Ban
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_flags import BASE_PERMISSIONS
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
BASE_PERMISSIONS = ["command.info.bans"]
COMMAND = "blacklists"
TABLE_NAME = Ban.__tablename__


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
        ("{guild_snowflake}", None, ["command.info.scope.guild"]),
        ("{member_snowflake}", None, ["command.info.scope.member"]),
        ("<@{member_snowflake}>", None, ["command.info.scope.member"]),
        (
            "{member_snowflake}",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_guilds"],
        ),
        ("{simplified_member_snowflake}", None, ["command.info.scope.member"]),
        (
            "{simplified_member_snowflake}>",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_guilds"],
        ),
    ],
)
async def test_blacklists_text_command(
    bot,
    prefix: str,
    target: str | None,
    other_guild: str | None,
    extra_permissions: list[str],
):
    docstring = """
    List blacklisted members which are registered in the PostgresSQL database
    'vyrtuous' in the table 'active_bans'.

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
    >>> !blacklists
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
            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
            )
            full = f"{prefix}{COMMAND}"
            if target is None:
                t = target
            else:
                t = target.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    guild_snowflake=GUILD_SNOWFLAKE,
                    member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                    simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
                )
                full += f" {t}"
            if other_guild is None:
                g = None
            else:
                g = other_guild.format(
                    other_guild_snowflake=OTHER_GUILD_SNOWFLAKE,
                )
                full += f" {g}"
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
            "<#{channel_snowflake}",
            "{other_guild_snowflake}",
            ["command.info.scope.channel", "other_guilds"],
        ),
        ("{guild_snowflake}", None, ["command.info.scope.guild"]),
        ("{member_snowflake}", None, ["command.info.scope.member"]),
        ("<@{member_snowflake}>", None, ["command.info.scope.member"]),
        (
            "{member_snowflake}",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_guilds"],
        ),
        (
            "<@{member_snowflake}",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_guilds"],
        ),
        ("{simplified_member_snowflake}", None, ["command.info.scope.member"]),
        (
            "{simplified_member_snowflake}>",
            "{other_guild_snowflake}",
            ["command.info.scope.member", "other_guilds"],
        ),
    ],
)
async def test_blacklists_app_command(
    bot,
    target: str | None,
    other_guild: str | None,
    extra_permissions: list[str],
):
    docstring = """
    List blacklisted members which are registered in the PostgresSQL database
    'vyrtuous' in the table 'active_bans'.

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
    >>> !blacklists
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
            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
            )
            cog = bot.get_cog("InfoAppCommands")
            command = cog.list_blacklists_app_command
            if target is None:
                t = target
            else:
                t = target.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    guild_snowflake=GUILD_SNOWFLAKE,
                    member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                    simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
                )
            if other_guild is None:
                g = None
            else:
                g = other_guild.format(
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
                if t:
                    resolved_target = await transformer.transform(inx, t)
                else:
                    resolved_target = None
                if g:
                    resolved_guild = await transformer.transform(inx, g)
                else:
                    resolved_guild = None
                await command.callback(
                    cog, interaction=inx, target=resolved_target, guild=resolved_guild
                )
            for kind, content in end_results:
                assert kind == "success"
