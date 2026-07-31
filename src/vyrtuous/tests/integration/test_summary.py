"""!/bin/python3
test_mutes.py The purpose of this program is to be the integration test for the mutes list command for Vyrtuous.

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
from vyrtuous.models.scope import AppScope
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_rmv import BASE_PERMISSIONS
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    send_message,
    setup,
)

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
OTHER_GUILD_SNOWFLAKE = 10000000000000501
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005

COMMAND = "summary"
BASE_PERMISSIONS = [
    "command.info.bans",
    "command.info.flags",
    "command.info.text-mutes",
    "command.info.voice-mutes",
]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "member, scope, other_guild, extra_permissions",
    [
        ("{member_snowflake}", None, None, []),
        (
            "{member_snowflake}",
            "all",
            None,
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
            ],
        ),
        ("{member_snowflake}", "auto", None, ["command.info.voice-mutes.auto"]),
        ("{member_snowflake}", "click", None, ["command.info.voice-mutes.click"]),
        ("{member_snowflake}", "command", None, ["command.info.voice-mutesi.command"]),
        ("{member_snowflake}", "server", None, ["command.info.voice-mutes.server"]),
        (
            "<@{member_snowflake}>",
            "all",
            None,
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
            ],
        ),
        ("<@{member_snowflake}>", "auto", None, ["command.info.voice-mutes.auto"]),
        ("<@{member_snowflake}>", "click", None, ["command.info.voice-mutes.click"]),
        (
            "<@{member_snowflake}>",
            "command",
            None,
            ["command.info.voice-mutes.command"],
        ),
        ("<@{member_snowflake}>", "server", None, ["command.info.voice-mutes.server"]),
        (
            "{member_snowflake}",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{member_snowflake}",
            "auto",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.auto", "other_guilds"],
        ),
        (
            "{member_snowflake}",
            "click",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.click", "other_guilds"],
        ),
        (
            "{member_snowflake}",
            "command",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.command", "other_guilds"],
        ),
        (
            "{member_snowflake}",
            "server",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.server", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "auto",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.auto", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "click",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.click", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "command",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.command", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "server",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.server", "other_guilds"],
        ),
        ("{simplified_member_snowflake}", None, None, []),
        (
            "{simplified_member_snowflake}",
            "all",
            None,
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "auto",
            None,
            ["command.info.voice-mutes.auto"],
        ),
        (
            "{simplified_member_snowflake}",
            "click",
            None,
            ["command.info.voice-mutes.click"],
        ),
        (
            "{simplified_member_snowflake}",
            "command",
            None,
            ["command.info.voice-mutes.command"],
        ),
        (
            "{simplified_member_snowflake}",
            "server",
            None,
            ["command.info.voice-mutes.server"],
        ),
        (
            "{simplified_member_snowflake}",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "auto",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.auto", "other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "click",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.click", "other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "command",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.command", "other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "server",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.server", "other_guilds"],
        ),
    ],
)
async def test_summary_text_command(
    bot,
    prefix: str,
    member: str,
    scope: str | None,
    other_guild: str | None,
    extra_permissions: list[str],
):
    docstring = """
    List bans, flags, text-mutes and voice-mutes on a member

    Parameters
    ----------
    member : str | int
        Resolves to: int | discord.Member
        Examples: 10000000000000010 | <@10000000000000010>

    guild (Optional) : str | int
        Resolves to: discord.Guild
        Examples: 10000000000000010

    Examples
    --------
    >>> !summary 10000000000000003
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
            m = member.format(
                member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
            )
            full = f"{prefix}{COMMAND} {m}"
            if scope is None:
                s = scope
            else:
                s = scope
                full += f" {s}"
            if other_guild is None:
                g = other_guild
            else:
                g = other_guild.format(other_guild_snowflake=OTHER_GUILD_SNOWFLAKE)
                full += f" {g}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "member, scope, other_guild, extra_permissions",
    [
        ("{member_snowflake}", None, None, []),
        (
            "{member_snowflake}",
            "all",
            None,
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
            ],
        ),
        ("{member_snowflake}", "auto", None, ["command.info.voice-mutes.auto"]),
        ("{member_snowflake}", "click", None, ["command.info.voice-mutes.click"]),
        ("{member_snowflake}", "command", None, ["command.info.voice-mutesi.command"]),
        ("{member_snowflake}", "server", None, ["command.info.voice-mutes.server"]),
        (
            "<@{member_snowflake}>",
            "all",
            None,
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
            ],
        ),
        ("<@{member_snowflake}>", "auto", None, ["command.info.voice-mutes.auto"]),
        ("<@{member_snowflake}>", "click", None, ["command.info.voice-mutes.click"]),
        (
            "<@{member_snowflake}>",
            "command",
            None,
            ["command.info.voice-mutes.command"],
        ),
        ("<@{member_snowflake}>", "server", None, ["command.info.voice-mutes.server"]),
        (
            "{member_snowflake}",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{member_snowflake}",
            "auto",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.auto", "other_guilds"],
        ),
        (
            "{member_snowflake}",
            "click",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.click", "other_guilds"],
        ),
        (
            "{member_snowflake}",
            "command",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.command", "other_guilds"],
        ),
        (
            "{member_snowflake}",
            "server",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.server", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "auto",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.auto", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "click",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.click", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "command",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.command", "other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "server",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.server", "other_guilds"],
        ),
        ("{simplified_member_snowflake}", None, None, []),
        (
            "{simplified_member_snowflake}",
            "all",
            None,
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "auto",
            None,
            ["command.info.voice-mutes.auto"],
        ),
        (
            "{simplified_member_snowflake}",
            "click",
            None,
            ["command.info.voice-mutes.click"],
        ),
        (
            "{simplified_member_snowflake}",
            "command",
            None,
            ["command.info.voice-mutes.command"],
        ),
        (
            "{simplified_member_snowflake}",
            "server",
            None,
            ["command.info.voice-mutes.server"],
        ),
        (
            "{simplified_member_snowflake}",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.command",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "auto",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.auto", "other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "click",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.click", "other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "command",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.command", "other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "server",
            "{other_guild_snowflake}",
            ["command.info.voice-mutes.server", "other_guilds"],
        ),
    ],
)
async def test_summary_app_command(
    bot,
    member: str,
    scope: str | None,
    other_guild: str | None,
    extra_permissions: list[str],
):
    docstring = """
    List bans, flags, text-mutes and voice-mutes on a member

    Parameters
    ----------
    member : str | int
        Resolves to: int | discord.Member
        Examples: 10000000000000010 | <@10000000000000010>

    guild (Optional) : str | int
        Resolves to: discord.Guild
        Examples: 10000000000000010

    Examples
    --------
    >>> !summary 10000000000000003
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
            command = cog.list_moderation_summary_app_command
            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
            )
            m = member.format(
                member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
            )
            if scope is None:
                s = scope
            else:
                s = scope
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
            target_transformer = AppTarget()
            scope_transformer = AppScope()
            resolved_member = await target_transformer.transform(inx, m)
            if g:
                resolved_guild = await target_transformer.transform(inx, g)
            else:
                resolved_guild = None
            if s:
                resolved_scope = await scope_transformer.transform(inx, s)
            else:
                resolved_scope = None
            async with capture_command() as end_results:
                await command.callback(
                    cog,
                    interaction=inx,
                    member=resolved_member,
                    scope=resolved_scope,
                    guild=resolved_guild,
                )
        for kind, content in end_results:
            assert kind == "success"
