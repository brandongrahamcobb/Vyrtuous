"""!/bin/python3
test_hero.py The purpose of this program is to be the integration test for the hero promotion command for Vyrtuous.

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
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005

BASE_PERMISSIONS = ["command.users.hero"]
COMMAND = "hero"


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "member, target, extra_permissions",
    [
        ("{member_snowflake}", None, []),
        ("<@{member_snowflake}>", None, []),
        ("{member_snowflake}", "all", ["other_guilds"]),
        ("<@{member_snowflake}>", "all", ["other_guilds"]),
        ("{member_snowflake}", "{guild_snowflake}", []),
        ("<@{member_snowflake}>", "{guild_snowflake}", []),
        ("{member_snowflake}", "{other_guild_snowflake}", ["other_guilds"]),
        ("<@{member_snowflake}>", "{other_guild_snowflake}", ["other_guilds"]),
        ("{simplified_member_snowflake}", None, []),
        ("{simplified_member_snowflake}", "all", ["other_guilds"]),
        ("{simplified_member_snowflake}", "{guild_snowflake}", []),
        ("{simplified_member_snowflake}", "{other_guild_snowflake}", ["other_guilds"]),
    ],
)
async def test_hero_text_command(
    bot,
    prefix: str,
    member: str,
    target: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Toggle a hero (removes all infractions and keeps them from being moderated)
    in the invincible MemberState cache.

    Parameters
    ----------
    member : str | int
        Resolves to: int | discord.Member
        Excamples: 10000000000000010 | <@10000000000000010>
    target (Optional) : str | int | None
        Resolves to: str | discord.Guild
        Examples: all | 10000000000000010

    Example
    --------
    >>> !hero
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
            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
            )
            m = member.format(
                member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
            )
            full = f"{prefix}{COMMAND} {m}"
            if target is None:
                t = target
            else:
                t = target.format(
                    guild_snowflake=GUILD_SNOWFLAKE,
                    other_guild_snowflake=OTHER_GUILD_SNOWFLAKE,
                )
                full += f" {t}"
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "member, target, extra_permissions",
    [
        ("{member_snowflake}", None, []),
        ("<@{member_snowflake}>", None, []),
        ("{member_snowflake}", "all", ["other_guilds"]),
        ("<@{member_snowflake}>", "all", ["other_guilds"]),
        ("{member_snowflake}", "{guild_snowflake}", []),
        ("<@{member_snowflake}>", "{guild_snowflake}", []),
        ("{member_snowflake}", "{other_guild_snowflake}", ["other_guilds"]),
        ("<@{member_snowflake}>", "{other_guild_snowflake}", ["other_guilds"]),
        ("{simplified_member_snowflake}", None, []),
        ("{simplified_member_snowflake}", "all", ["other_guilds"]),
        ("{simplified_member_snowflake}", "{guild_snowflake}", []),
        ("{simplified_member_snowflake}", "{other_guild_snowflake}", ["other_guilds"]),
    ],
)
async def test_hero_app_command(
    bot,
    member: str,
    target: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Toggle a hero (removes all infractions and keeps them from being moderated)
    in the invincible MemberState cache.

    Parameters
    ----------
    member : str | int
        Resolves to: int | discord.Member
        Excamples: 10000000000000010 | <@10000000000000010>
    target (Optional) : str | int | None
        Resolves to: str | discord.Guild
        Examples: all | 10000000000000010

    Example
    --------
    >>> /hero
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
            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
            )
            cog = bot.get_cog("UserManagementAppCommands")
            command = cog.toggle_hero_app_command
            m = member.format(
                member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
            )
            if target is None:
                t = target
            else:
                t = target.format(
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
                resolved_member = await transformer.transform(inx, m)
                if t:
                    resolved_target = await transformer.transform(inx, t)
                else:
                    resolved_target = None
                await command.callback(
                    cog, interaction=inx, member=resolved_member, target=resolved_target
                )
            for kind, content in end_results:
                assert kind == "success"
