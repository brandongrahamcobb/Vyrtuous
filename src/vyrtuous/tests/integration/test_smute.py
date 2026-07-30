"""!/bin/python3
test_smute.py The purpose of this program is to be the integration test for the smute command for Vyrtuous.

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

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
OTHER_GUILD_SNOWFLAKE = 10000000000000501
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005

BASE_PERMISSIONS = [
    "command.moderation.voice-mute.server",
]
COMMAND = "smute"


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "member, reason, other_guild, extra_permissions",
    [
        ("{member_snowflake}", None, None, []),
        ("<@{member_snowflake}>", None, None, []),
        (
            "{member_snowflake}",
            "test reason",
            "{other_guild_snowflake}",
            ["other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "test reason",
            "{other_guild_snowflake}",
            ["other_guilds"],
        ),
        ("{simplified_member_snowflake}", None, None, []),
        ("<#{simplified_member_snowflake}>", None, None, []),
        (
            "{simplified_member_snowflake}",
            "test reason",
            "{other_guild_snowflake}",
            ["other_guilds"],
        ),
        (
            "<#{simplified_member_snowflake}>",
            "test reason",
            "{other_guild_snowflake}",
            ["other_guilds"],
        ),
    ],
)
async def test_smute_app_command(
    bot,
    prefix: str,
    member: str | None,
    other_guild: str | None,
    reason: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Server mute a member. 

    Parameters
    ----------
    member : str | int
        Resolves to: int | discord.Member
        Examples: 10000000000000010 | <@10000000000000010>

    reason (Optional) : str
        Examples: test reason

    Example
    --------
    >>> !smute 10000000000000010 
    Successfully server muted Spawd in The Vyrtuous Network.
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
            full = f"{prefix}{COMMAND}"
            if member is None:
                m = member
            else:
                m = member.format(
                    member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                    simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
                )
                full += f" {m}"
            if other_guild is None:
                g = other_guild
            else:
                g = other_guild.format(other_guild_snowflake=OTHER_GUILD_SNOWFLAKE)
                full += f" {g}"
            if reason is None:
                r = reason
            else:
                r = reason
                full += f" {r}"

            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
            )
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "member, reason, other_guild, extra_permissions",
    [
        ("{member_snowflake}", None, None, []),
        ("<@{member_snowflake}>", None, None, []),
        ("{member_snowflake}", "test reason", None, []),
        ("<@{member_snowflake}>", "test reason", None, []),
        (
            "{member_snowflake}",
            "test reason",
            "{other_guild_snowflake}",
            ["other_guilds"],
        ),
        (
            "<@{member_snowflake}>",
            "test reason",
            "{other_guild_snowflake}",
            ["other_guilds"],
        ),
        ("{simplified_member_snowflake}", None, None, []),
        ("<#{simplified_member_snowflake}>", None, None, []),
        ("{simplified_member_snowflake}", "test reason", None, []),
        ("<#{simplified_member_snowflake}>", "test reason", None, []),
        (
            "{simplified_member_snowflake}",
            "test reason",
            "{other_guild_snowflake}",
            ["other_guilds"],
        ),
        (
            "<#{simplified_member_snowflake}>",
            "test reason",
            "{other_guild_snowflake}",
            ["other_guilds"],
        ),
    ],
)
async def test_smute_text_command(
    bot,
    member: str | None,
    other_guild: str | None,
    reason: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Server mute a member. 

    Parameters
    ----------
    member : str | int
        Resolves to: int | discord.Member
        Examples: 10000000000000010 | <@10000000000000010>

    reason (Optional) : str
        Examples: test reason

    Example
    --------
    >>> !smute 10000000000000010 
    Successfully server muted Spawd in The Vyrtuous Network.
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
            cog = bot.get_cog("ModerationAppCommands")
            command = cog.toggle_server_mute_app_command
            if member is None:
                m = member
            else:
                m = member.format(
                    member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                    simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
                )
            if other_guild is None:
                g = other_guild
            else:
                g = other_guild.format(other_guild_snowflake=OTHER_GUILD_SNOWFLAKE)
            if reason is None:
                r = reason
            else:
                r = reason
            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
            )
            objects = setup(bot)
            msg = build_message(
                author=objects.get("author", None),
                channel=objects.get("voice_channel", None),
                content="",
                guild=objects.get("guild", None),
                state=objects.get("state", None),
            )
            inx = interaction(
                bot=bot,
                channel=objects.get("voice_channel", None),
                guild=objects.get("guild", None),
                message=msg,
            )
            async with capture_command() as end_results:
                target_transformer = AppTarget()
                if m:
                    resolved_member = await target_transformer.transform(inx, m)
                else:
                    resolved_member = m
                if g:
                    resolved_guild = await target_transformer.transform(inx, g)
                else:
                    resolved_guild = m

                await command.callback(
                    cog,
                    interaction=inx,
                    member=resolved_member,
                    reason=r,
                    guild=resolved_guild,
                )
            for kind, content in end_results:
                assert kind == "success"
