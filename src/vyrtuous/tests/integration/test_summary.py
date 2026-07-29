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
from unittest.mock import AsyncMock, patch

import pytest

from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.models.scope import AppScope
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
OTHER_GUILD_SNOWFLAKE = 10000000000000501
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, member, scope, guild",
    [
        ("Moderator", "summary", "{member_snowflake}", None, None),
        ("Moderator", "summary", "{member_snowflake}", "all", None),
        ("Moderator", "summary", "{member_snowflake}", "click", None),
        ("Moderator", "summary", "{member_snowflake}", "command", None),
        (
            "Moderator",
            "summary",
            "<@{member_snowflake}>",
            "click",
            "{other_guild_snowflake}",
        ),
        ("Moderator", "summary", "{simplified_member_snowflake}", None, None),
        ("Moderator", "summary", "{simplified_member_snowflake}", "all", None),
        ("Moderator", "summary", "{simplified_member_snowflake}", "click", None),
        ("Moderator", "summary", "{simplified_member_snowflake}", "command", None),
        (
            "Moderator",
            "summary",
            "<@{simplified_member_snowflake}>",
            "click",
            "{other_guild_snowflake}",
        ),
    ],
)
async def test_summary(
    bot, command: str, prefix: str, member, guild, scope, permission_role
):
    """
    List voice-mutes on members which are registered in the PostgresSQL database
    'vyrtuous' in the table 'active_voice_mutes'.

    Parameters
    ----------
    all : str, optional
        Generic showing all voice mutes in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with voice-mutes on members
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where mutes are present.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who has been voice-muted
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !summary <@10000000000000003>
    [{emoji} Infractions for Member1\n Guild1\n Guild2]

    >>> !summary 10000000000000003
    [{emoji} Infractions for Member1\n Guild1\n Guild2]
    """
    permission_state = bot.registry.get(PermissionState)
    bot.registry.get(MemberState).active.update(
        {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
    )
    m = member.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
        simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
    )
    full = f"{prefix}{command} {m}"
    if guild is None and not scope:
        g = None
        s = None
    elif scope and not guild:
        g = None
        s = scope
        full = f"{prefix}{command} {m} {s}"
    elif guild:
        g = guild.format(other_guild_snowflake=OTHER_GUILD_SNOWFLAKE)
        s = scope
        full = f"{prefix}{command} {m} {s} {g}"
    else:
        g = None
        s = None
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.utils.permissions.permission_service.resolve_effective_group",
                    new=AsyncMock(
                        return_value=permission_state.groups.get(
                            permission_role.lower()
                        )
                    ),
                )
            )
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]
    if (
        os.environ["TEST_MODE"].lower() == "app"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        objects = setup(bot)
        msg = build_message(
            author=objects.get("author", None),
            channel=objects.get("text_channel", None),
            content=full,
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
            cog = bot.get_cog("InfoAppCommands")
            command = cog.list_moderation_summary_app_command
            transformer = AppTarget()
            resolved_member = await transformer.transform(inx, m)
            if g:
                resolved_guild = await transformer.transform(inx, g)
            else:
                resolved_guild = None
            scope_transformer = AppScope()
            if s:
                resolved_scope = await scope_transformer.transform(inx, s)
            else:
                resolved_scope = None
            with ExitStack() as stack:
                stack.enter_context(
                    patch(
                        "vyrtuous.utils.permissions.permission_service.resolve_effective_group",
                        new=AsyncMock(
                            return_value=permission_state.groups.get(
                                permission_role.lower()
                            )
                        ),
                    )
                )
                await command.callback(
                    cog,
                    interaction=inx,
                    member=resolved_member,
                    scope=resolved_scope,
                    guild=resolved_guild,
                )
        for kind, content in end_results:
            assert kind == "success"
