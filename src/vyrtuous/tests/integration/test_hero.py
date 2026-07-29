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
from unittest.mock import AsyncMock, patch

import pytest

from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)

GUILD_SNOWFLAKE = 10000000000000500
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, member, target",
    [
        ("Guild_Owner", "hero", "{member_snowflake}", None),
        ("Guild_Owner", "hero", "{member_snowflake}", "all"),
        ("Guild_Owner", "hero", "<@{member_snowflake}>", "{guild_snowflake}"),
        ("Guild_Owner", "hero", "{simplified_member_snowflake}", "all"),
        (
            "Guild_Owner",
            "hero",
            "<@{simplified_member_snowflake}>",
            "{guild_snowflake}",
        ),
    ],
)
async def test_hero(bot, command: str, prefix: str, member, target, permission_role):
    """
    Promote or demote member to 'Hero' in memory (lost on reload).

    Parameters
    ----------
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is will be or is a hero
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------

    >>> !hero <@10000000000000003>
    [{emoji} Invincibility granted for Member1]

    >>> !hero 10000000000000003
    [{emoji} Invincibility granted for Member1]
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
    t = None
    if target:
        t = target.format(guild_snowflake=GUILD_SNOWFLAKE)
        full = f"{prefix}{command} {m} {t}"
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
            stack.enter_context(
                patch(
                    "vyrtuous.utils.permissions.permission_service.has_equal_or_lower_role",
                    new=AsyncMock(return_value=True),
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
            cog = bot.get_cog("UserManagementAppCommands")
            command = cog.toggle_invincibility_app_command
            transformer = AppTarget()
            resolved_member = await transformer.transform(inx, m)
            if t:
                resolved_target = await transformer.transform(inx, t)
            else:
                resolved_target = None
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
                stack.enter_context(
                    patch(
                        "vyrtuous.utils.permissions.permission_service.has_equal_or_lower_role",
                        new=AsyncMock(return_value=True),
                    )
                )

                await command.callback(
                    cog, interaction=inx, member=resolved_member, target=resolved_target
                )
        for kind, content in end_results:
            assert kind == "success"
