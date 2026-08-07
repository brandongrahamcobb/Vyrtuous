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

from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.autoassign import AutoAssignRole
from vyrtuous.models.group import AppGroup
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    send_message,
    setup,
)

VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
GUILD_SNOWFLAKE = 10000000000000500
OTHER_GUILD_SNOWFLAKE = 10000000000000501
ROLE_SNOWFLAKE = 10000000000000200


COMMAND = "autoassign"
BASE_PERMISSIONS = ["command.users.autoassign"]
TABLE_NAME = AutoAssignRole.__tablename__


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "elevated_group, group, role, channel, guild",
    [
        ("guild_owner", "administrator", "{role_snowflake}", None, None),
        (
            "guild_owner",
            "administrator",
            "{role_snowflake}",
            "{channel_snowflake}",
            None,
        ),
        (
            "guild_owner",
            "administrator",
            "{role_snowflake}",
            "{channel_snowflake}",
            "{guild_snowflake}",
        ),
        (
            "guild_owner",
            "administrator",
            "{role_snowflake}",
            "{channel_snowflake}",
            "{other_guild_snowflake}",
        ),
    ],
)
async def test_autoassign_text_command(
    bot,
    prefix: str,
    elevated_group: str,
    group: str,
    role: str,
    guild: str | None,
    channel: str | None,
):
    docstring = """
    Toggle autoassignment roles which are registered in the PostgreSQL database
    'vyrtuous' in the table 'autoassign_roles'.

    Parameters
    ----------
    group : str
        Resolved to: PermissionGroup
        Examples: administrator

    role : str | int
        Resolves to: discord.Role
        Examples: 10000000000000010 | <@&10000000000000010>

    channel (Optional) : str | int
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010> 

    guild (Optional) : int
        Resolves to: discord.Guild
        Examples: 10000000000000010 

    Example
    --------
    >>> !autoassign 10000000000000010
    Embed
    """
    assert TABLE_NAME in docstring
    assert COMMAND in docstring
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        permission_state = bot.registry.get(PermissionState)
        elevated_group = permission_state.groups.get(elevated_group)
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions",
                    side_effect=check_permissions(BASE_PERMISSIONS),
                )
            )
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.resolve_effective_group",
                    return_value=elevated_group,
                )
            )
            r = role.format(role_snowflake=ROLE_SNOWFLAKE)
            full = f"{prefix}{COMMAND} {group} {r}"
            if channel is None:
                c = channel
            else:
                c = channel.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                )
                full += f" {c}"
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
    "elevated_group, group, role, channel, guild",
    [
        ("guild_owner", "administrator", "{role_snowflake}", None, None),
        (
            "guild_owner",
            "administrator",
            "{role_snowflake}",
            "{channel_snowflake}",
            None,
        ),
        (
            "guild_owner",
            "administrator",
            "{role_snowflake}",
            "{channel_snowflake}",
            "{guild_snowflake}",
        ),
        (
            "guild_owner",
            "administrator",
            "{role_snowflake}",
            None,
            "{other_guild_snowflake}",
        ),
    ],
)
async def test_autoassign_app_command(
    bot,
    elevated_group: str,
    group: str,
    role: str,
    channel: str | None,
    guild: str | None,
):
    docstring = """
    Toggle autoassignment roles which are registered in the PostgreSQL database
    'vyrtuous' in the table 'autoassign_roles'.

    Parameters
    ----------
    group : str
        Resolved to: PermissionGroup
        Examples: administrator

    role : str | int
        Resolves to: discord.Role
        Examples: 10000000000000010 | <@&10000000000000010>

    guild (Optional) : int
        Resolves to: discord.Guild
        Examples: 10000000000000010 

    Example
    --------
    >>> /autoassign 10000000000000010
    Embed
    """
    assert TABLE_NAME in docstring
    assert COMMAND in docstring
    if (
        os.environ["TEST_MODE"].lower() == "app"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        permission_state = bot.registry.get(PermissionState)
        elevated_group = permission_state.groups.get(elevated_group)
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions",
                    side_effect=check_permissions(BASE_PERMISSIONS),
                )
            )
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.resolve_effective_group",
                    return_value=elevated_group,
                )
            )
            cog = bot.get_cog("UserManagementAppCommands")
            command = cog.toggle_autoassign_role_app_command
            r = role.format(role_snowflake=ROLE_SNOWFLAKE)
            if channel is None:
                c = channel
            else:
                c = channel.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                )
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
                group_transformer = AppGroup()
                resolved_group = await group_transformer.transform(inx, group)
                target_transformer = AppTarget()
                resolved_role = await target_transformer.transform(inx, r)
                if c:
                    resolved_channel = await target_transformer.transform(inx, c)
                else:
                    resolved_channel = None
                if g:
                    resolved_guild = await target_transformer.transform(inx, g)
                else:
                    resolved_guild = None
                await command.callback(
                    cog,
                    interaction=inx,
                    group=resolved_group,
                    role=resolved_role,
                    channel=resolved_channel,
                    guild=resolved_guild,
                )
            for kind, content in end_results:
                assert kind == "success"
