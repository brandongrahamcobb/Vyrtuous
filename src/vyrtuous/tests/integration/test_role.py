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

from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    check_permissions,
    setup,
)
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.view.role_view import RoleView

VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
OTHER_GUILD_CHANNEL_SNOWFLAKE = 10000000000000013
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005
OTHER_GUILD_SNOWFLAKE = 10000000000000501
GUILD_SNOWFLAKE = 10000000000000500

COMMAND = "role"
BASE_PERMISSIONS = ["command.moderation.role"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "member, channel, group_alias, extra_permissions",
    [
        ("{member_snowflake}", None, "moderator", []),
        ("{member_snowflake}", None, "coordinator", []),
        ("{member_snowflake}", None, "administrator", []),
        ("{member_snowflake}", None, "guild_owner", []),
        ("{member_snowflake}", None, "developer", []),
        ("{member_snowflake}", None, "sysadmin", []),
        ("{member_snowflake}", "{channel_snowflake}", "moderator", []),
        ("{member_snowflake}", "{channel_snowflake}", "coordinator", []),
        ("{member_snowflake}", "{channel_snowflake}", "administrator", []),
        ("{member_snowflake}", "{channel_snowflake}", "guild_owner", []),
        ("{member_snowflake}", "{channel_snowflake}", "developer", []),
        ("{member_snowflake}", "{channel_snowflake}", "sysadmin", []),
        (
            "{simplified_member_snowflake}",
            "{channel_snowflake}",
            "moderator",
            [],
        ),
        (
            "{simplified_member_snowflake}",
            "{channel_snowflake}",
            "coordinator",
            [],
        ),
        (
            "{simplified_member_snowflake}",
            "{channel_snowflake}",
            "administrator",
            [],
        ),
        (
            "{simplified_member_snowflake}",
            "{channel_snowflake}",
            "guild_owner",
            [],
        ),
        (
            "{simplified_member_snowflake}",
            "{channel_snowflake}",
            "developer",
            [],
        ),
        (
            "{simplified_member_snowflake}",
            "{channel_snowflake}",
            "sysadmin",
            [],
        ),
        (
            "{simplified_member_snowflake}",
            "{other_guild_channel_snowflake}",
            "moderator",
            ["other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "{other_guild_channel_snowflake}",
            "coordinator",
            ["other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "{other_guild_channel_snowflake}",
            "administrator",
            ["other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "{other_guild_channel_snowflake}",
            "guild_owner",
            ["other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "{other_guild_channel_snowflake}",
            "developer",
            ["other_guilds"],
        ),
        (
            "{simplified_member_snowflake}",
            "{other_guild_channel_snowflake}",
            "sysadmin",
            ["other_guilds"],
        ),
    ],
)
async def test_role_app_command(
    bot,
    member: str,
    channel: str | None,
    group_alias: str,
    extra_permissions: list[str],
):
    docstring = """
    Role/unrole a member

    Parameters
    ----------
    member : str | int
        Resolves to: int | discord.Member
        Examples: 10000000000000010 | <@10000000000000010>

    channel (Optional) : str | int
        Resolves to: discord.VoiceChannel | discord.TextChannel | discord.StageChannel
        Examples: 10000000000000010 | <#10000000000000010>

    Examples
    --------
    >>> /role 10000000000000003
    Embed
    """
    assert COMMAND in docstring
    if (
        os.environ["TEST_MODE"].lower() == "app"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        extra_permissions.extend(BASE_PERMISSIONS)
        database_factory: DatabaseFactory = DatabaseFactory(VoiceMute)
        group = bot.registry.get(PermissionState).groups.get(group_alias)
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
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.resolve_effective_group",
                    return_value=group,
                )
            )
            cog = bot.get_cog("UserManagementAppCommands")
            command = cog.toggle_role_app_command
            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
            )
            m = member.format(
                member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
            )
            if channel is None:
                c = channel
                g = GUILD_SNOWFLAKE
                voice_mute = VoiceMute(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    guild_snowflake=g,
                    member_snowflake=int(m),
                    target="command",
                )
                await database_factory.create(voice_mute)
                if group_alias in ["moderator", "coordinator"]:
                    return_value = [
                        (
                            group,
                            GUILD_SNOWFLAKE,
                            VOICE_CHANNEL_SNOWFLAKE,
                        )
                    ]
                elif group_alias in ["administrator", "guild_owner"]:
                    return_value = [
                        (
                            group,
                            GUILD_SNOWFLAKE,
                            None,
                        )
                    ]
                elif group_alias in ["developer", "sysadmin"]:
                    return_value = [
                        (
                            group,
                            None,
                            None,
                        )
                    ]
                else:
                    assert False
                stack.enter_context(
                    patch(
                        "vyrtuous.permissions.permission_service.resolve_all_assigned_groups",
                        return_value=return_value,
                    )
                )

            else:
                c = channel.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    other_guild_channel_snowflake=OTHER_GUILD_CHANNEL_SNOWFLAKE,
                )
                g = (
                    GUILD_SNOWFLAKE
                    if c == str(VOICE_CHANNEL_SNOWFLAKE)
                    else OTHER_GUILD_SNOWFLAKE
                )
                voice_mute = VoiceMute(
                    channel_snowflake=int(c),
                    guild_snowflake=int(g),
                    member_snowflake=int(m),
                    target="command",
                )
                await database_factory.create(voice_mute)
                if group_alias in ["moderator", "coordinator"]:
                    return_value = [
                        (
                            group,
                            g,
                            c,
                        )
                    ]
                elif group_alias in ["administrator", "guild_owner"]:
                    return_value = [
                        (
                            group,
                            g,
                            None,
                        )
                    ]
                elif group_alias in ["developer", "sysadmin"]:
                    return_value = [
                        (
                            group,
                            None,
                            None,
                        )
                    ]
                else:
                    assert False
                stack.enter_context(
                    patch(
                        "vyrtuous.permissions.permission_service.resolve_all_assigned_groups",
                        return_value=return_value,
                    )
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
            target_transformer = AppTarget()
            resolved_member = await target_transformer.transform(inx, m)
            if c:
                resolved_channel = await target_transformer.transform(inx, c)
            else:
                resolved_channel = None
            async with capture_command() as end_results:
                await command.callback(
                    cog,
                    interaction=inx,
                    member=resolved_member,
                    channel=resolved_channel,
                )
            tick = Tick(bot=bot, interaction=inx)
            ctx = SnowflakeContext(
                channel_snowflake=int(c or VOICE_CHANNEL_SNOWFLAKE),
                guild_snowflake=g,
                member_snowflake=int(m),
            )
            author = objects.get("author", None)
            if author is not None:
                view = RoleView(
                    author_snowflake=author.id,
                    ctx=ctx,
                    interaction=inx,
                    tick=tick,
                )
                await view.setup()
        for kind, content in end_results:
            assert kind == "success"
