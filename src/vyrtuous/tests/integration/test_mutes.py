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
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.models.scope import AppScope
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
COMMAND = "mutes"
BASE_PERMISSIONS = ["command.info.voice-mutes", "command.info.voice-mutes.command"]
TABLE_NAME = VoiceMute.__tablename__


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target, scope, other_guild, extra_permissions",
    [
        ("{channel_snowflake}", None, None, ["command.info.scope.channel"]),
        ("<#{channel_snowflake}>", None, None, ["command.info.scope.channel"]),
        (
            "{channel_snowflake}",
            "all",
            None,
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "<#{channel_snowflake}>",
            "all",
            None,
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{channel_snowflake}",
            "auto",
            None,
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "<#{channel_snowflake}>",
            "auto",
            None,
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "{channel_snowflake}",
            "click",
            None,
            ["command.info.scope.channel", "command.info.voice-mutes.click"],
        ),
        (
            "<#{channel_snowflake}>",
            "click",
            None,
            ["command.info.scope.channel", "command.info.voice-mutes.click"],
        ),
        (
            "{channel_snowflake}",
            "command",
            None,
            ["command.info.scope.channel"],
        ),
        (
            "<#{channel_snowflake}>",
            "command",
            None,
            ["command.info.scope.channel"],
        ),
        (
            "{channel_snowflake}",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "<#{channel_snowflake}>",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{channel_snowflake}",
            "click",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "<#{channel_snowflake}>",
            "click",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "{channel_snowflake}",
            "command",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "<#{channel_snowflake}>",
            "command",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "{member_snowflake}",
            None,
            None,
            ["command.info.scope.member", "other_channels"],
        ),
        (
            "<@{member_snowflake}>",
            None,
            None,
            ["command.info.scope.member", "other_channels"],
        ),
        (
            "{member_snowflake}",
            "all",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "all",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{member_snowflake}",
            "auto",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "auto",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "{member_snowflake}",
            "click",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "click",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "{member_snowflake}",
            "command",
            None,
            [
                "command.info.scope.member",
                "other_channels",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "command",
            None,
            [
                "command.info.scope.member",
                "other_channels",
            ],
        ),
        (
            "{member_snowflake}",
            "server",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "server",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{member_snowflake}",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{member_snowflake}",
            "click",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "click",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "{member_snowflake}",
            "command",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "command",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "{guild_snowflake}",
            None,
            None,
            ["command.info.scope.guild", "other_channels"],
        ),
        (
            "{guild_snowflake}",
            "all",
            None,
            [
                "command.info.scope.guild",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{guild_snowflake}",
            "auto",
            None,
            [
                "command.info.scope.guild",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "{guild_snowflake}",
            "click",
            None,
            [
                "command.info.scope.guild",
                "other_channels",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "{guild_snowflake}",
            "command",
            None,
            [
                "command.info.scope.guild",
                "other_channels",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.info.scope.member", "other_channels"],
        ),
        (
            "{simplified_member_snowflake}",
            "all",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "auto",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "click",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "command",
            None,
            [
                "command.info.scope.member",
                "other_channels",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "command",
            None,
            [
                "command.info.scope.member",
                "other_channels",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "click",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "command",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "other_guilds",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "server",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
    ],
)
async def test_mutes_text_command(
    bot,
    prefix: str,
    target: str | None,
    other_guild: str | None,
    scope: str | None,
    extra_permissions: list[str],
):
    docstring = """
    List voice-mutes on members which are registered in the PostgresSQL database
    'vyrtuous' in the table 'active_voice_mutes'.

    Parameters
    ----------
    target (Optional) : str | int
        Defaults to context channel
        Resolves to: int | discord.VoiceChannel | discord.TextChannel | discord.StageChannel | discord.Member | discord.Guild
        Examples: 10000000000000010 | <#10000000000000010> | <@10000000000000010>

    scope (Optiona) : str
        Default to all
        Resolves to: ScopeObject
        Examples: all | auto | click | command | server

    guild (Optional) : str | int
        Resolves to: discord.Guild
        Examples: 10000000000000010

    Example
    --------
    >>> !mutes
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
                    member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                    simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
                )
                full += f" {t}"
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
            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
            )
            captured = await send_message(bot=bot, content=full)
            assert captured["end_result"] == ["success"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "target, scope, other_guild, extra_permissions",
    [
        ("{channel_snowflake}", None, None, ["command.info.scope.channel"]),
        ("<#{channel_snowflake}>", None, None, ["command.info.scope.channel"]),
        (
            "{channel_snowflake}",
            "all",
            None,
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "<#{channel_snowflake}>",
            "all",
            None,
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{channel_snowflake}",
            "auto",
            None,
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
            ],
        ),
        (
            "<#{channel_snowflake}>",
            "auto",
            None,
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
            ],
        ),
        (
            "{channel_snowflake}",
            "click",
            None,
            ["command.info.scope.channel", "command.info.voice-mutes.click"],
        ),
        (
            "<#{channel_snowflake}>",
            "click",
            None,
            ["command.info.scope.channel", "command.info.voice-mutes.click"],
        ),
        (
            "{channel_snowflake}",
            "command",
            None,
            ["command.info.scope.channel"],
        ),
        (
            "<#{channel_snowflake}>",
            "command",
            None,
            ["command.info.scope.channel"],
        ),
        (
            "{channel_snowflake}",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "<#{channel_snowflake}>",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{channel_snowflake}",
            "click",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "<#{channel_snowflake}>",
            "click",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "{channel_snowflake}",
            "command",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "<#{channel_snowflake}>",
            "command",
            "{other_guild_snowflake}",
            [
                "command.info.scope.channel",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "{member_snowflake}",
            None,
            None,
            ["command.info.scope.member", "other_channels"],
        ),
        (
            "<@{member_snowflake}>",
            None,
            None,
            ["command.info.scope.member", "other_channels"],
        ),
        (
            "{member_snowflake}",
            "all",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "all",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{member_snowflake}",
            "auto",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "auto",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
            ],
        ),
        (
            "{member_snowflake}",
            "click",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "click",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "{member_snowflake}",
            "command",
            None,
            [
                "command.info.scope.member",
                "other_channels",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "command",
            None,
            [
                "command.info.scope.member",
                "other_channels",
            ],
        ),
        (
            "{member_snowflake}",
            "server",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "server",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{member_snowflake}",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{member_snowflake}",
            "click",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "click",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "{member_snowflake}",
            "command",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.command",
                "other_guilds",
            ],
        ),
        (
            "<@{member_snowflake}>",
            "command",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.command",
                "other_guilds",
            ],
        ),
        (
            "{guild_snowflake}",
            None,
            None,
            ["command.info.scope.guild", "other_channels"],
        ),
        (
            "{guild_snowflake}",
            "all",
            None,
            [
                "command.info.scope.guild",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{guild_snowflake}",
            "auto",
            None,
            [
                "command.info.scope.guild",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{guild_snowflake}",
            "click",
            None,
            [
                "command.info.scope.guild",
                "other_channels",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "{guild_snowflake}",
            "command",
            None,
            [
                "command.info.scope.guild",
                "other_channels",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.info.scope.member", "other_channels"],
        ),
        (
            "{simplified_member_snowflake}",
            "all",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "auto",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "click",
            None,
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "command",
            None,
            [
                "command.info.scope.member",
                "other_channels",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "all",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.auto",
                "command.info.voice-mutes.click",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "click",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.click",
                "other_guilds",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "command",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "other_guilds",
            ],
        ),
        (
            "{simplified_member_snowflake}",
            "server",
            "{other_guild_snowflake}",
            [
                "command.info.scope.member",
                "other_channels",
                "command.info.voice-mutes.server",
                "other_guilds",
            ],
        ),
    ],
)
async def test_mutes_app_command(
    bot,
    target: str | None,
    other_guild: str | None,
    scope: str | None,
    extra_permissions: list[str],
):
    docstring = """
    List voice-mutes on members which are registered in the PostgresSQL database
    'vyrtuous' in the table 'active_voice_mutes'.

    Parameters
    ----------
    target (Optional) : str | int
        Defaults to context channel
        Resolves to: int | discord.VoiceChannel | discord.TextChannel | discord.StageChannel | discord.Member | discord.Guild
        Examples: 10000000000000010 | <#10000000000000010> | <@10000000000000010>

    scope (Optiona) : str
        Default to all
        Resolves to: ScopeObject
        Examples: all | auto | click | command | server

    guild (Optional) : str | int
        Resolves to: discord.Guild
        Examples: 10000000000000010

    Example
    --------
    >>> !mutes
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
            command = cog.list_mutes_app_command
            if target is None:
                t = target
            else:
                t = target.format(
                    channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
                    guild_snowflake=GUILD_SNOWFLAKE,
                    member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                    simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
                )
            s = scope
            if other_guild is None:
                g = other_guild
            else:
                g = other_guild.format(other_guild_snowflake=OTHER_GUILD_SNOWFLAKE)
            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
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
            transformer = AppTarget()
            if t:
                resolved_target = await transformer.transform(inx, t)
            else:
                resolved_target = None
            if g:
                resolved_guild = await transformer.transform(inx, g)
            else:
                resolved_guild = None
            scope_transformer = AppScope()
            if s:
                resolved_scope = await scope_transformer.transform(inx, s)
            else:
                resolved_scope = None
            async with capture_command() as end_results:
                await command.callback(
                    cog,
                    interaction=inx,
                    target=resolved_target,
                    scope=resolved_scope,
                    guild=resolved_guild,
                )
        for kind, content in end_results:
            assert kind == "success"


COLUMNS = [
    ("channel_snowflake", "bigint", True),
    ("guild_snowflake", "bigint", False),
    ("member_snowflake", "bigint", False),
    ("expires_in", "timestamp with time zone", True),
    ("target", "text", False),
    ("created_at", "timestamp with time zone", True),
    ("updated_at", "timestamp with time zone", True),
    ("reason", "text", False),
]


@pytest.mark.asyncio
@pytest.mark.parametrize("field, datatype, nullable", COLUMNS)
async def test_active_voice_mutes_database_table(
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
