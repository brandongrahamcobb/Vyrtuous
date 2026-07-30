"""!/bin/python3
test_cogs.py The purpose of this program is to be the integration test for the cogs list command for Vyrtuous.

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

from vyrtuous.cache.registry import MemberState
from vyrtuous.tests.integration.test_suite import (
    build_message,
    check_permissions,
    send_message,
    setup,
)

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command, member, duration, reason, extra_permissions",
    [
        ("testban", "{member_snowflake}", None, None, ["command.moderation.ban"]),
        ("testban", "{member_snowflake}", None, None, ["command.moderation.ban"]),
        ("testban", "{member_snowflake}", "1h", None, ["command.moderation.ban"]),
        ("testban", "{member_snowflake}", None, None, ["command.moderation.ban"]),
        (
            "testban",
            "{member_snowflake}",
            "1h",
            "test reason",
            ["command.moderation.ban"],
        ),
        ("testban", "{member_snowflake}", None, None, ["command.moderation.ban"]),
        (
            "testban",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.ban"],
        ),
        (
            "testban",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.ban"],
        ),
        (
            "testban",
            "{simplified_member_snowflake}",
            "1h",
            None,
            ["command.moderation.ban"],
        ),
        (
            "testban",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.ban"],
        ),
        (
            "testban",
            "{simplified_member_snowflake}",
            "1h",
            "test reason",
            ["command.moderation.ban"],
        ),
        (
            "testban",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.ban"],
        ),
        ("testban", "<@{member_snowflake}>", None, None, ["command.moderation.ban"]),
        ("testban", "<@{member_snowflake}>", None, None, ["command.moderation.ban"]),
        ("testban", "<@{member_snowflake}>", "1h", None, ["command.moderation.ban"]),
        ("testban", "<@{member_snowflake}>", None, None, ["command.moderation.ban"]),
        (
            "testban",
            "<@{member_snowflake}>",
            "1h",
            "test reason",
            ["command.moderation.ban"],
        ),
        ("testban", "<@{member_snowflake}>", None, None, ["command.moderation.ban"]),
        ("testflag", "{member_snowflake}", None, None, ["command.moderation.flag"]),
        ("testflag", "{member_snowflake}", None, None, ["command.moderation.flag"]),
        (
            "testflag",
            "{member_snowflake}",
            None,
            "test reason",
            ["command.moderation.flag"],
        ),
        ("testflag", "{member_snowflake}", None, None, ["command.moderation.flag"]),
        (
            "testflag",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.flag"],
        ),
        (
            "testflag",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.flag"],
        ),
        (
            "testflag",
            "{simplified_member_snowflake}",
            None,
            "test reason",
            ["command.moderation.flag"],
        ),
        (
            "testflag",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.flag"],
        ),
        ("testflag", "<@{member_snowflake}>", None, None, ["command.moderation.flag"]),
        ("testflag", "<@{member_snowflake}>", None, None, ["command.moderation.flag"]),
        (
            "testflag",
            "<@{member_snowflake}>",
            None,
            "test reason",
            ["command.moderation.flag"],
        ),
        ("testflag", "<@{member_snowflake}>", None, None, ["command.moderation.flag"]),
        (
            "testmute",
            "{member_snowflake}",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{member_snowflake}",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{member_snowflake}",
            "1h",
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{member_snowflake}",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{member_snowflake}",
            "1h",
            "test reason",
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{member_snowflake}",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{simplified_member_snowflake}",
            "1h",
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{simplified_member_snowflake}",
            "1h",
            "test reason",
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "<@{member_snowflake}>",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "<@{member_snowflake}>",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "<@{member_snowflake}>",
            "1h",
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "<@{member_snowflake}>",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "<@{member_snowflake}>",
            "1h",
            "test reason",
            ["command.moderation.voice-mute"],
        ),
        (
            "testmute",
            "<@{member_snowflake}>",
            None,
            None,
            ["command.moderation.voice-mute"],
        ),
        (
            "testtmute",
            "{member_snowflake}",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{member_snowflake}",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{member_snowflake}",
            "1h",
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{member_snowflake}",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{member_snowflake}",
            "1h",
            "test reason",
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{member_snowflake}",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{simplified_member_snowflake}",
            "1h",
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{simplified_member_snowflake}",
            "1h",
            "test reason",
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "{simplified_member_snowflake}",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "<@{member_snowflake}>",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "<@{member_snowflake}>",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "<@{member_snowflake}>",
            "1h",
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "<@{member_snowflake}>",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "<@{member_snowflake}>",
            "1h",
            "test reason",
            ["command.moderation.text-mute"],
        ),
        (
            "testtmute",
            "<@{member_snowflake}>",
            None,
            None,
            ["command.moderation.text-mute"],
        ),
        ("testrole", "{member_snowflake}", None, None, ["command.moderation.role"]),
        ("testrole", "{member_snowflake}", None, None, ["command.moderation.role"]),
        ("testrole", "<@{member_snowflake}>", None, None, ["command.moderation.role"]),
        ("testrole", "<@{member_snowflake}>", None, None, ["command.moderation.role"]),
    ],
)
async def test_alias_text_commands(
    bot,
    prefix: str,
    command: str,
    member: str,
    duration: str | None,
    reason: str | None,
    extra_permissions: list[str],
):
    docstring = """
    Toggle a ban, flag, role, text-mute or voice-mute.

    Parameters
    ----------
    member : str | int
        Resolves to: int | discord.Member
        Examples: 10000000000000010 | <@10000000000000010>

    duration (Optional): str
        Resolves to: DurationObject
        Examples: 0 | 1 | 30m | 1h | 1d | 1w 

    reason (Optional): str
        Examples: test reason

    Examples
    --------
    >>> !testban 10000000000000010
    Embed
    >>> !testflag 10000000000000010
    Embed
    >>> !testrole 10000000000000010
    Embed
    >>> !testtmute 10000000000000010
    Embed
    >>> !testmute 10000000000000010
    Embed
    """
    assert command in docstring
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        objects = setup(bot)
        msg = build_message(
            author=objects.get("author", None),
            channel=objects.get("text_channel", None),
            content="",
            guild=objects.get("guild", None),
            state=objects.get("state", None),
        )
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.permissions.permission_service.has_permissions",
                    side_effect=check_permissions(extra_permissions),
                )
            )
            stack.enter_context(
                patch(
                    "vyrtuous.tests.integration.mock_discord_channel.MockVoiceChannel.fetch_message",
                    new=AsyncMock(return_value=msg),
                )
            )
            bot.registry.get(MemberState).active.update(
                {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
            )
            m = member.format(
                member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
                simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
            )
            full = f"{prefix}{command} {m}"
            if duration is None:
                d = duration
            else:
                d = duration
                full += f" {d}"
            if reason is None:
                r = reason
            else:
                r = reason
                full += f" {r}"
            captured = await send_message(bot=bot, content=full)
            assert captured == ["success"]
