"""!/bin/python3
test_mods.py The purpose of this program is to be the integration test for the mods list command for Vyrtuous.

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
from datetime import datetime, timezone

import pytest

from vyrtuous.cache.registry import MemberState
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
VOICE_CHANNEL_SNOWFLAKE = 10000000000000011
DUMMY_MEMBER_SNOWFLAKE_TWO = 10000000000000005


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, target, guild",
    [
        ("Guild Owner", "heroes", "{guild_snowflake}", None),
        ("Guild Owner", "heroes", "{member_snowflake}", "{guild_snowflake}"),
        ("Guild Owner", "heroes", "<@{member_snowflake}>", None),
        ("Guild Owner", "heroes", "{simplified_member_snowflake}", "{guild_snowflake}"),
        ("Guild Owner", "heroes", "<@{simplified_member_snowflake}>", None),
    ],
)
async def test_mods(bot, command: str, prefix: str, target, guild, permission_role):
    """
    List members who are registered in the PostgresSQL database
    'vyrtuous' in the table 'moderators'.

    Parameters
    ----------
    all : str, optional
        Generic showing all moderators in all guilds
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with moderators
        in any of the guilds Vyrtuous has access inside.
    guild_snowflake : int | str, optional
        Snowflake of a guild where moderators are present.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is an moderator
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !mods "all"
    [{emoji} Heroes\n Guild1\n Guild2]

    >>> !mods <#10000000000000010>
    [{emoji} Heroes for Channel1\n Member1\n Member2]

    >>> !mods 10000000000000010
    [{emoji} Heroes for Channel1\n Member1\n Member2]

    >>> !mods 10000000000000500
    [{emoji} Heroes\n Guild1]

    >>> !mods <@10000000000000003>
    [{emoji} Heroes for Member1\n Guild1\n Guild2]

    >>> !mods 10000000000000003
    [{emoji} Heroes for Member1\n Guild1\n Guild2]
    """
    bot.registry.get(MemberState).active.update(
        {DUMMY_MEMBER_SNOWFLAKE_TWO: ("DUMMY", datetime.now(timezone.utc))}
    )
    t = target.format(
        guild_snowflake=GUILD_SNOWFLAKE,
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
        simplified_member_snowflake=DUMMY_MEMBER_SNOWFLAKE_TWO,
    )
    g = None
    if guild is None:
        full = f"{prefix}{command} {t}"
    else:
        g = guild.format(
            channel_snowflake=VOICE_CHANNEL_SNOWFLAKE,
            guild_snowflake=GUILD_SNOWFLAKE,
        )
        full = f"{prefix}{command} {t} {g}"
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
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
            cog = bot.get_cog("HiddenGuildOwnerAppCommands")
            command = cog.list_heroes_app_command
            transformer = AppTarget()
            if t:
                resolved = await transformer.transform(inx, t)
            else:
                resolved = None
            if g:
                resolved_guild = await transformer.transform(inx, g)
            else:
                resolved_guild = None
            await command.callback(
                cog, interaction=inx, target=resolved, guild=resolved_guild
            )
        for kind, content in end_results:
            assert kind == "success"
