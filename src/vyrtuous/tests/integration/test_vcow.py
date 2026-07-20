"""!/bin/python3
test_coord.py The purpose of this program is to be the integration test for the coord demotion/promotion command for Vyrtuous.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import os
from contextlib import ExitStack
from unittest.mock import AsyncMock, patch

import pytest

from vyrtuous.models.target import AppTarget
from vyrtuous.tests.conftest import interaction
from vyrtuous.tests.integration.test_suite import (
    build_message,
    capture_command,
    send_message,
    setup,
)

DUMMY_MEMBER_SNOWFLAKE = 10000000000000003
GUILD_SNOWFLAKE = 10000000000000500


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, member, guild",
    [
        ("Moderator", "!vcow", "{member_snowflake}", None),
        ("Moderator", "!vcow", "<@{member_snowflake}>", "{guild_snowflake}"),
    ],
)
async def test_vcow(bot, command: str, member, guild, permission_role):
    """
    Promote or demote a member with 'Coordinator' by registering them in the PostgresSQL database
    'vyrtuous' in the table 'coordinators'.

    Parameters
    ----------
    channel_snowflake : int | str, optional
        Mention or snowflake of a channel with coordinators
        in any of the guilds Vyrtuous has access inside.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is an coordinator
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------

    >>> !coord <@10000000000000003> <@10000000000000010>
    [{emoji} Coordinator granted for Member1]

    >>> !coord 10000000000000003 10000000000000010
    [{emoji} Coordinator granted for Member1]
    """
    m = member.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE,
    )
    g = None
    full = f"{command} {m} {GUILD_SNOWFLAKE}"
    if guild:
        g = guild.format(guild_snowflake=GUILD_SNOWFLAKE)
    if (
        os.environ["TEST_MODE"].lower() == "text"
        or os.environ["TEST_MODE"].lower() == "all"
    ):
        captured = await send_message(bot=bot, content=full)
        assert captured == ["success"]
    elif (
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
        with ExitStack() as stack:
            stack.enter_context(
                patch(
                    "vyrtuous.tests.integration.mock_discord_channel.MockVoiceChannel.fetch_message",
                    new=AsyncMock(return_value=msg),
                )
            )
            async with capture_command() as end_results:
                cog = bot.get_cog("HiddenModeratorAppCommands")
                command = cog.toggle_vegan_app_command
                transformer = AppTarget()
                resolved_member = await transformer.transform(inx, m)
                if g:
                    resolved_guild = await transformer.transform(inx, g)
                else:
                    resolved_guild = None
                await command.callback(
                    cog,
                    interaction=inx,
                    member=resolved_member,
                    guild=resolved_guild,
                    notes="",
                )
            for kind, content in end_results:
                print(content)
                assert kind == "success"
