"""!/bin/python3
test_admins.py The purpose of this program is to be the integration test for the admins list command for Vyrtuous.

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

import pytest

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
NOT_PRIVILEGED_AUTHOR_SNOWFLAKE_ONE = 10000000000000002
NOT_PRIVILEGED_AUTHOR_NAME_ONE = "Not Privileged Author Name One"


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission_role, command, target, guild",
    [
        ("Administrator", "admins", "{guild_snowflake}", None),
        ("Moderator", "admins", "{member_snowflake}", None),
        ("Moderator", "admins", "<@{member_snowflake}>", None),
        ("Moderator", "admins", "<@{member_snowflake}>", "{guild_snowflake}"),
    ],
)
async def test_admins(bot, command: str, prefix: str, target, guild, permission_role):
    """
    List members who are registered in the PostgresSQL database
    'vyrtuous' in the table 'administrators'.

    Parameters
    ----------
    all : str, optional
        Generic showing all administrators in all guilds
    guild_snowflake : int | str, optional
        Snowflake of a guild where administrators are present.
    member_snowflake : int | str, optional
        Mention or snowflake of a member who is an administrator
        in any of the guilds Vyrtuous has access inside.

    Examples
    --------
    >>> !admins "all"
    [{emoji} Administrators\n Guild1\n Guild2]

    >>> !admins 10000000000000500
    [{emoji} Administrators\n Guild1]

    >>> !admins <@10000000000000003>
    [{emoji} Administrators for Member1\n Guild1\n Guild2]

    >>> !admins 10000000000000003
    [{emoji} Administrators for Member1\n Guild1\n Guild2]
    """
    t = target.format(
        member_snowflake=DUMMY_MEMBER_SNOWFLAKE, guild_snowflake=GUILD_SNOWFLAKE
    )
    if guild is None:
        g = None
        full = f"{prefix}{command} {t}"
    else:
        g = guild.format(guild_snowflake=GUILD_SNOWFLAKE)
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
            cog = bot.get_cog("HiddenModeratorAppCommands")
            command = cog.list_administrators_app_command
            transformer = AppTarget()
            resolved_target = await transformer.transform(inx, t)
            if g:
                resolved_guild = await transformer.transform(inx, g)
            else:
                resolved_guild = None
            await command.callback(
                cog, interaction=inx, target=resolved_target, guild=resolved_guild
            )
        for kind, content in end_results:
            assert kind == "success"
