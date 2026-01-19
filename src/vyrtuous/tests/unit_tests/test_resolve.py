"""test_resolve.py The purpose of this program is to test the DiscordObject resolution module.
Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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

import discord
import pytest

from vyrtuous.inc.helpers import (
    GUILD_ID,
    MESSAGE_ID,
    PRIVILEGED_AUTHOR_ID,
    ROLE_ID,
    TEXT_CHANNEL_ID,
)
from vyrtuous.tests.black_box.make_mock_objects import create_message
from vyrtuous.tests.black_box.test_suite import (
    bot,
    guild,
    prepare_context,
    privileged_author,
    role,
    text_channel,
)
from vyrtuous.service.logging_service import logger
from vyrtuous.service.resolution.discord_object_service import DiscordObject


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "discord_object_type,dictionary",
    [
        # (discord.abc.GuildChannel, {"channel_mention": f'<#{TEXT_CHANNEL_ID}>'}),
        # (discord.abc.GuildChannel, {"channel_snowflake": TEXT_CHANNEL_ID}),
        # (discord.Guild, {"guild_snowflake": GUILD_ID}),
        # (discord.Member, {"member_mention": f'<@{PRIVILEGED_AUTHOR_ID}>'}),
        # (discord.Member, {"member_snowflake": PRIVILEGED_AUTHOR_ID}),
        (discord.Role, {"role_mention": f"<@&{ROLE_ID}>"}),
        (discord.Role, {"role_snowflake": ROLE_ID}),
    ],
)
async def test_resolve_target(
    bot, dictionary, discord_object_type, guild, privileged_author, role, text_channel
):
    bot.get_guild = lambda guild_id: guild if guild_id == GUILD_ID else None
    # guild.get_role = lambda role_id: role if role_id == ROLE_ID else None

    message = create_message(
        allowed_mentions=True,
        author=privileged_author,
        bot=bot,
        content="",
        guild=guild,
        channel=text_channel,
        id=MESSAGE_ID,
    )
    for key, value in dictionary.items():
        ctx = await prepare_context(bot, message, "!")
        do = DiscordObject(source=ctx)
        obj = await do.determine_from_target(ctx, value)
        assert isinstance(obj, discord_object_type)
