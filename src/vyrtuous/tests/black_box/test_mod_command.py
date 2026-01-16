"""test_mod_command.py The purpose of this program is to black box test the moderator toggle command.

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

from typing import Optional

import pytest

from vyrtuous.tests.black_box.test_suite import (
    extract_embed_text,
    prepared_command_handling,
    RESET,
    YELLOW,
    RED,
    GREEN,
)
from vyrtuous.utils.emojis import EMOJIS


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission,command,ref_channel,ref_guild, ref_member,should_warn",
    [
        (
            "Coordinator",
            "mod {not_privileged_author_id} {voice_channel_one_id}",
            True,
            False,
            True,
            False,
        ),
        (
            "Coordinator",
            "mod {not_privileged_author_id} {voice_channel_one_id}",
            True,
            False,
            True,
            False,
        ),
    ],
    indirect=["permission"],
)
async def test_mod_command(
    bot,
    command: Optional[str],
    guild,
    not_privileged_author,
    permission,
    prefix: Optional[str],
    privileged_author,
    ref_channel,
    ref_guild,
    ref_member,
    text_channel,
    voice_channel_one,
):
    channel_values = (voice_channel_one.mention, voice_channel_one.id)
    guild_values = (guild.name, guild.id)
    member_values = (not_privileged_author.mention, not_privileged_author.id)
    formatted = command.format(
        voice_channel_one_id=voice_channel_one.id,
        not_privileged_author_id=not_privileged_author.id,
    )
    captured = await prepared_command_handling(
        author=privileged_author,
        bot=bot,
        channel=text_channel,
        content=formatted,
        guild=guild,
        highest_role=permission,
        prefix=prefix,
    )
    message = captured[0]["message"]
    message_type = captured[0]["type"]
    if message.embeds:
        embed = message.embeds[0]
        content = extract_embed_text(embed)
    elif message.embed:
        content = extract_embed_text(message.embed)
    else:
        content = message.content
    if message_type == "error":
        print(f"{RED}Error:{RESET} {content}")
    if message_type == "warning":
        print(f"{YELLOW}Warning:{RESET} {content}")
    if message_type == "success":
        print(f"{GREEN}Success:{RESET} {content}")
        if ref_channel:
            assert any(
                str(channel_value) in content for channel_value in channel_values
            )
        if ref_guild:
            assert any(str(guild_value) in content for guild_value in guild_values)
        if ref_member:
            assert any(str(member_value) in content for member_value in member_values)
        assert any(emoji in content for emoji in EMOJIS)
