"""test_del_commands.py The purpose of this program is to black box test the message delete command. FIX NO MESSAGE FOUND ERROR

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
from vyrtuous.inc.helpers import RESET, YELLOW, RED, GREEN
from vyrtuous.tests.black_box.test_suite import (
    extract_embed_text,
    prepared_command_handling,
    MESSAGE_ID,
)
from vyrtuous.utils.emojis import EMOJIS
import pytest


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission,command,needs_message_id, ref_channel,ref_guild,ref_member,should_warn",
    [
        ("Developer", "ping", False, False, False, False, False),
        ("Moderator", "del {message_id}", True, False, False, False, False),
    ],
    indirect=["permission"],
)
async def test_del_command(
    bot,
    command: Optional[str],
    needs_message_id,
    not_privileged_author,
    permission,
    prefix: Optional[str],
    privileged_author,
    ref_channel,
    ref_guild,
    ref_member,
    text_channel,
    voice_channel_one,
    guild,
):
    channel_values = (voice_channel_one.mention, voice_channel_one.id)
    guild_values = (guild.name, guild.id)
    member_values = (not_privileged_author.mention, not_privileged_author.id)

    if needs_message_id:
        formatted = command.format(message_id=MESSAGE_ID)
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
    else:
        formatted = command
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
