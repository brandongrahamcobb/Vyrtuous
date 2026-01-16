"""test_rmv_command.py The purpose of this program is to black box test the room move command.
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
    RESET, YELLOW, RED, GREEN,
)
from vyrtuous.utils.emojis import EMOJIS


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission,command,ref_channel_one,ref_channel_two,should_warn",
    [("Administrator", "rmv {source_id} {target_id}", True, True, False)],
    indirect=["permission"],
)
async def test_rmv_command(
    bot,
    command: Optional[str],
    guild,
    permission,
    prefix: Optional[str],
    privileged_author,
    ref_channel_one,
    ref_channel_two,
    text_channel,
    voice_channel_one,
    voice_channel_two,
):
    channel_one_values = (voice_channel_one.mention, voice_channel_one.id)
    channel_two_values = (voice_channel_two.mention, voice_channel_two.id)
    source_id = voice_channel_one.id
    target_id = voice_channel_two.id
    formatted = command.format(source_id=source_id, target_id=target_id)
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
        if ref_channel_one:
            assert any(
                str(channel_one_value) in content
                for channel_one_value in channel_one_values
            )
        if ref_channel_two:
            assert any(
                str(channel_two_value) in content
                for channel_two_value in channel_two_values
            )
        assert any(emoji in content for emoji in EMOJIS)
