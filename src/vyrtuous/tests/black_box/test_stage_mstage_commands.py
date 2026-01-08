''' test_cstage_mstage_xstage_commands.py The purpose of this program is to black box test the stage related commands.
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
'''
from typing import Optional
from vyrtuous.inc.helpers import *
from vyrtuous.tests.black_box.make_mock_objects import *
from vyrtuous.tests.black_box.test_suite import *
from vyrtuous.enhanced_members.administrator import Administrator
from vyrtuous.enhanced_members.moderator import Moderator
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission,command,duration,ref_channel,ref_member,should_warn",
    [
        ("Administrator", "stage {voice_channel_one_id}", '1m', True, False, False),
        ("Administrator", "mstage {not_privileged_author_id}", None, False, True, False),
        ("Administrator", "stage {voice_channel_one_id}", None, True, False, False),
        ("Administrator", "stage {voice_channel_one_id}", '1h', True, False, False),
        ("Administrator", "stage {voice_channel_one_id}", None, True, False, False),
        ("Administrator", "stage {voice_channel_one_id}", '1d', True, False, False),
        ("Administrator", "stage {voice_channel_one_id}", None, True, False, False)
    ],
    indirect=['permission']
)

async def test_stage_mstage_command(
    bot,
    command: Optional[str],
    duration,
    guild,
    not_privileged_author,
    permission,
    prefix: Optional[str],
    privileged_author,
    ref_channel,
    ref_member,
    should_warn,
    text_channel,
    voice_channel_one
):
    channel_values = (voice_channel_one.mention, voice_channel_one.id)
    member_values = (not_privileged_author.mention, not_privileged_author.id)
    formatted = command.format(
            not_privileged_author_id=not_privileged_author.id,
            voice_channel_one_id=voice_channel_one.id,
            duration=duration
        )
    captured = await prepared_command_handling(author=privileged_author, bot=bot, channel=text_channel, content=formatted, guild=guild, highest_role=permission, prefix=prefix)
    message = captured[0]['message']
    message_type = captured[0]['type']
    if message.embeds:
        embed = message.embeds[0]
        content = extract_embed_text(embed)
    elif message.embed:
        content = extract_embed_text(message.embed)
    else:
        content = message.content
    if message_type == "error":
        print(f"{RED}Error:{RESET} {content}")
    # if message_type == "warning":
        # print(f"{YELLOW}Warning:{RESET} {content}")
    if message_type == "success":
        # print(f"{GREEN}Success:{RESET} {content}")
        if ref_channel:
            assert any(str(channel_value) in content for channel_value in channel_values)
        if ref_member:
            assert any(str(member_value) in content for member_value in member_values)
        assert any(emoji in content for emoji in Emojis.EMOJIS)