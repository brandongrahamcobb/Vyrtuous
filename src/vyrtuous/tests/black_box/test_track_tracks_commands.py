''' test_log_logs_mlog_commands.py The purpose of this program is to black box test the logging related commands.

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
from vyrtuous.tests.black_box.test_suite import *
from vyrtuous.enhanced_member.administrator import Administrator
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission,command,action,target_type,ref_channel,ref_member,ref_text,should_warn",
    [
        ('Administrator', "track", "create", "all", False, False, True, False),
        ('Administrator', "track", None, None, False, False, True, False),
        ('Administrator', "track", None, None, False, False, True, False),
        ('Developer', "tracks", "all", None, False, False, False, False),
        ('Administrator', "tracks", None, None, False, False, False, False),
        ('Administrator', "track", "modify", "all", False, False, True, False),
        ('Administrator', "track", "delete", "all", False, False, True, False),
        ('Administrator', "track", "create", "channel", True, False, True, False),
        ('Administrator', "track", "modify", "channel", True, False, True, False),
        ('Administrator', "track", "delete", None, False, False, True, False),
        ('Administrator', "track", "create", "member", False, True, True, False),
        ('Administrator', "track", "modify", "member", False, True, True, False),
        ('Administrator', "track", "delete", None, False, False, True, False)
    ],
    indirect=['permission']
)

async def test_log_logs_mlog_command(
    action: Optional[str],
    bot,
    command: Optional[str],
    guild,
    not_privileged_author,
    permission,
    prefix: Optional[str],
    privileged_author,
    ref_channel,
    ref_member,
    ref_text,
    role,
    should_warn,
    target_type: Optional[str],
    text_channel,
    voice_channel_one
):
    text_channel_values = (text_channel.mention, text_channel.id)
    voice_channel_values = (voice_channel_one.mention, voice_channel_one.id)
    privileged_author_values = (privileged_author.mention, privileged_author.id)
    not_privileged_author_values = (not_privileged_author.mention, not_privileged_author.id)
    if command == "track":
        match target_type:
            case "all":
                formatted = f"{command} {text_channel.id} {action} {target_type}".strip()
            case "channel":
                formatted = f"{command} {text_channel.id} {action} {target_type} {voice_channel_one.id}".strip()
            case "member":
                formatted = f"{command} {text_channel.id} {action} {target_type} {not_privileged_author.id}".strip()
            case None:
                formatted = f"{command} {text_channel.id}".strip()
    if command == "tracks" and action:
        formatted = f"{command} {action}".strip()
    elif command == "tracks":
        formatted = f"{command}".strip()
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
    if message_type == "warning":
        print(f"{YELLOW}Warning:{RESET} {content}")
    if message_type == "success":
        # print(f"{GREEN}Success:{RESET} {content}")
        assert any(emoji in content for emoji in Emojis.EMOJIS)
        if ref_text:
            assert any(str(text_channel_value) in content for text_channel_value in text_channel_values)
        if ref_channel:
            assert any(str(voice_channel_value) in content for voice_channel_value in voice_channel_values)
        if ref_member:
            assert any(str(not_privileged_author_value) in content for not_privileged_author_value in not_privileged_author_values)