''' test_backup_command.py The purpose of this program is to black box test the backup command.
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
from vyrtuous.enhanced_member.developer import Developer
from vyrtuous.tests.black_box.test_suite import *
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission,command,should_warn",
    [
        ('Developer', "backup", False)
    ],
    indirect=['permission']
)

async def test_backup_command(
    bot,
    command: Optional[str],
    guild,
    permission,
    prefix: Optional[str],
    privileged_author,
    should_warn,
    text_channel,
    voice_channel_one
):
    captured = await prepared_command_handling(author=privileged_author, bot=bot, channel=text_channel, content=command, guild=guild, highest_role=permission, prefix=prefix)
    message = captured[0]['message']
    if message.embeds:
        embed = message.embeds[0]
        content = extract_embed_text(embed)
    elif message.embed:
        content = extract_embed_text(message.embed)
    elif message.content:
        content = message.content
    elif message.file:
        content = message.file.filename
    message_type = captured[0]['type']
    if message_type == "error":
        print(f"{RED}Error:{RESET} {content}")
    if message_type == "warning":
        print(f"{YELLOW}Warning:{RESET} {content}")
        if should_warn:
            assert True
    if message_type == "success":
        # print(f"{GREEN}Success:{RESET} {content}")
        assert content.startswith("backup")