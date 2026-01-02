''' test_rmv_command.py The purpose of this program is to black box test the room move command.
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
from vyrtuous.tests.black_box.test_suite import *
from vyrtuous.utils.administrator import Administrator
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref_one,channel_ref_two",
    [
        ("rmv {source_id} {target_id}", True, True)
    ]
)

async def test_rmv_command(bot, voice_channel_one, voice_channel_two, guild, privileged_author, prefix: Optional[str], role, command: Optional[str], channel_ref_one, channel_ref_two):
    administrator = Administrator(guild_snowflake=guild.id, member_snowflake=privileged_author.id, role_snowflakes=[role.id])
    await administrator.grant()
    try:
        voice_channel_one.messages.clear() 
        channel_value_one = voice_channel_one.mention if channel_ref_one else voice_channel_one.name
        channel_value_two = voice_channel_two.mention if channel_ref_two else voice_channel_two.name
        source_id = voice_channel_one.id
        target_id = voice_channel_two.id
        formatted = command.format(
            source_id=source_id,
            target_id=target_id
        )
        captured = await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
        message = captured['message']
        message_type = captured['type']
        if isinstance(message, discord.Embed):
            content = extract_embed_text(message)
        elif isinstance(message, discord.File):
            content = message.filename
        else:
            content = message
        if message_type == "error":
            print(f"{RED}Error:{RESET} {content}")
        if message_type == "warning":
            print(f"{YELLOW}Warning:{RESET} {content}")
        if message_type == "success":
            # print(f"{GREEN}Success:{RESET} {content}")
            assert any(emoji in content for emoji in Emojis.EMOJIS)
            assert channel_value_one in content
            assert channel_value_two in content
    finally:
        await administrator.revoke()