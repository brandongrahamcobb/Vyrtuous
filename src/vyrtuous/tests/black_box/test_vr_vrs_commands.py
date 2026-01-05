''' test_chown_temp_temps_xtemp_commands.py The purpose of this program is to black box test the temporary room commands.
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
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref,should_warn",
    [
        ("vr {voice_channel_one_id}", True, False),
        ("vrs", False, False),
        ("vr {voice_channel_one_id}", True, False),
        ("vrs", False, True)
    ]
)

async def test_chown_temp_xtemp_commands(bot, voice_channel_one, guild, not_privileged_author, privileged_author, prefix: Optional[str], role, command: Optional[str], channel_ref, should_warn):    
    administrator = Administrator(guild_snowflake=guild.id, member_snowflake=privileged_author.id, role_snowflakes=[role.id])
    await administrator.grant()
    try:
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        voice_channel_one.messages.clear() 
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id
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
            assert should_warn
        if message_type == "success":
            print(f"{GREEN}Success:{RESET} {content}")
            if channel_ref:
                assert channel_value in content
            assert any(emoji in content for emoji in Emojis.EMOJIS)
    finally:
        await administrator.revoke()