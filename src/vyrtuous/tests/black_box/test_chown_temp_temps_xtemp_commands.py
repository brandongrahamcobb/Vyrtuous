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
from vyrtuous.tests.black_box.test_admin_helpers import admin_cleanup, admin_initiation
from vyrtuous.tests.black_box.test_suite import bot, config, guild, not_privileged_author, prepared_command_handling, prefix, privileged_author, voice_channel_one
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref,member_ref",
    [
        ("temp {voice_channel_one_id} {member_id}", True, True),
        ("temps all", False, False),
        ("temps {voice_channel_one_id}", True, False),
        ("chown {voice_channel_one_id} {member_id}", True, True),
        ("temp {voice_channel_one_id}", True, False),
        ("temp {channel_mention} {member_mention}", True, True),
        ("chown {channel_mention} {member_mention}", True, True),
        ("migrate \"{channel_name}\" {voice_channel_one_id}", True, False),
        ("temp {voice_channel_one_id}", True, False)
    ]
)

async def test_chown_temp_xtemp_commands(bot, voice_channel_one, guild, not_privileged_author, privileged_author, prefix: Optional[str], command: Optional[str], channel_ref, member_ref):    
    await admin_initiation(guild.id, privileged_author.id)
    try:
        voice_channel_one.messages.clear() 
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id,
            channel_name=voice_channel_one.name,
            channel_mention=voice_channel_one.mention,
            member_id=not_privileged_author.id,
            member_mention=not_privileged_author.mention
        )
        bot.wait_for = mock_wait_for
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]
        channel_values = (voice_channel_one.mention, voice_channel_one.name)
        member_values = (not_privileged_author.mention, not_privileged_author.name)
        if "temps" in command:
            assert any(emoji in response["embed"].title for emoji in Emojis.EMOJIS) 
        else:
            assert any(emoji in response["content"] for emoji in Emojis.EMOJIS) 
            assert any(val in response["content"] for val in channel_values)
        if member_ref:
            assert any(val in response["content"] for val in member_values)
    finally:
        await admin_cleanup(guild.id, privileged_author.id)