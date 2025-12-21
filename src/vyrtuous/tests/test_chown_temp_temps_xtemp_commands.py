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
from discord.ext.commands import Context, view as cmd_view
from types import SimpleNamespace
from typing import Optional
from unittest.mock import PropertyMock, patch
from vyrtuous.inc.helpers import *
from vyrtuous.tests.make_mock_objects import *
from vyrtuous.tests.test_admin_helpers import *
from vyrtuous.tests.test_suite import *
from vyrtuous.utils.emojis import Emojis
import asyncio
import discord
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref,member_ref",
    [
        ("temp {voice_channel_one_id} {member_id}", "channel", "member"),
        ("temps all", None, None),
        ("temps {voice_channel_one_id}", None, None),
        ("chown {voice_channel_one_id} {member_id}", "channel", "member"),
        ("temp {voice_channel_one_id}", "channel", None),
        ("temp {channel_mention} {member_mention}", "channel", "member"),
        ("chown {channel_mention} {member_mention}", "channel", "member"),
        ("temp {voice_channel_one_id}", "channel", None)
    ]
)

async def test_chown_temp_xtemp_commands(bot, voice_channel_one, guild, privileged_author, prefix: Optional[str], command: Optional[str], channel_ref, member_ref):    
    await admin_initiation(guild.id, privileged_author.id)
    try:
        voice_channel_one.messages.clear() 
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id,
            channel_mention=voice_channel_one.mention,
            member_id=privileged_author.id,
            member_mention=privileged_author.mention,
        )
        bot.wait_for = mock_wait_for
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, content=formatted, guild=guild, prefix=prefix)
        response = voice_channel_one.messages[0]
        print(response)
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        member_value = privileged_author.mention if member_ref else privileged_author.name
        if "temps" in command:
            assert any(emoji in response["embed"].title for emoji in Emojis.EMOJIS) 
        else:
            assert any(emoji in response["content"] for emoji in Emojis.EMOJIS) 
            assert any(val in response["content"] for val in [channel_value])
        if member_ref:
            assert any(val in response["content"] for val in [member_value])
    finally:
        await admin_cleanup(guild.id, privileged_author.id)