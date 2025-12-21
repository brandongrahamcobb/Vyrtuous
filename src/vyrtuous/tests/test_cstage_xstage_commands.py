''' test_cstage_xstage_commands.py The purpose of this program is to black box test the stages commands.
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
import discord
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,duration,channel_ref",
    [
        ("cstage {voice_channel_one_id}", '1m', True),
        ("xstage {voice_channel_one_id}", None, True),
        ("cstage {voice_channel_one_id}", '1h', True),
        ("xstage {voice_channel_one_id}", None, True),
        ("cstage {voice_channel_one_id}", '1d', True),
        ("xstage {voice_channel_one_id}", None, True)
    ]
)

async def test_cstage_xstage_command(bot, voice_channel_one, guild, privileged_author, prefix: Optional[str], command: Optional[str], duration, channel_ref):
    await admin_initiation(guild.id, privileged_author.id)
    try:
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id,
            duration=duration
        )
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, content=formatted, guild=guild, prefix=prefix)
        response = voice_channel_one.messages[0]["content"]
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        assert any(emoji in response for emoji in Emojis.EMOJIS)
        assert any(val in response for val in [channel_value])
        voice_channel_one.messages.clear() 
    finally:
        await admin_cleanup(guild.id, privileged_author.id)