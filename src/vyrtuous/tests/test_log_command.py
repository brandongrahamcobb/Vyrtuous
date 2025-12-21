''' test_log_command.py The purpose of this program is to black box test the log command.
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
from vyrtuous.utils.setup_logging import logger, setup_logging

import asyncio
import discord
import pytest
import pytest_asyncio
import uuid

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref",
    [
        ("log", True),
        ("log", True),
        ("logs all", False),
        ("logs", True)
    ]
)

async def test_log_command(bot, voice_channel_one, text_channel, guild, privileged_author, not_privileged_author, prefix: Optional[str], command: Optional[str], channel_ref):
    await admin_initiation(guild.id, privileged_author.id)
    try:
        if not channel_ref:
            formatted = f"{command}".strip()
        else:
            formatted = f"{command} {text_channel.id}".strip()
        await prepared_command_handling(author=privileged_author, bot=bot, content=formatted, data=None, guild=guild, prefix=prefix, text_channel=text_channel, voice_channel=voice_channel_one)
        response = voice_channel_one.messages[0]
        print(response)
        if not channel_ref:
            embed = response["embed"]
            assert any(emoji in embed for emoji in Emojis.EMOJIS)
        else:
            channel_value = text_channel.mention if channel_ref else text_channel.name
            assert any(emoji in response for emoji in Emojis.EMOJIS)
            assert any(val in response for val in [channel_value])
        voice_channel_one.messages.clear() 
    finally:
        await admin_cleanup(guild.id, privileged_author.id)