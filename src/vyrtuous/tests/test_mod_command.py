''' test_coord_command.py The purpose of this program is to black box test the coord command.
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
from vyrtuous.tests.test_coord_helpers import coord_cleanup, coord_initiation
from vyrtuous.tests.test_suite import bot, config, guild, not_privileged_author, prepared_command_handling, prefix, privileged_author, voice_channel_one
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref,member_ref",
    [
        ("mod {not_privileged_author_id} {voice_channel_one_id}", True, True),
        ("mod {not_privileged_author_id} {voice_channel_one_id}", True, True)
    ]
)

async def test_mod_command(bot, voice_channel_one, guild, privileged_author, not_privileged_author, prefix: Optional[str], command: Optional[str], channel_ref, member_ref):
    await coord_initiation(voice_channel_one.id, guild.id, privileged_author.id)
    try:
        voice_channel_one.messages.clear() 
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id,
            not_privileged_author_id=not_privileged_author.id
        )
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="CoordinatorCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.coordinator_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]["content"]
        assert any(emoji in response for emoji in Emojis.EMOJIS)
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        member_value = not_privileged_author.mention if member_ref else not_privileged_author.name
        assert any(val in response for val in [channel_value])
        assert any(val in response for val in [member_value])
    finally:
        await coord_cleanup(voice_channel_one.id, guild.id, privileged_author.id)