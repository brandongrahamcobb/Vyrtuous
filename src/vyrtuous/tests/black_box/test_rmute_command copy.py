''' test_smute_command.py The purpose of this program is to black box test the smute command.
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
from vyrtuous.tests.black_box.test_coord_helpers import coord_cleanup, coord_initiation
from vyrtuous.tests.black_box.test_suite import bot, config, guild, not_privileged_author, prepared_command_handling, prefix, privileged_author, voice_channel_one
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,member_ref",
    [
        ("rmute", "True"),
        ("rmute", "True")
    ]
)

async def test_rmute_command(bot, voice_channel_one, guild, privileged_author, not_privileged_author, prefix: Optional[str], command: Optional[str], member_ref):
    await coord_initiation(voice_channel_one.id, guild.id, privileged_author.id)
    try:
        voice_channel_one.messages.clear() 
        formatted = f"{command} {voice_channel_one.id}"
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="CoordinatorCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.coordinator_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]["content"]
        assert any(emoji in response for emoji in Emojis.EMOJIS)
    finally:
        await coord_cleanup(voice_channel_one.id, guild.id, privileged_author.id)