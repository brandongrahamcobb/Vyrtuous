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
from vyrtuous.tests.black_box.test_mod_helpers import mod_cleanup, mod_initiation
from vyrtuous.tests.black_box.test_suite import bot, config, guild, not_privileged_author, prepared_command_handling, prefix, privileged_author, voice_channel_one
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        ("survey {voice_channel_one_id}")
    ]
)

async def test_survey_command(bot, voice_channel_one, guild, not_privileged_author, privileged_author, prefix: Optional[str], command: Optional[str]):    
    await mod_initiation(voice_channel_one.id, guild.id, privileged_author.id)
    try:
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id
        )
        voice_channel_one.messages.clear() 
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="EveryoneCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.everyone_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]
        assert any(emoji in response["content"] for emoji in Emojis.EMOJIS)
    finally:
        await mod_cleanup(voice_channel_one.id, guild.id, privileged_author.id)