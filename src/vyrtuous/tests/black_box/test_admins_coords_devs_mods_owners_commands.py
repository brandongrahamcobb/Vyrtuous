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
from vyrtuous.utils.moderator import Moderator
from vyrtuous.tests.black_box.test_suite import bot, config, guild, not_privileged_author, prepared_command_handling, prefix, privileged_author, voice_channel_one
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref,member_ref",
    [
        ("admins", False, False),
        ("coords {voice_channel_one_id}", True, False),
        ("coords {member_id}", False, True),
        ("coords all", False, False),
        ("devs", False, False),
        ("devs {member_id}", False, True),
        ("devs all", False, False),
        ("mods {voice_channel_one_id}", True, False),
        ("mods {member_id}", False, True),
        ("mods all", False, False),
        ("owners {member_id}", False, True),
        ("owners {member_id}", False, True),
        ("owners all", False, False)
    ]
)

async def test_admins_coords_devs_mods_owners_commands(bot, voice_channel_one, guild, not_privileged_author, privileged_author, prefix: Optional[str], command: Optional[str], channel_ref, member_ref):    
    moderator = Moderator(channel_snowflake=voice_channel_one.id, guild_snowflake=guild.id, member_snowflake=privileged_author.id)
    await moderator.grant()
    try:
        voice_channel_one.messages.clear() 
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id,
            member_id=not_privileged_author.id
        )
        bot.wait_for = mock_wait_for
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="EveryoneCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.everyone_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]
        member_values = (not_privileged_author.mention, not_privileged_author.name)
        if response["embed"]:
            assert any(emoji in response["embed"].title for emoji in Emojis.EMOJIS) 
        else:
            assert any(emoji in response["content"] for emoji in Emojis.EMOJIS) 
    finally:
        await moderator.revoke()