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
        ("bans {voice_channel_one_id}", True, False),
        ("bans {member_id}", False, True),
        ("bans all", False, False),
        ("caps {voice_channel_one_id}", True, False),
        ("caps all", False, False),
        ("cmds {voice_channel_one_id} {member_id}", True, True),
        ("cmds all", False, False),
        ("flags {voice_channel_one_id}", True, False),
        ("flags {member_id}", False, True),
        ("flags all", False, False),
        ("ls {voice_channel_one_id}", True, False),
        ("ls {member_id}", False, True),
        ("ls all", False, False),
        ("mutes {voice_channel_one_id}", True, False),
        ("mutes {member_id}", False, True),
        ("mutes all", False, False),
        ("stages {voice_channel_one_id}", True, False),
        ("stages all", False, False),
        ("tmutes {voice_channel_one_id}", True, False),
        ("tmutes {member_id}", False, True),
        ("tmutes all", False, False),
    ]
)

async def test_bans_caps_cmds_flags_ls_mutes_stages_tmutes_commands(bot, voice_channel_one, guild, not_privileged_author, privileged_author, prefix: Optional[str], command: Optional[str], channel_ref, member_ref):    
    await Moderator.grant(channel_id=voice_channel_one.id, guild_id=guild.id, member_id=privileged_author.id)
    try:
        voice_channel_one.messages.clear() 
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id,
            member_id=not_privileged_author.id
        )
        bot.wait_for = mock_wait_for
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="ModeratorCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.moderator_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]
        member_values = (not_privileged_author.mention, not_privileged_author.name)
        if response["embed"]:
            assert any(emoji in response["embed"].title for emoji in Emojis.EMOJIS) 
        else:
            assert any(emoji in response["content"] for emoji in Emojis.EMOJIS) 
        if member_ref and response["embed"]:
            assert any(val in response["embed"].title for val in member_values)
    finally:
        await Moderator.revoke(channel_id=voice_channel_one.id, member_id=privileged_author.id)