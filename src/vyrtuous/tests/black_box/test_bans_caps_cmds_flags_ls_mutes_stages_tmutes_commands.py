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
from vyrtuous.tests.black_box.test_suite import *
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref,member_ref,should_warn",
    [
        ("alias ban testban {voice_channel_one_id}", True, False, False),
        ("alias unban testunban {voice_channel_one_id}", True, False, False),
        ("testban {member_id}", False, True, False),
        # ("bans {voice_channel_one_id}", True, False, False),
        # ("bans {member_id}", False, True, False),
        # ("bans all", False, False, False),
        # ("testunban {member_id}", False, True, False),
        # ("bans {voice_channel_one_id}", True, False, True),
        # ("bans {member_id}", False, True, True),
        # ("bans all", False, False, False),
        # ("xalias testban", True, False, False),
        # ("xalias testunban", True, False, False),
        # ("caps {voice_channel_one_id}", True, False, False),
        # ("caps all", False, False, True),
        # ("cmds {voice_channel_one_id} {member_id}", True, True, False),
        # ("cmds all", False, False, True),
        # ("flags {voice_channel_one_id}", True, False, False),
        # ("flags {member_id}", False, True, False),
        # ("flags all", False, False, True),
        # ("ls {voice_channel_one_id}", True, False, False),
        # ("ls {member_id}", False, True, False),
        # ("ls all", False, False, True),
        # ("mutes {voice_channel_one_id}", True, False, False),
        # ("mutes {member_id}", False, True, False),
        # ("mutes all", False, False, True),
        # ("stages {voice_channel_one_id}", True, False, False),
        # ("stages all", False, False, True),
        # ("tmutes {voice_channel_one_id}", True, False, False),
        # ("tmutes {member_id}", False, True, False),
        # ("tmutes all", False, False, True),
    ]
)

async def test_bans_caps_cmds_flags_ls_mutes_stages_tmutes_commands(bot, voice_channel_one, guild, not_privileged_author, privileged_author, prefix: Optional[str], command: Optional[str], channel_ref, member_ref, should_warn):    
    moderator = Moderator(channel_snowflake=voice_channel_one.id, guild_snowflake=guild.id, member_snowflake=privileged_author.id)
    await moderator.grant()
    try:
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        member_value = not_privileged_author.mention if member_ref else not_privileged_author.name
        voice_channel_one.messages.clear() 
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id,
            member_id=not_privileged_author.id
        )
        bot.wait_for = mock_wait_for
        captured = await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="EventListeners", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
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
        if message_type == "success":
            print(f"{GREEN}Success:{RESET} {content}")
            assert any(emoji in content for emoji in Emojis.EMOJIS) 
            if member_ref:
                assert member_value in content
            if channel_ref:
                assert channel_value in content
    finally:
        await moderator.revoke()