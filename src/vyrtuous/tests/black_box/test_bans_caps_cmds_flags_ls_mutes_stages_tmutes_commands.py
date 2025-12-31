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

def generate_cap_test_cases():
    durations = ['0', '1', '24', '48']
    moderation_types = ['ban', 'vmute', 'tmute']
    pre_cases = []
    post_cases = []
    for mod_type in moderation_types:
        for duration in durations:
            pre_cases.append((f"cap {{voice_channel_one_id}} {mod_type} {duration}", mod_type, True, False, False))
            post_cases.append((f"xcap {{voice_channel_one_id}} {mod_type}", mod_type, True, False, False))
    return pre_cases, post_cases

base_cases = [
    ("alias ban testban {voice_channel_one_id}", None, True, False, False),
    ("alias unban testunban {voice_channel_one_id}", None, True, False, False),
    ("testban {member_id}", None, False, True, False),
    ("bans {voice_channel_one_id}", None, True, False, False),
    ("bans {member_id}", None, False, True, False),
    ("bans all", None, False, False, False),
    ("testunban {member_id}", None, False, True, False),
    ("bans {voice_channel_one_id}", None, True, False, True),
    ("bans {member_id}", None, False, True, True),
    ("bans all", None, False, False, False),
    ("caps {voice_channel_one_id}", None, True, False, False),
    ("caps all", None, False, False, True),
    ("alias flag testflag {voice_channel_one_id}", None, True, False, False),
    ("alias unflag testunflag {voice_channel_one_id}", None, True, False, False),
    ("testflag {member_id}", None, False, True, False),
    ("flags {voice_channel_one_id}", None, True, False, False),
    ("flags {member_id}", None, False, True, False),
    ("flags all", None, False, False, False),
    ("testunflag {member_id}", None, False, True, False),
    ("flags {voice_channel_one_id}", None, True, False, True),
    ("flags {member_id}", None, False, True, True),
    ("flags all", None, False, False, False),
    ("alias vegan testvegan {voice_channel_one_id}", None, True, False, False),
    ("alias carnist testcarnist {voice_channel_one_id}", None, True, False, False),
    ("testvegan {member_id}", None, True, False, False),
    ("ls {voice_channel_one_id}", None, True, False, False),
    ("ls {member_id}", None, False, True, False),
    ("ls all", None, False, False, False),
    ("testcarnist {member_id}", None, True, False, False),
    ("ls {voice_channel_one_id}", None, True, False, True),
    ("ls {member_id}", None, False, True, True),
    ("ls all", None, False, False, True),
    ("alias vmute testvm {voice_channel_one_id}", None, True, False, False),
    ("alias unvmute testunvm {voice_channel_one_id}", None, True, False, False),
    ("testvm {member_id}", None, True, False, False),
    ("mutes {voice_channel_one_id}", None, True, False, False),
    ("mutes {member_id}", None, False, True, False),
    ("mutes all", None, False, False, False),
    ("testunvm {member_id}", None, True, False, False),
    ("mutes {voice_channel_one_id}", None, True, False, True),
    ("mutes {member_id}", None, False, True, True),
    ("mutes all", None, False, False, True),
    ("stages {voice_channel_one_id}", None, True, False, False),
    ("stages all", None, False, False, True),
    ("alias tmute testtm {voice_channel_one_id}", None, True, False, False),
    ("alias untmute testuntm {voice_channel_one_id}", None, True, False, False),
    ("testtm {member_id}", None, True, False, False),
    ("tmutes {voice_channel_one_id}", None, True, False, False),
    ("tmutes {member_id}", None, False, True, False),
    ("tmutes all", None, False, False, False),
    ("testuntm {member_id}", None, True, False, False),
    ("tmutes {voice_channel_one_id}", None, True, False, True),
    ("tmutes {member_id}", None, False, True, True),
    ("tmutes all", None, False, False, True),
    ("cmds {voice_channel_one_id}", None, True, False, False),
    ("cmds all", None, False, False, False),
    # ("xalias testban", None, True, False, False),
    # ("xalias testunban", None, True, False, False),
    # ("xalias testflag", None, True, False, False),
    # ("xalias testunflag", None, True, False, False),
    # ("xalias testvegan", None, True, False, False),
    # ("xalias testcarnist", None, True, False, False),
    # ("xalias testvm", None, True, False, False),
    # ("xalias testunvm", None, True, False, False),
    # ("xalias testtm", None, True, False, False),
    # ("xalias testuntm", None, True, False, False),
    ("cmds {voice_channel_one_id}", None, True, False, True),
    ("cmds all", None, False, False, True),
]

pre_caps, post_caps = generate_cap_test_cases()
all_cases = pre_caps + base_cases + post_caps

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,moderation_type,channel_ref,member_ref,should_warn",
    all_cases
)
async def test_bans_caps_cmds_flags_ls_mutes_stages_tmutes_commands(bot, voice_channel_one, guild, not_privileged_author, privileged_author, prefix: Optional[str], command: Optional[str], channel_ref, member_ref, should_warn, moderation_type):    
    moderator = Moderator(channel_snowflake=voice_channel_one.id, guild_snowflake=guild.id, member_snowflake=privileged_author.id)
    await moderator.grant()
    try:
        match moderation_type:
            case "tmute":
                moderation_type = "text_mute"
            case "vmute":
                moderation_type = "voice_mute"
            case "untmute":
                moderation_type = "untext_mute"
            case "unvmute":
                moderation_type = "unvoice_mute"
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        member_value = not_privileged_author.name
        voice_channel_one.messages.clear() 
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id,
            member_id=not_privileged_author.id
        )
        bot.wait_for = mock_wait_for
        captured = await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
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
            # print(f"{YELLOW}Warning:{RESET} {content}")
            if should_warn:
                assert True
        if message_type == "success":
            # print(f"{GREEN}Success:{RESET} {content}")
            assert any(emoji in content for emoji in Emojis.EMOJIS) 
            if member_ref:
                assert member_value in content
            if channel_ref:
                assert channel_value in content
            if moderation_type:
                assert moderation_type in content 
    finally:
        await moderator.revoke()