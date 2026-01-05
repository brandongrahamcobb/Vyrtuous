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
from vyrtuous.enhanced_members.moderator import Moderator
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
            pre_cases.append((f"Administrator", f"cap {{voice_channel_one_id}} {mod_type} {duration}", mod_type, True, False, False, False))
            post_cases.append((f"Administrator", f"cap {{voice_channel_one_id}} {mod_type} {duration}", mod_type, True, False, False, False))
    return pre_cases, post_cases

base_cases = [
    ("Administrator", "alias ban testban {voice_channel_one_id}", None, True, False, False, False),
    ("Administrator", "alias unban testunban {voice_channel_one_id}", None, True, False, False, False),
    ("Moderator", "testban {member_id}", None, True, False, True, False),
    ("Moderator", "bans {voice_channel_one_id}", None, False, True, False, False),
    ("Moderator", "bans {member_id}", None, False, False, True, False),
    ("Administrator", "bans {guild_id}", None, False, True, False, False),
    ("Developer", "bans all", None, False, False, False, False),
    ("Moderator", "testunban {member_id}", None, True, False, True, False),
    ("Moderator", "bans {voice_channel_one_id}", None, False, True, False, True),
    ("Moderator", "bans {member_id}", None, False, False, True, True),
    ("Administrator", "bans {guild_id}", None, False, True, False, True),
    ("Developer", "bans all", None, False, False, False, True),
    ("Moderator", "caps {voice_channel_one_id}", None, True, False, False, False),
    ("Administrator", "caps {guild_id}", None, False, True, False, False),
    ("Developer", "caps all", None, False, False, False, False),
    ("Administrator", "alias flag testflag {voice_channel_one_id}", None, True, False, False, False),
    ("Administrator", "alias unflag testunflag {voice_channel_one_id}", None, True, False, False, False),
    ("Moderator", "testflag {member_id}", None, True, False, True, False),
    ("Moderator", "flags {voice_channel_one_id}", None, False, True, False, False),
    ("Moderator", "flags {member_id}", None, False, False, True, False),
    ("Administrator", "flags {guild_id}", None, False, True, False, False),
    ("Developer", "flags all", None, False, False, False, False),
    ("Moderator", "testunflag {member_id}", None, True, False, True, False),
    ("Moderator", "flags {voice_channel_one_id}", None, False, True, False, True),
    ("Moderator", "flags {member_id}", None, False, False, True, True),
    ("Administrator", "flags {guild_id}", None, False, True, False, True),
    ("Developer", "flags all", None, False, False, False, True),
    ("Administrator", "alias vegan testvegan {voice_channel_one_id}", None, True, False, False, False),
    ("Administrator", "alias carnist testcarnist {voice_channel_one_id}", None, True, False, False, False),
    ("Moderator", "testvegan {member_id}", None, True, False, True, True),
    ("Moderator", "ls {voice_channel_one_id}", None, True, False, False, False),
    ("Moderator", "ls {member_id}", None, False, False, True, False),
    ("Administrator", "ls {guild_id}", None, False, True, False, False),
    ("Developer", "ls all", None, False, False, False, False),
    ("Moderator", "testcarnist {member_id}", None, True, False, True, True),
    ("Moderator", "ls {voice_channel_one_id}", None, True, False, False, True),
    ("Moderator", "ls {member_id}", None, False, False, True, True),
    ("Administrator", "ls {guild_id}", None, False, True, False, True),
    ("Developer", "ls all", None, False, False, False, True),
    ("Administrator", "alias vmute testvm {voice_channel_one_id}", None, True, False, False, False),
    ("Administrator", "alias unvmute testunvm {voice_channel_one_id}", None, True, False, False, False),
    ("Moderator", "testvm {member_id}", None, True, False, False, False),
    ("Moderator", "mutes {voice_channel_one_id}", None, False, True, False, False),
    ("Moderator", "mutes {member_id}", None, False, False, True, False),
    ("Administrator", "mutes {guild_id}", None, False, True, False, False),
    ("Developer", "mutes all", None, False, False, False, False),
    ("Moderator", "testunvm {member_id}", None, True, False, True, False),
    ("Moderator", "mutes {voice_channel_one_id}", None, False, True, False, True),
    ("Moderator", "mutes {member_id}", None, False, False, True, True),
    ("Administrator", "mutes {guild_id}", None, False, True, False, True),
    ("Developer", "mutes all", None, False, False, False, True),
    ("Moderator", "stages {voice_channel_one_id}", None, False, True, False, True),
    ("Administrator", "stages {guild_id}", None, False, True, False, True),
    ("Developer", "stages all", None, False, False, False, True),
    ("Administrator", "alias tmute testtm {voice_channel_one_id}", None, True, False, False, False),
    ("Administrator", "alias untmute testuntm {voice_channel_one_id}", None, True, False, False, False),
    ("Moderator", "testtm {member_id}", None, True, False, True, False),
    ("Moderator", "tmutes {voice_channel_one_id}", None, False, False, False, False),
    ("Moderator", "tmutes {member_id}", None, False, False, True, False),
    ("Administrator", "tmutes {guild_id}", None, False, True, False, False),
    ("Developer", "tmutes all", None, False, False, False, False),
    ("Moderator", "testuntm {member_id}", None, True, False, True, False),
    ("Moderator", "tmutes {voice_channel_one_id}", None, False, True, False, True),
    ("Moderator", "tmutes {member_id}", None, False, False, True, True),
    ("Administrator", "tmutes {guild_id}", None, False, True, False, True),
    ("Developer", "tmutes all", None, False, False, False, True),
    ("Moderator", "cmds {voice_channel_one_id}", None, False, True, False, False),
    ("Developer", "cmds {guild_id}", None, False, True, False, False),
    ("Developer", "cmds all", None, False, False, False, False),
    ("Administrator", "xalias testban", None, True, False, False, False),
    ("Administrator", "xalias testunban", None, True, False, False, False),
    ("Administrator", "xalias testflag", None, True, False, False, False),
    ("Administrator", "xalias testunflag", None, True, False, False, False),
    ("Administrator", "xalias testvegan", None, True, False, False, False),
    ("Administrator", "xalias testcarnist", None, True, False, False, False),
    ("Administrator", "xalias testvm", None, True, False, False, False),
    ("Administrator", "xalias testunvm", None, True, False, False, False),
    ("Administrator", "xalias testtm", None, True, False, False, False),
    ("Administrator", "xalias testuntm", None, True, False, False, False),
    ("Moderator", "cmds {voice_channel_one_id}", None, False, True, False, True),
    ("Administrator", "cmds {guild_id}", None, False, True, False, True),
    ("Developer", "cmds all", None, False, False, False, True)
]

pre_caps, post_caps = generate_cap_test_cases()
all_cases = pre_caps + base_cases + post_caps

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission,command,moderation_type,ref_channel,ref_guild,ref_member,should_warn",
    all_cases,
    indirect=['permission']
)
async def test_bans_caps_cmds_flags_ls_mutes_stages_tmutes_commands(
    bot,
    command: Optional[str],
    guild,
    moderation_type,
    not_privileged_author,
    permission,
    prefix: Optional[str],
    privileged_author,
    ref_channel,
    ref_guild,
    ref_member,
    should_warn,
    text_channel,
    voice_channel_one
):
    match moderation_type:
        case "tmute":
            moderation_type = "text_mute"
        case "vmute":
            moderation_type = "voice_mute"
        case "untmute":
            moderation_type = "untext_mute"
        case "unvmute":
            moderation_type = "unvoice_mute"
    channel_values = (voice_channel_one.mention, voice_channel_one.id)
    guild_values = (guild.name, guild.id)
    member_values = (not_privileged_author.mention, not_privileged_author.id)
    formatted = command.format(
        voice_channel_one_id=voice_channel_one.id,
        member_id=not_privileged_author.id,
        guild_id=guild.id
    )
    captured = await prepared_command_handling(author=privileged_author, bot=bot, channel=text_channel, content=formatted, guild=guild, highest_role=permission, prefix=prefix)
    message = captured[0]['message']
    if message.embeds:
        embed = message.embeds[0]
        content = extract_embed_text(embed)
    elif message.embed:
        content = extract_embed_text(message.embed)
    else:
        content = message.content
    message_type = captured[0]['type']
    if message_type == "error":
        print(f"{RED}Error:{RESET} {content}")
    if message_type == "warning":
        # print(f"{YELLOW}Warning:{RESET} {content}")
        if should_warn:
            assert True
    if message_type == "success":
        # print(f"{GREEN}Success:{RESET} {content}")
        if ref_channel:
            assert any(str(channel_value) in content for channel_value in channel_values)
        if ref_guild:
            assert any(str(guild_value) in content for guild_value in guild_values)
        if ref_member:
            assert any(str(member_value) in content for member_value in member_values)
        if moderation_type:
            assert moderation_type in content 