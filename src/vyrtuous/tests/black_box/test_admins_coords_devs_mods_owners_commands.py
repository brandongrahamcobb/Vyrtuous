"""test_chown_temp_temps_xtemp_commands.py The purpose of this program is to black box test the temporary room commands.
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
"""

from typing import Optional
from vyrtuous.inc.helpers import *
from vyrtuous.tests.black_box.make_mock_objects import *
from vyrtuous.database.roles.moderator import Moderator
from vyrtuous.tests.black_box.test_suite import *
from vyrtuous.utils.emojis import get_random_emoji, EMOJIS
import pytest


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "permission,command,ref_channel,ref_guild,ref_member,ref_role,should_warn",
    [
        ("Developer", "arole {role_id}", False, False, False, True, False),
        ("Administrator", "aroles {guild_id}", False, True, False, False, False),
        ("Developer", "aroles all", False, False, False, False, False),
        (None, "admins", False, True, False, False, False),
        (None, "admins {member_id}", False, True, True, False, False),
        ("Administrator", "admins {guild_id}", False, True, False, False, False),
        ("Developer", "admins all", False, False, False, False, False),
        ("Developer", "arole {role_id}", False, False, False, True, False),
        ("Administrator", "aroles {guild_id}", False, True, False, False, True),
        ("Developer", "aroles all", False, False, False, False, True),
        (None, "admins", False, True, False, False, True),
        (None, "admins {member_id}", False, True, True, False, True),
        ("Administrator", "admins {guild_id}", False, True, False, False, True),
        ("Developer", "admins all", False, False, False, False, True),
        (
            "Administrator",
            "coord {member_id} {voice_channel_one_id}",
            True,
            False,
            True,
            False,
            False,
        ),
        (None, "coords {voice_channel_one_id}", True, True, False, False, False),
        (None, "coords {member_id}", False, True, True, False, False),
        ("Administrator", "coords {guild_id}", False, True, False, False, False),
        ("Developer", "coords all", False, False, False, False, False),
        (
            "Administrator",
            "coord {member_id} {voice_channel_one_id}",
            True,
            False,
            True,
            False,
            False,
        ),
        (None, "coords {voice_channel_one_id}", True, True, False, False, True),
        (None, "coords {member_id}", False, True, True, False, True),
        ("Administrator", "coords {guild_id}", False, True, False, False, True),
        ("Developer", "coords all", False, False, False, False, True),
        # ("System Owner", "dev {member_id}", False, True, True, False, False)
        (None, "devs", False, True, False, False, True),
        ("Administrator", "devs {guild_id}", False, True, False, False, True),
        (None, "devs {member_id}", False, True, True, False, True),
        ("Developer", "devs all", False, False, False, False, True),
        # ("System Owner", "dev {member_id}", False, True, True, False, False)
        (None, "devs", False, True, False, False, True),
        ("Administrator", "devs {guild_id}", False, True, False, False, True),
        (None, "devs {member_id}", False, True, True, False, True),
        ("Developer", "devs all", False, False, False, False, True),
        (
            "Administrator",
            "mod {member_id} {voice_channel_one_id}",
            True,
            False,
            True,
            False,
            False,
        ),
        (None, "mods {voice_channel_one_id}", True, True, False, False, False),
        (None, "mods {member_id}", False, True, True, False, False),
        ("Administrator", "mods {guild_id}", False, True, False, False, False),
        ("Developer", "mods all", False, False, False, False, False),
        (
            "Administrator",
            "mod {member_id} {voice_channel_one_id}",
            True,
            False,
            True,
            False,
            False,
        ),
        (None, "mods {voice_channel_one_id}", True, True, False, False, True),
        (None, "mods {member_id}", False, True, True, False, True),
    ],
    indirect=["permission"],
)
async def test_admins_coords_devs_mods_owners_commands(
    bot,
    command: Optional[str],
    guild,
    not_privileged_author,
    permission,
    privileged_author,
    prefix: Optional[str],
    ref_channel,
    ref_guild,
    ref_member,
    ref_role,
    should_warn,
    text_channel,
    voice_channel_one,
):
    channel_values = (
        voice_channel_one.mention,
        voice_channel_one.id,
        voice_channel_one.name,
    )
    guild_values = (guild.id, guild.name)
    role_values = (ROLE_ID, ROLE_NAME)
    formatted = command.format(
        voice_channel_one_id=voice_channel_one.id,
        guild_id=guild.id,
        member_id=not_privileged_author.id,
        role_id=ROLE_ID,
    )
    captured = await prepared_command_handling(
        author=privileged_author,
        bot=bot,
        channel=text_channel,
        content=formatted,
        guild=guild,
        highest_role=permission,
        prefix=prefix,
    )
    message = captured[0]["message"]
    if message.embeds:
        embed = message.embeds[0]
        content = extract_embed_text(embed)
    elif message.embed:
        content = extract_embed_text(message.embed)
    else:
        content = message.content
    message_type = captured[0]["type"]
    if message_type == "error":
        print(f"{RED}Error:{RESET} {content}")
    if message_type == "warning":
        # print(f"{YELLOW}Warning:{RESET} {content}")
        if should_warn:
            assert True
        else:
            assert False
    if message_type == "success":
        # print(f"{GREEN}Success:{RESET} {content}")
        if ref_channel:
            assert any(
                str(channel_value) in content for channel_value in channel_values
            )
        if ref_guild:
            assert any(str(guild_value) in content for guild_value in guild_values)
        if ref_role:
            assert any(str(role_value) in content for role_value in role_values)
        assert any(emoji in content for emoji in EMOJIS)
