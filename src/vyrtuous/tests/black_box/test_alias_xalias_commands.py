''' test_alias_xalias_commands.py The purpose of this program is to black box test the Aliases module.
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
from vyrtuous.tests.black_box.test_suite import bot, config, guild, prepared_command_handling, prefix, privileged_author, role, voice_channel_one
from vyrtuous.utils.administrator import Administrator
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,alias_type,alias_name,channel_ref,role_ref",
    [
        ("alias", "ban", "testban", True, False),
        ("xalias", None, "testban", False, False),
        ("alias", "unban", "testunban", True, False),
        ("xalias", None, "testunban", False, False),
        ("alias", "mute", "testmute", True, False),
        ("xalias", None, "testmute", False, False),
        ("alias", "unmute", "testunmute", True, False),
        ("xalias", None, "testunmute", False, False),
        ("alias", "flag", "testflag", True, False),
        ("xalias", None, "testflag",  False, False),
        ("alias", "unflag", "testunflag", True, False),
        ("xalias", None, "testunflag", False, False),
        ("alias", "cow", "testcow", True, False),
        ("xalias", None, "testcow", False, False),
        ("alias", "uncow", "testuncow", True, False),
        ("xalias", None, "testuncow", False, False),
        ("alias", "tmute", "testtmute", True, False),
        ("xalias", None, "testtmute", False, False),
        ("alias", "untmute", "testuntmute", True, False),
        ("xalias", None, "testuntmute", False, False),
        ("alias", "role", "testrole", True, True),
        ("xalias", None, "testrole", False, False),
        ("alias", "unrole", "testunrole", True, True),
        ("xalias", None, "testunrole", False, False)
    ]
)

async def test_alias_xalias_command(bot, voice_channel_one, guild, privileged_author, prefix: Optional[str], command: Optional[str], role, alias_type, alias_name, channel_ref, role_ref):
    administrator = Administrator(guild_snowflake=guild.id, member_snowflake=privileged_author.id, role_snowflake=role.id)
    await administrator.grant()
    try:
        voice_channel_one.messages.clear()
        if channel_ref and role_ref:
            formatted = f"{command} {alias_name} {alias_type} {voice_channel_one.id} {ROLE_ID}"
        elif channel_ref:
            formatted = f"{command} {alias_name} {alias_type} {voice_channel_one.id}"
        else:
            formatted = f"{command} {alias_name}"
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]["content"]
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        assert any(emoji in response for emoji in Emojis.EMOJIS)
    finally:
        await administrator.revoke()