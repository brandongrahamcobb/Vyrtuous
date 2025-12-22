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
from discord.ext.commands import Context, view as cmd_view
from types import SimpleNamespace
from typing import Optional
from vyrtuous.inc.helpers import *
from vyrtuous.tests.test_admin_helpers import admin_cleanup, admin_initiation
from vyrtuous.tests.test_suite import bot, config, guild, prepared_command_handling, prefix, privileged_author, voice_channel_one
from vyrtuous.utils.emojis import Emojis
import discord
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,alias_type,alias_name,channel_ref,role_ref",
    [
        ("alias", "ban", "testban", True, False),
        ("xalias", None, "testban", None, None),
        ("alias", "unban", "testunban", True, False),
        ("xalias", None, "testunban", None, None),
        ("alias", "mute", "testmute", True, False),
        ("xalias", None, "testmute", None, None),
        ("alias", "unmute", "testunmute", True, False),
        ("xalias", None, "testunmute", None, None),
        ("alias", "flag", "testflag", True, False),
        ("xalias", None, "testflag",  None, None),
        ("alias", "unflag", "testunflag", True, False),
        ("xalias", None, "testunflag", None, None),
        ("alias", "cow", "testcow", True, False),
        ("xalias", None, "testcow", None, None),
        ("alias", "uncow", "testuncow", True, False),
        ("xalias", None, "testuncow", None, None),
        ("alias", "tmute", "testtmute", True, False),
        ("xalias", None, "testtmute", None, None),
        ("alias", "untmute", "testuntmute", True, False),
        ("xalias", None, "testuntmute", None, None),
        ("alias", "role", "testrole", True, True),
        ("xalias", None, "testrole", None, None),
        ("alias", "unrole", "testunrole", True, True),
        ("xalias", None, "testunrole", None, None)
    ]
)

async def test_alias_xalias_command(bot, voice_channel_one, guild, privileged_author, prefix: Optional[str], command: Optional[str], alias_type, alias_name, channel_ref, role_ref):
    await admin_initiation(guild.id, privileged_author.id)
    try:
        channel_token = voice_channel_one.mention
        formatted = "{command} {alias_type} {alias_name} {channel} {role}".format(
            alias_name=alias_name,
            alias_type=alias_type or "",
            channel=voice_channel_one.mention if channel_ref else "",
            command=command,
            role=str(ROLE_ID) if role_ref else ""
        ).strip()
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]["content"]
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        assert any(emoji in response for emoji in Emojis.EMOJIS)
        if command == "alias":
            if alias_type in ('role', 'unrole'):
                assert str(ROLE_ID) in response or f"<@&{ROLE_ID}>" in response
            else:
                assert any(val in response for val in [channel_value])
        voice_channel_one.messages.clear() 
    finally:
        await admin_cleanup(guild.id, PRIVILEGED_AUTHOR_ID)