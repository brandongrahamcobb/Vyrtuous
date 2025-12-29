''' test_coord_command.py The purpose of this program is to black box test the coord command.
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
from vyrtuous.tests.black_box.test_suite import *
from vyrtuous.utils.administrator import Administrator
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref,member_ref",
    [
        ("coord {not_privileged_author_id} {voice_channel_one_id}", True, True),
        ("coord {not_privileged_author_id} {voice_channel_one_id}", True, True)
    ]
)

async def test_coord_command(bot, voice_channel_one, guild, privileged_author, not_privileged_author, prefix: Optional[str], role, command: Optional[str], channel_ref, member_ref):
    administrator = Administrator(guild_snowflake=guild.id, member_snowflake=privileged_author.id, role_snowflake=role.id)
    await administrator.grant()
    try:
        voice_channel_one.messages.clear()
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id,
            not_privileged_author_id=not_privileged_author.id
        )
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]["content"]
        assert any(emoji in response for emoji in Emojis.EMOJIS)
    finally:
        await administrator.revoke()