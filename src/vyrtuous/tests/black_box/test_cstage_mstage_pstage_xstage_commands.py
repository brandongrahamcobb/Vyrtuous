''' test_cstage_xstage_commands.py The purpose of this program is to black box test the stages commands.
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
from vyrtuous.tests.black_box.test_admin_helpers import admin_cleanup, admin_initiation
from vyrtuous.utils.moderator import Moderator
from vyrtuous.tests.black_box.test_suite import bot, config, guild, not_privileged_author, prepared_command_handling, prefix, privileged_author, voice_channel_one
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,duration,channel_ref,member_ref",
    [
        ("cstage {voice_channel_one_id}", '1m', True, False),
        ("mstage {not_privileged_author_id}", None, False, True),
        ("pstage {not_privileged_author_id}", None, False, True),
        ("xstage {voice_channel_one_id}", None, True, False),
        ("cstage {voice_channel_one_id}", '1h', True, False),
        ("xstage {voice_channel_one_id}", None, True, False),
        ("cstage {voice_channel_one_id}", '1d', True, False),
        ("xstage {voice_channel_one_id}", None, True, False)
    ]
)

async def test_cstage_mstage_pstage_xstage_command(bot, voice_channel_one, guild, not_privileged_author, privileged_author, prefix: Optional[str], command: Optional[str], duration, channel_ref, member_ref):
    try:
        voice_channel_one.messages.clear() 
        formatted = command.format(
                not_privileged_author_id=not_privileged_author.id,
                voice_channel_one_id=voice_channel_one.id,
                duration=duration
            )
        if "cstage" in command or "xstage" in command:
            await admin_initiation(guild.id, privileged_author.id)
            await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
        else:
            await Moderator.grant(channel_snowflake=voice_channel_one.id, guild_snowflake=guild.id, member_snowflake=privileged_author.id)
            await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="ModeratorCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.moderator_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]["content"]
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        member_value = not_privileged_author.mention if member_ref else not_privileged_author.name
        assert any(emoji in response for emoji in Emojis.EMOJIS)
        if channel_ref:
            assert any(val in response for val in [channel_value])
        elif member_ref:
            assert any(val in response for val in [member_value])
    finally:
        if "cstage" in command or "xstage" in command:
            await admin_cleanup(guild.id, privileged_author.id)
        else:
            await Moderator.revoke(channel_snowflake=voice_channel_one.id, member_snowflake=privileged_author.id)