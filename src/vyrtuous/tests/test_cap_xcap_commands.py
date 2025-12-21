''' test_cap_xcap_commands.py The purpose of this program is to black box test the Cap commands.
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
from unittest.mock import PropertyMock, patch
from vyrtuous.inc.helpers import *
from vyrtuous.tests.test_admin_helpers import *
from vyrtuous.tests.test_suite import *
from vyrtuous.utils.emojis import Emojis
import pytest

def generate_cap_test_cases():
    durations = ['1m', '1h', '1d']
    moderation_types = ['ban', 'mute', 'tmute']
    cases = []
    for mod_type in moderation_types:
        for duration in durations:
            cases.append((f"cap {{voice_channel_one_id}} {mod_type} {duration}", "channel", f"{mod_type} {duration}"))
            cases.append((f"xcap {{voice_channel_one_id}} {mod_type}", "channel", mod_type))
    return cases

@pytest.mark.asyncio
@pytest.mark.parametrize("command,channel_ref,arg", generate_cap_test_cases())
async def test_cap_commands(bot, voice_channel_one, guild, privileged_author, prefix: Optional[str], command: Optional[str], channel_ref, arg):
    await admin_initiation(guild.id, privileged_author.id)
    try:
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id
        )
        mock_message = make_mock_message(allowed_mentions=True, author=privileged_author, channel=voice_channel_one, content=f"{prefix}{formatted}", embeds=[], guild=guild, id=MESSAGE_ID)
        view = cmd_view.StringView(mock_message.content)
        view.skip_string(prefix)
        ctx = Context(
            message=mock_message,
            bot=bot,
            prefix=prefix,
            view=view
        )
        command_name = formatted.split()[0]
        ctx.command = bot.get_command(command_name)
        ctx.invoked_with = command_name
        view.skip_string(command_name)
        view.skip_ws() 
        mock_bot_user = make_mock_member(id=PRIVILEGED_AUTHOR_ID, name=PRIVILEGED_AUTHOR_NAME)
        capturing_send = make_capturing_send(voice_channel_one, privileged_author)
        with patch.object(type(bot), "user", new_callable=PropertyMock) as mock_user:
            mock_user.return_value = mock_bot_user
            cog_instance = bot.get_cog("AdminCommands")
            cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
            await bot.invoke(ctx)
        response = voice_channel_one.messages[0]
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        assert any(emoji in response for emoji in Emojis.EMOJIS)
        if channel_ref:
            assert any(val in response for val in [channel_value])
        voice_channel_one.messages.clear() 
    finally:
        await admin_cleanup(guild.id, privileged_author.id)