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
from discord.ext.commands import Context, view as cmd_view
from types import SimpleNamespace
from typing import Optional
from unittest.mock import PropertyMock, patch
from vyrtuous.inc.helpers import *
from vyrtuous.tests.make_mock_objects import *
from vyrtuous.tests.test_admin_helpers import *
from vyrtuous.tests.test_suite import *
from vyrtuous.utils.emojis import Emojis
import asyncio
import discord
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref,member_ref",
    [
        ("temp {voice_channel_one_id} {member_id}", "channel", "member"),
        ("temps all", None, None),
        ("temps {voice_channel_one_id}", None, None),
        ("chown {voice_channel_one_id} {member_id}", "channel", "member"),
        ("temp {voice_channel_one_id}", "channel", None),
        ("temp {channel_mention} {member_mention}", "channel", "member"),
        ("chown {channel_mention} {member_mention}", "channel", "member"),
        ("temp {voice_channel_one_id}", "channel", None)
    ]
)

async def test_chown_temp_xtemp_commands(bot, voice_channel_one, guild, privileged_author, prefix: Optional[str], command: Optional[str], channel_ref, member_ref):    
    await admin_initiation(guild.id, privileged_author.id)
    try:
        voice_channel_one.messages.clear() 
        formatted = command.format(
            voice_channel_one_id=voice_channel_one.id,
            channel_mention=voice_channel_one.mention,
            member_id=privileged_author.id,
            member_mention=privileged_author.mention,
        )
        mock_message = make_mock_message(allowed_mentions=True, author=privileged_author, channel=voice_channel_one, content=f"{prefix}{formatted}", embeds=[], guild=guild, id=MESSAGE_ID)
        view = cmd_view.StringView(mock_message.content)
        view.skip_string(prefix) 
        bot.wait_for = mock_wait_for
        async def mock_channel_send(content=None, embed=None, embeds=None, **kwargs):
            print(f"Channel.send called - content: {content}, embed: {embed}")
            if content:
                voice_channel_one.messages.append(content)
            elif embed:
                voice_channel_one.messages.append(f"[Paginator embed]")
            return make_mock_message(allowed_mentions=True, author=privileged_author, channel=voice_channel_one, content=content, embeds=[embed], guild=guild, id=MESSAGE_ID)
        original_channel_send = voice_channel_one.send
        voice_channel_one.send = mock_channel_send
        mock_bot_user = make_mock_member(id=PRIVILEGED_AUTHOR_ID, name=PRIVILEGED_AUTHOR_NAME)
        capturing_send = make_capturing_send(voice_channel_one, privileged_author)
        with patch.object(type(bot), "user", new_callable=PropertyMock) as mock_user:
            mock_user.return_value = mock_bot_user
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
            ctx.send = mock_channel_send
            cog_instance = bot.get_cog("AdminCommands")
            cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
            with patch.object(cog_instance.channel_service, "resolve_channel", return_value=voice_channel_one):
                voice_channel_one.type = discord.ChannelType.voice
                with patch("vyrtuous.cogs.admin_commands.isinstance", side_effect=lambda obj, cls: True if cls == discord.VoiceChannel else isinstance(obj, cls)):
                    await bot.invoke(ctx)
        response = voice_channel_one.messages[0]
        print(response)
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        member_value = privileged_author.mention if member_ref else privileged_author.name
        if "temps" in command:
            assert "Paginator embed" in response
        else:
            assert any(emoji in response for emoji in Emojis.EMOJIS) 
            assert any(val in response for val in [channel_value])
        if member_ref:
            assert any(val in response for val in [member_value])
    finally:
        await admin_cleanup(guild.id, privileged_author.id)