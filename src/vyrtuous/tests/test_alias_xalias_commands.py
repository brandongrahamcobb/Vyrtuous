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
from unittest.mock import PropertyMock, patch
from vyrtuous.inc.helpers import *
from vyrtuous.tests.test_admin_helpers import *
from vyrtuous.tests.test_suite import *
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
    voice_channel_one.messages.clear() 
    try:
        channel_token = voice_channel_one.mention if channel_ref else ""
        if role_ref:
            formatted = f"{command} {alias_type} {alias_name} {channel_token} {ROLE_ID}".strip()
        elif command == 'alias':
            formatted = f"{command} {alias_type} {alias_name} {channel_token}".strip()
        elif command == 'xalias':
            formatted = f"{command} {alias_name}".strip()
        mock_message = make_mock_message(allowed_mentions=True, author=privileged_author, channel=voice_channel_one, content=f"{prefix}{formatted}", embeds=[], guild=guild, id=MESSAGE_ID)
        view = cmd_view.StringView(mock_message.content)
        view.skip_string(prefix) 
        mock_bot_user = make_mock_member(id=PRIVILEGED_AUTHOR_ID, name=PRIVILEGED_AUTHOR_NAME)
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
            cog_instance = bot.get_cog("AdminCommands")
            capturing_send = make_capturing_send(voice_channel_one, privileged_author)
            cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
            with patch.object(cog_instance.channel_service, "resolve_channel", return_value=voice_channel_one):
                voice_channel_one.type = discord.ChannelType.voice
                await bot.invoke(ctx)
        response = voice_channel_one.messages[0]
        channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        assert any(emoji in response for emoji in Emojis.EMOJIS)
        if command == "alias":
            if alias_type in ('role', 'unrole'):
                assert str(ROLE_ID) in response or f"<@&{ROLE_ID}>" in response
            else:
                assert any(val in response for val in [channel_value])
    finally:
        await admin_cleanup(guild.id, PRIVILEGED_AUTHOR_ID)