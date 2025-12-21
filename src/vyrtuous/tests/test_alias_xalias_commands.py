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
from discord.ext.commands import view as cmd_view
from types import SimpleNamespace
from typing import Optional
from unittest.mock import PropertyMock, patch
from vyrtuous.tests.test_suite import *
from vyrtuous.utils.emojis import Emojis
import discord
import pytest
import uuid

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

async def test_alias_xalias_command(bot, bot_channel, client_channel, guild, self_member, prefix: Optional[str], command: Optional[str], alias_type, alias_name, channel_ref, role_ref):
    await admin_initiation(guild.id, self_member.id)
    client_channel.messages.clear() 
    try:
        channel_token = client_channel.mention if channel_ref else ""
        role_token = "987654321" if role_ref else ""
        if role_token:
            formatted = f"{command} {alias_type} {alias_name} {channel_token} {role_token}".strip()
        elif command == 'alias':
            formatted = f"{command} {alias_type} {alias_name} {channel_token}".strip()
        elif command == 'xalias':
            formatted = f"{command} {alias_name}".strip()
        print(formatted)
        id = str(uuid.uuid4())
        mock_message = MockMessage(content=f"{prefix}{formatted}", channel=client_channel, guild=guild, id=id, author=self_member)
        view = cmd_view.StringView(mock_message.content)
        view.skip_string(prefix) 
        async def capturing_send(self, ctx, content=None, allowed_mentions=None, **kwargs):
            client_channel.messages.append(content)
            return MockMessage(
                content=content,
                channel=ctx.channel,
                guild=ctx.guild,
                id=id,
                author=self_member
            )
        mock_bot_user = SimpleNamespace(id='123456789', bot=True)
        with patch.object(type(bot), "user", new_callable=PropertyMock) as mock_user:
            mock_user.return_value = mock_bot_user
            ctx = await bot.get_context(mock_message)
            cog_instance = bot.get_cog("AdminCommands")
            cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
            with patch.object(cog_instance.channel_service, "resolve_channel", return_value=client_channel):
                client_channel.type = discord.ChannelType.voice
                await bot.invoke(ctx)
        response = client_channel.messages[0]
        channel_value = client_channel.mention if channel_ref else client_channel.name
        assert any(emoji in response for emoji in Emojis.EMOJIS)
        if command == "alias":
            if alias_type in ('role', 'unrole'):
                assert role_token in response or f"<@&{role_token}>" in response
            else:
                assert any(val in response for val in [channel_value])
    finally:
        await admin_cleanup(guild.id, self_member.id)