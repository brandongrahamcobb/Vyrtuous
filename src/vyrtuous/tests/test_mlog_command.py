''' test_mlog_command.py The purpose of this program is to black box test the mlog command.
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

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,action,target_type,target_id,channel_ref,member_ref",
    [
        ("mlog", "create", "general", None, None, None),
        ("mlog", "modify", "general", None, None, None),
        ("mlog", "delete", "general", None, None, None),
        ("mlog", "create", "channel", "client_channel", "client_channel", None),
        ("mlog", "modify", "channel", "client_channel", "client_channel", None),
        ("mlog", "delete", None, "client_channel", "client_channel", None),
        ("mlog", "create", "member", "client_channel", None, None),
        ("mlog", "modify", "member", "client_channel", None, None),
        ("mlog", "delete", None, "client_channel", None, None)
    ]
)

async def test_mlog_command(bot, bot_channel, client_channel, text_channel, guild, self_member, dummy_member, prefix: Optional[str], command: Optional[str], action: Optional[str], target_type: Optional[str], target_id: Optional[str], channel_ref, member_ref):
    await admin_initiation(guild.id, self_member.id)
    try:
        client_channel.messages.clear() 
        target_obj = {
            "client_channel": client_channel,
            "self_member": self_member,
            "dummy_member": dummy_member
        }.get(target_id, None)
        if target_obj:
            target_str = target_obj.id
            formatted = f"{command} {text_channel.id} {action} {target_type} {target_str}".strip()
        else:
            formatted = f"{command} {text_channel.id} {action} {target_type}".strip()
        mock_message = MockMessage(content=f"{prefix}{formatted}", channel=client_channel, guild=guild, id='123456789', author=self_member)
        view = cmd_view.StringView(mock_message.content)
        view.skip_string(prefix) 
        async def capturing_send(self, ctx, content=None, allowed_mentions=None, **kwargs):
            client_channel.messages.append(content)
            return MockMessage(
                content=content,
                channel=ctx.channel,
                guild=ctx.guild,
                id='123456789',
                author=self_member
            )
        mock_bot_user = SimpleNamespace(id='123456789', bot=True)
        with patch.object(type(bot), "user", new_callable=PropertyMock) as mock_user:
            mock_user.return_value = mock_bot_user
            ctx = await bot.get_context(mock_message)
            cog_instance = bot.get_cog("AdminCommands")
            cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
            with patch.object(cog_instance.channel_service, "resolve_channel", return_value=text_channel):
                text_channel.type = discord.ChannelType.text
                await bot.invoke(ctx)
        response = client_channel.messages[0]
        print(response)
        channel_value = text_channel.mention if channel_ref else text_channel.name
        self_member_value = self_member.mention if member_ref else self_member.name
        dummy_member_value = dummy_member.mention if member_ref else dummy_member.name
        assert any(emoji in response for emoji in Emojis.EMOJIS)
        if channel_ref:
            assert any(val in response for val in [channel_value])
        if member_ref == "dummy_member":
            assert any(val in response for val in [dummy_member_value])
        if member_ref == "self_member":
            assert any(val in response for val in [self_member_value])
        client_channel.messages.clear() 
    finally:
        await admin_cleanup(guild.id, self_member.id)