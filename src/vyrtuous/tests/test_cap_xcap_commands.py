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
from discord.ext.commands import view as cmd_view
from types import SimpleNamespace
from typing import Optional
from unittest.mock import PropertyMock, patch
from vyrtuous.tests.test_suite import *
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref,member_ref",
    [
        ("cap {client_channel_id} {arg}", "channel", "member"),
        ("xcap {client_channel_id} {arg}", "channel", None),
        ("cap {client_channel_id} {arg}", None, "member"),
        ("xcap {client_channel_id} {arg}", None, None)
    ]
)

async def test_cap_commands(bot, bot_channel, client_channel, guild, self_member, prefix: Optional[str], command: Optional[str], channel_ref, member_ref):
    await admin_initiation(guild.id, self_member.id)
    try:
        durations = ['1m', '1h', '1d']
        moderation_types = ['ban', 'mute', 'tmute']
        args = [f'{m} {d}' for m in moderation_types for d in durations]
        for arg in args:
            formatted = command.format(
                bot=bot,
                bot_channel_id=bot_channel.id,
                client_channel_id=client_channel.id,
                channel_mention=client_channel.mention,
                member_id=self_member.id,
                member_mention=self_member.mention,
                arg=arg
            )
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
                await bot.invoke(ctx)
            response = client_channel.messages[0]
            channel_value = client_channel.mention if channel_ref else client_channel.name
            member_value = self_member.mention if member_ref else self_member.name
            assert any(emoji in response for emoji in Emojis.EMOJIS)
            assert any(duration in response for duration in durations)
            assert any(moderation_type in response for moderation_type in moderation_types)
            assert any(val in response for val in [channel_value])
            assert any(val in response for val in [member_value])
            client_channel.messages.clear() 
    finally:
        await admin_cleanup(guild.id, self_member.id)