''' test_logs_command.py The purpose of this program is to black box test the logs command.
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
import pytest_asyncio
import uuid

@pytest_asyncio.fixture(scope="function")
async def active_mlog(bot, client_channel, text_channel, guild, self_member, prefix: str):
    mock_bot_user = SimpleNamespace(id=123456789, bot=True)
    async def capturing_send(self, ctx, content=None, **kwargs):
        client_channel.messages.append(content)
        return MockMessage(
            content=content,
            channel=ctx.channel,
            guild=ctx.guild,
            id=str(uuid.uuid4()),
            author=self_member
        )
    with patch.object(type(bot), "user", new_callable=PropertyMock) as mock_user:
        mock_user.return_value = mock_bot_user
        ctx = await bot.get_context(MockMessage(
            content=f"{prefix}mlog {text_channel.id} create general",
            channel=client_channel,
            guild=guild,
            id=str(uuid.uuid4()),
            author=self_member
        ))
        cog_instance = bot.get_cog("AdminCommands")
        cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
        with patch.object(cog_instance.channel_service, "resolve_channel", return_value=text_channel):
            await bot.invoke(ctx)
        client_channel.messages.clear() 
    yield text_channel

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref",
    [
        ("logs", "all"),
        ("logs", True)
    ]
)

async def test_logs_command(bot, bot_channel, client_channel, text_channel, guild, self_member, dummy_member, prefix: Optional[str], command: Optional[str], channel_ref, active_mlog):
    await admin_initiation(guild.id, self_member.id)
    try:
        if channel_ref == "all":
            arg = "all"
        else:
            arg = active_mlog.id
        formatted = f"{command} {arg}".strip()
        mock_message = MockMessage(content=f"{prefix}{formatted}", channel=client_channel, guild=guild, id='123456789', author=self_member)
        view = cmd_view.StringView(mock_message.content)
        view.skip_string(prefix) 
        async def capturing_send(self, ctx, content=None, embed=None, allowed_mentions=None, **kwargs):
            client_channel.messages.append({'content': content, 'embed': embed})
            return MockMessage(
                content=content,
                channel=ctx.channel,
                embeds=[embed] if embed else [],
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
            fake_channels = {guild.id: [text_channel]}
            with patch("vyrtuous.utils.statistics.Statistics.get_statistic_channels", return_value=fake_channels):
                with patch.object(cog_instance.channel_service, "resolve_channel", return_value=text_channel_obj):
                    await bot.invoke(ctx)
        response = client_channel.messages[0]
        embed = response["embed"]
        assert any(emoji in embed.title for emoji in Emojis.EMOJIS)
        field_text = "".join(f.name + f.value for f in embed.fields)
        assert str(text_channel.id) in field_text
        client_channel.messages.clear() 
    finally:
        await admin_cleanup(guild.id, self_member.id)