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
from discord.ext.commands import Context, view as cmd_view
from types import SimpleNamespace
from typing import Optional
from unittest.mock import PropertyMock, patch
from vyrtuous.tests.make_mock_objects import *
from vyrtuous.tests.test_admin_helpers import *
from vyrtuous.tests.test_suite import *
from vyrtuous.utils.emojis import Emojis
import pytest
import pytest_asyncio
import uuid

# @pytest_asyncio.fixture(scope="function")
# async def active_mlog(bot, voice_channel_one, text_channel, guild, privileged_author, prefix: str):
#     mock_bot_user = SimpleNamespace(id=123456789, bot=True)
#     async def capturing_send(self, ctx, content=None, **kwargs):
#         voice_channel_one.messages.append(content)
#         return MockMessage(
#             content=content,
#             channel=ctx.channel,
#             guild=ctx.guild,
#             id=str(uuid.uuid4()),
#             author=privileged_author
#         )
#     with patch.object(type(bot), "user", new_callable=PropertyMock) as mock_user:
#         mock_user.return_value = mock_bot_user
#         ctx = await bot.get_context(MockMessage(
#             content=f"{prefix}mlog {text_channel.id} create general",
#             channel=voice_channel_one,
#             guild=guild,
#             id=str(uuid.uuid4()),
#             author=privileged_author
#         ))
#         cog_instance = bot.get_cog("AdminCommands")
#         cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
#         with patch.object(cog_instance.channel_service, "resolve_channel", return_value=text_channel):
#             await bot.invoke(ctx)
#         voice_channel_one.messages.clear() 
#     yield text_channel

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref",
    [
        ("logs all", False),
        ("logs", True)
    ]
)

async def test_logs_command(bot, voice_channel_one, text_channel, guild, privileged_author,  prefix: Optional[str], command: Optional[str], channel_ref):
    await admin_initiation(guild.id, privileged_author.id)
    try:
        formatted = f"{command} {arg}".strip()
        mock_message = make_mock_message(allowed_mentions=True, author=privileged_author, channel=voice_channel_one, content=f"{prefix}{formatted}", embeds=[], guild=guild, id=MESSAGE_ID)
        view = cmd_view.StringView(mock_message.content)
        view.skip_string(prefix) 
        capturing_send = make_capturing_send(voice_channel_one, privileged_author)
        mock_bot_user = make_mock_member(id=PRIVILEGED_AUTHOR_ID, name=PRIVILEGED_AUTHOR_NAME)
        with patch.object(type(bot), "user", new_callable=PropertyMock) as mock_user:
            mock_user.return_value = mock_bot_user
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
            cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
            fake_channels = {guild.id: [text_channel]}
            with patch("vyrtuous.utils.statistics.Statistics.get_statistic_channels", return_value=fake_channels):
                with patch.object(cog_instance.channel_service, "resolve_channel", return_value=text_channel_obj):
                    await bot.invoke(ctx)
        response = voice_channel_one.messages[0]
        embed = response["embed"]
        assert any(emoji in embed.title for emoji in Emojis.EMOJIS)
        field_text = "".join(f.name + f.value for f in embed.fields)
        assert str(text_channel.id) in field_text
        voice_channel_one.messages.clear() 
    finally:
        await admin_cleanup(guild.id, privileged_author.id)