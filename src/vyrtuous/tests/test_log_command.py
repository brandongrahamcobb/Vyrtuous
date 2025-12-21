''' test_log_command.py The purpose of this program is to black box test the log command.
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
from vyrtuous.utils.setup_logging import logger, setup_logging

import asyncio
import discord
import pytest
import pytest_asyncio
import uuid

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref",
    [
        ("log", True),
        ("log", True),
        ("logs all", False),
        ("logs", True)
    ]
)

async def test_log_command(bot, voice_channel_one, text_channel, guild, privileged_author, not_privileged_author, prefix: Optional[str], command: Optional[str], channel_ref):
    await admin_initiation(guild.id, privileged_author.id)
    try:
        if not channel_ref:
            formatted = f"{command}".strip()
        else:
            formatted = f"{command} {text_channel.id}".strip()
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
            fake_channels = {
                guild.id: [
                    {
                        "channel_id": text_channel.id,
                        "channel_name": text_channel.name,
                        "enabled": False
                    }
                ]
            }
            with patch("vyrtuous.utils.statistics.Statistics.get_statistic_channels", return_value=fake_channels):
                with patch.object(cog_instance.channel_service, "resolve_channel", return_value=text_channel_obj):
                    await bot.invoke(ctx)
        response = voice_channel_one.messages[0]
        if not channel_ref:
            embed = response["embed"]
            assert any(emoji in embed for emoji in Emojis.EMOJIS)
        else:
            channel_value = text_channel.mention if channel_ref else text_channel.name
            assert any(emoji in response for emoji in Emojis.EMOJIS)
            assert any(val in response for val in [channel_value])
        voice_channel_one.messages.clear() 
    finally:
        await admin_cleanup(guild.id, privileged_author.id)