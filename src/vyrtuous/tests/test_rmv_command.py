from discord.ext import commands
from discord.ext.commands import Context, view as cmd_view
from types import SimpleNamespace
from typing import Optional
from unittest.mock import AsyncMock, PropertyMock, patch
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bot.discord_client import DiscordClient
from vyrtuous.cogs.admin_commands import AdminCommands
from vyrtuous.config import Config
from vyrtuous.inc.helpers import *
from vyrtuous.service.discord_message_service import DiscordMessageService
from vyrtuous.tests.test_suite import *
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.database import Database
from vyrtuous.utils.duration import Duration
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.stage import Stage
from vyrtuous.utils.statistics import Statistics
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.utils.setup_logging import logger, setup_logging

import asyncio
import discord
import pytest
import pytest_asyncio
import vyrtuous

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,source_id,target_id,channel_ref_one,channel_ref_two",
    [
        ("rmv {source_id} {target_id}", None, None, "client_channel", "bot_channel")
    ]
)

async def test_rmv_command(bot, bot_channel, client_channel, guild, self_member, prefix: Optional[str], command: Optional[str], source_id, target_id, channel_ref_one, channel_ref_two):
    await admin_initiation(guild.id, self_member.id)
    try:
        source_id = client_channel.id
        target_id = bot_channel.id
        formatted = command.format(
            source_id=source_id,
            target_id=target_id
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
            with patch.object(cog_instance.channel_service, "resolve_channel", return_value=client_channel):
                client_channel.type = discord.ChannelType.voice
                with patch("vyrtuous.cogs.admin_commands.isinstance", side_effect=lambda obj, cls: True if cls == discord.VoiceChannel else isinstance(obj, cls)):
                    await bot.invoke(ctx)
        response = client_channel.messages[0]
        assert any(emoji in response for emoji in Emojis.EMOJIS)
        channel_value_one = client_channel.mention if channel_ref_one else client_channel.name
        channel_value_two = bot_channel.mention if channel_ref_two else bot_channel.name
        assert any(val in response for val in [channel_value_one])
        assert any(val in response for val in [channel_value_two])
        client_channel.messages.clear() 
    finally:
        await admin_cleanup(guild.id, self_member.id)