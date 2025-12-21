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
    yield text_channel


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref",
    [
        ("log", True),
        ("log", True)
    ]
)

async def test_log_command(bot, bot_channel, client_channel, text_channel, guild, self_member, dummy_member, prefix: Optional[str], command: Optional[str], channel_ref, active_mlog):
    await admin_initiation(guild.id, self_member.id)
    try:
        channel = active_mlog
        formatted = f"{command} {channel.id}".strip()
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
            fake_channels = {guild.id: [text_channel]}
            with patch("vyrtuous.utils.statistics.Statistics.get_statistic_channels", return_value=fake_channels):
                with patch.object(cog_instance.channel_service, "resolve_channel", return_value=text_channel_obj):
                    await bot.invoke(ctx)
        response = client_channel.messages[0]
        print(response)
        channel_value = text_channel.mention if channel_ref else text_channel.name
        assert any(emoji in response for emoji in Emojis.EMOJIS)
        assert any(val in response for val in [channel_value])
        client_channel.messages.clear() 
    finally:
        await admin_cleanup(guild.id, self_member.id)