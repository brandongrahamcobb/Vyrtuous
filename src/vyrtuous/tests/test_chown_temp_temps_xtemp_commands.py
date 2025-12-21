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

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,channel_ref,member_ref",
    [
        ("temp {client_channel_id} {member_id}", "channel", "member"),
        ("temps all", None, None),
        ("temps {client_channel_id}", None, None),
        ("chown {client_channel_id} {member_id}", "channel", "member"),
        ("temp {client_channel_id}", "channel", None),
        ("temp {channel_mention} {member_mention}", "channel", "member"),
        ("chown {channel_mention} {member_mention}", "channel", "member"),
        ("temp {client_channel_id}", "channel", None)
    ]
)

async def test_chown_temp_xtemp_commands(bot, bot_channel, client_channel, guild, self_member, prefix: Optional[str], command: Optional[str], channel_ref, member_ref):    
    await admin_initiation(guild.id, self_member.id)
    try:
        client_channel.messages.clear() 
        formatted = command.format(
            bot=bot,
            bot_channel_id=bot_channel.id,
            client_channel_id=client_channel.id,
            channel_mention=client_channel.mention,
            member_id=self_member.id,
            member_mention=self_member.mention,
        )
        mock_message = MockMessage(content=f"{prefix}{formatted}", channel=client_channel, guild=guild, id='123456789', author=self_member)
        view = cmd_view.StringView(mock_message.content)
        view.skip_string(prefix) 
        async def mock_wait_for(event, timeout=None, check=None):
            raise asyncio.TimeoutError()
        bot.wait_for = mock_wait_for
        async def mock_channel_send(content=None, embed=None, embeds=None, **kwargs):
            print(f"Channel.send called - content: {content}, embed: {embed}")
            if content:
                client_channel.messages.append(content)
            elif embed:
                client_channel.messages.append(f"[Paginator embed]")
            return MockMessage(
                content=content or "",
                channel=client_channel,
                guild=guild,
                id='123456789',
                author=self_member,
                embeds=[embed] if embed else []
            )
        original_channel_send = client_channel.send
        client_channel.send = mock_channel_send
        async def capturing_send(self, ctx, content=None, embed=None, allowed_mentions=None, **kwargs):
            client_channel.messages.append({'content': content, 'embed': embed})
            return MockMessage(
                content=content,
                channel=ctx.channel,
                guild=ctx.guild,
                id='123456789',
                author=self_member,
                allowed_mentions=allowed_mentions,
                embeds=[embed]
            )
        mock_bot_user = SimpleNamespace(id='123456789', bot=True)
        with patch.object(type(bot), "user", new_callable=PropertyMock) as mock_user:
            mock_user.return_value = mock_bot_user
            ctx = await bot.get_context(mock_message)
            ctx.send = mock_channel_send
            cog_instance = bot.get_cog("AdminCommands")
            cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
            with patch.object(cog_instance.channel_service, "resolve_channel", return_value=client_channel):
                client_channel.type = discord.ChannelType.voice
                with patch("vyrtuous.cogs.admin_commands.isinstance", side_effect=lambda obj, cls: True if cls == discord.VoiceChannel else isinstance(obj, cls)):
                    await bot.invoke(ctx)
        response = client_channel.messages[0]
        print(response)
        channel_value = client_channel.mention if channel_ref else client_channel.name
        member_value = self_member.mention if member_ref else self_member.name
        if "temps" in command:
            assert "Paginator embed" in response
        else:
            assert any(emoji in response['content'] for emoji in Emojis.EMOJIS) 
            assert any(val in response['content'] for val in [channel_value])
        if member_ref:
            assert any(val in response['content'] for val in [member_value])
        client_channel.messages.clear() 
    finally:
        await admin_cleanup(guild.id, self_member.id)