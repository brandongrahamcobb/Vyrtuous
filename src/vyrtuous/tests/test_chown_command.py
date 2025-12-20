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
from vyrtuous.tests.test_suite import MockMessage, bot, bot_channel, client_channel, client, config, self_member, dummy_member, guild
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
        ("temp {client_channel_id} {member_id}", None, "member"),
        ("chown {client_channel_id} {member_id}", "channel", "member"),
        ("xtemp {client_channel_id} {client_channel_id}", "channel", None),
        ("temp {channel_mention} {member_mention}", None, "member"),
        ("chown {channel_mention} {member_mention}", "channel", "member")
    ]
)

async def test_chown_commands(bot, bot_channel, client_channel, guild, self_member, prefix: Optional[str], command: Optional[str], channel_ref, member_ref):
    formatted = command.format(
        bot=bot,
        bot_channel_id=bot_channel.id,
        client_channel_id=client_channel.id,
        channel_mention=client_channel.mention,
        member_id=self_member.id,
        member_mention=self_member.mention,
    )
    await admin_initiation(guild.id, self_member.id)
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
    channel_value = client_channel.mention
    member_value = self_member.mention
    assert any(emoji in response for emoji in Emojis.EMOJIS) 
    assert any(val in response for val in [channel_value])
    assert any(val in response for val in [member_value])
    await admin_cleanup(guild.id, self_member.id)