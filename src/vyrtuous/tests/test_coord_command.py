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
        ("coord", "channel", "member"),
        ("coord", "channel", None),
        ("coord", None, "member"),
        ("coord", None, None),
    ]
)

async def test_coord_command(bot, bot_channel, client_channel, guild, self_member, prefix: Optional[str], command: Optional[str], channel_ref, member_ref):
    await admin_initiation(guild.id, self_member.id)
    formatted = command.format(
        bot=bot,
        bot_channel_id=bot_channel.id,
        client_channel_id=client_channel.id,
        channel_mention=client_channel.mention,
        member_id=self_member.id,
        member_mention=self_member.mention
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
    channel_value = client_channel.mention
    member_value = self_member.mention
    assert any(emoji in response for emoji in Emojis.EMOJIS)
    assert any(val in response for val in [channel_value])
    assert any(val in response for val in [member_value])
    await admin_cleanup(guild.id, self_member.id)
    #     bot = DiscordBot.get_instance()
    #     try:
    #         await channel.send(f'{prefix}coord {member.id} {channel.id}')
    #         await asyncio.sleep(10)
    #         report.pass_('coord.send')
    #     except Exception as e:
    #         report.error('coord.send', e)
    #         return
    #     try:
    #         async with bot.db_pool.acquire() as conn:
    #             row = await conn.fetchrow('SELECT coordinator_channel_ids FROM users WHERE discord_snowflake=$1', member.id)
    #         if row is not None:
    #             report.pass_('coord.insert')
    #         else:
    #             report.fail('coord.insert')
    #     except Exception as e:
    #         report.error('coord.insert', e)
    #         return
    #     try:
    #         if row['coordinator_channel_ids'] and channel.id in row['coordinator_channel_ids']:
    #             report.pass_('coord.validate')
    #         else:
    #             report.fail('coord.validate')
    #     except Exception as e:
    #         report.error('coord.validate', e)
            
    # # async def cstage_command_test(self, channel: discord.abc.GuildChannel, prefix: Optional[str]):
    #     bot = DiscordBot.get_instance()
    #     durations = ['1m', '1h', '1d'] #'1s', '1y'
    #     for duration in durations:
    #         try:
    #             await channel.send(f'{prefix}cstage {channel.id} {duration}')
    #             await asyncio.sleep(10)
    #         except Exception as e:
    #             report.error('cstage.send', e)
    #             continue
    #         try:
    #             stages = await Stage.fetch_stage_by_guild_id_channel_id_and_channel_name(guild_id=channel.guild.id, channel_id=channel.id, channel_name=channel.name)
    #             if not stages:
    #                 report.fail('cstage.create')
    #                 continue
    #             report.pass_('cstage.create')
    #         except Exception as e:
    #             report.error('cstage.fetch', e)
    #             continue
    #         await channel.send(f'{prefix}cstage {channel.id} ')
    # # async def log_command_test(self, channel: discord.abc.GuildChannel, member: discord.Member, prefix: Optional[str]):
    # # async def logs_command_test(self, channel: discord.abc.GuildChannel, prefix: Optional[str]):
    # async def mlog_command_test(self, first_channel: discord.abc.GuildChannel, second_channel: discord.abc.GuildChannel, text_channel: discord.abc.GuildChannel, self_member: discord.Member, dummy_member: discord.Member, prefix: Optional[str], report: TestReport):
    #     bot = DiscordBot.get_instance()
    #     try:
    #         await first_channel.send(f'{prefix}mlog {text_channel.id} create general')
    #         await asyncio.sleep(10)
    #         async with bot.db_pool.acquire() as conn:
    #             rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2', first_channel.guild.id, text_channel.id)
    #         if rows:
    #             report.pass_('mlog.create.general')
    #         else:
    #             report.fail('mlog.create.general', 'No rows found after create')
    #     except Exception as e:
    #         report.error('mlog.create.general', e)
    #     try:
    #         await first_channel.send(f'{prefix}mlog {text_channel.id} modify general')
    #         await asyncio.sleep(10)
    #         async with bot.db_pool.acquire() as conn:
    #             rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2', first_channel.guild.id, text_channel.id)
    #         if rows:
    #             report.pass_('mlog.modify.general')
    #         else:
    #             report.fail('mlog.modify.general', 'No rows found after modify')
    #     except Exception as e:
    #         report.error('mlog.modify.general', e)
    #     try:
    #         await first_channel.send(f'{prefix}mlog {text_channel.id} delete general')
    #         await asyncio.sleep(10)
    #         async with bot.db_pool.acquire() as conn:
    #             rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2', first_channel.guild.id, text_channel.id)
    #         if not rows:
    #             report.pass_('mlog.delete.general')
    #         else:
    #             report.fail('mlog.delete.general', 'Rows still exist after delete')
    #     except Exception as e:
    #         report.error('mlog.delete.general', e)
    #     try:
    #         await first_channel.send(f'{prefix}mlog {text_channel.id} create channel {first_channel.id}')
    #         await asyncio.sleep(10)
    #         async with bot.db_pool.acquire() as conn:
    #             rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2', first_channel.guild.id, text_channel.id)
    #         if rows:
    #             report.pass_('mlog.create.channel')
    #         else:
    #             report.fail('mlog.create.channel', 'No rows found after create with channel')
    #     except Exception as e:
    #         report.error('mlog.create.channel', e)
    #     try:
    #         await first_channel.send(f'{prefix}mlog {text_channel.id} modify channel {second_channel.id}')
    #         await asyncio.sleep(10)
    #         async with bot.db_pool.acquire() as conn:
    #             rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2', first_channel.guild.id, text_channel.id)
    #         if rows and rows[0]['snowflakes'] and second_channel.id in rows[0]['snowflakes']:
    #             report.pass_('mlog.modify.channel')
    #         else:
    #             report.fail('mlog.modify.channel', f"Expected {second_channel.id} in snowflakes, got {rows[0]['snowflakes'] if rows else None}")
    #     except Exception as e:
    #         report.error('mlog.modify.channel', e)
    #     try:
    #         await first_channel.send(f'{prefix}mlog {text_channel.id} delete channel {first_channel.id}')
    #         await asyncio.sleep(10)
    #         async with bot.db_pool.acquire() as conn:
    #             rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2', first_channel.guild.id, text_channel.id)
    #         if not rows:
    #             report.pass_('mlog.delete.channel')
    #         else:
    #             report.fail('mlog.delete.channel', 'Rows still exist after delete')
    #     except Exception as e:
    #         report.error('mlog.delete.channel', e)
    #     try:
    #         await first_channel.send(f'{prefix}mlog {text_channel.id} create member {self_member.id}')
    #         await asyncio.sleep(10)
    #         async with bot.db_pool.acquire() as conn:
    #             rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2', first_channel.guild.id, text_channel.id)
    #         if rows:
    #             report.pass_('mlog.create.member')
    #         else:
    #             report.fail('mlog.create.member', 'No rows found after create with member')
    #     except Exception as e:
    #         report.error('mlog.create.member', e)
    #     try:
    #         await first_channel.send(f'{prefix}mlog {text_channel.id} modify member {dummy_member.id}')
    #         await asyncio.sleep(10)
    #         async with bot.db_pool.acquire() as conn:
    #             rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2', first_channel.guild.id, text_channel.id)
    #         if rows and rows[0]['snowflakes'] and dummy_member.id in rows[0]['snowflakes']:
    #             report.pass_('mlog.modify.member')
    #         else:
    #             report.fail('mlog.modify.member', f"Expected {dummy_member.id} in snowflakes, got {rows[0]['snowflakes'] if rows else None}")
    #     except Exception as e:
    #         report.error('mlog.modify.member', e)
    #     try:
    #         await first_channel.send(f'{prefix}mlog {text_channel.id} delete member {self_member.id}')
    #         await asyncio.sleep(10)
    #         async with bot.db_pool.acquire() as conn:
    #             rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2', first_channel.guild.id, text_channel.id)
    #         if not rows:
    #             report.pass_('mlog.delete.member')
    #         else:
    #             report.fail('mlog.delete.member', 'Rows still exist after delete')
    #     except Exception as e:
    #         report.error('mlog.delete.member', e)

    # async def rmv_command_test(self, channel: discord.abc.GuildChannel, prefix: Optional[str]):
    # async def smute_command_test(self, channel: discord.abc.GuildChannel, prefix: Optional[str]):
    # async def temps_command_test(self, channel: discord.abc.GuildChannel, prefix: Optional[str]):
    # async def xcap_command_test(self, channel: discord.abc.GuildChannel, prefix: Optional[str]):
    # async def xstage_command_test(self, channel: discord.abc.GuildChannel, prefix: Optional[str]):