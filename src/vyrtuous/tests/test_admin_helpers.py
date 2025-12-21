
''' test_admin_helpers.py The purpose of this program is to provide the shared test variables for tests using Discord objects.
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
from typing import Optional
from unittest.mock import PropertyMock, patch
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.tests.make_mock_objects import make_mock_member, make_mock_message
from vyrtuous.tests.test_suite import make_capturing_send
import discord

async def admin_cleanup(guild_id: Optional[int], privileged_author_id: Optional[int]):
    bot = DiscordBot.get_instance()
    async with bot.db_pool.acquire() as conn:
        await conn.execute('''
            UPDATE users SET administrator_guild_ids=$2 WHERE discord_snowflake=$1
        ''', int(privileged_author_id), [int(guild_id)])

async def admin_initiation(guild_id: Optional[int], privileged_author_id: Optional[int]):
    bot = DiscordBot.get_instance()
    async with bot.db_pool.acquire() as conn:
        await conn.execute('''
            INSERT INTO users (discord_snowflake, developer_guild_ids, updated_at, created_at)
            VALUES ($1, $2, NOW(), NOW())
            ON CONFLICT (discord_snowflake) 
            DO UPDATE SET developer_guild_ids = $2, updated_at = NOW()
        ''', int(privileged_author_id), [int(guild_id)])
    
async def prepared_command_handling(author, bot, channel, content, data, guild, prefix):
    mock_message = make_mock_message(allowed_mentions=True, author=author, channel=channel, content=content, data=data, embeds=[], guild=guild, id=MESSAGE_ID)
    view = cmd_view.StringView(f"{prefix}{mock_message.content}")
    view.skip_string(prefix)
    mock_bot_user = make_mock_member(id=PRIVILEGED_AUTHOR_ID, name=PRIVILEGED_AUTHOR_NAME)
    with patch.object(type(bot), "user", new_callable=PropertyMock) as mock_user:
        mock_user.return_value = mock_bot_user
        view = cmd_view.StringView(mock_message.content)
        view.skip_string(prefix)
        ctx = Context(
            bot=bot,
            message=mock_message,
            prefix=prefix,
            view=view
        )
        command_name = content.lstrip(prefix).split()[0]
        ctx.command = bot.get_command(command_name)
        ctx.send = channel.send.__get__(channel, type(channel))
        ctx.invoked_with = command_name
        view.skip_string(command_name)
        view.skip_ws()
        fake_channels = {}
        if channel.type == discord.ChannelType.text:
            fake_channels = {
                guild.id: [
                    {"channel_id": channel.id, "channel_name": channel.name, "enabled": False}
                ]
            }
        cog_instance = bot.get_cog("AdminCommands")
        capturing_send = make_capturing_send(channel, author)
        cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
        with patch.object(cog_instance.channel_service, "resolve_channel", return_value=channel):
            def mock_isinstance(obj, cls):
                if cls == discord.VoiceChannel:
                    return hasattr(obj, 'type') and obj.type == discord.ChannelType.voice
                elif cls == discord.TextChannel:
                    return hasattr(obj, 'type') and obj.type == discord.ChannelType.text
                else:
                    return isinstance(obj, cls)
            with patch("vyrtuous.cogs.admin_commands.isinstance", side_effect=mock_isinstance):
                if channel.type == discord.ChannelType.text:
                    with patch("vyrtuous.utils.statistics.Statistics.get_statistic_voice_channels", return_value=fake_channels):
                        await bot.invoke(ctx)
                else:
                    await bot.invoke(ctx)