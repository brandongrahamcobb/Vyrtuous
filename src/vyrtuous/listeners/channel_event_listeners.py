"""!/bin/python3
channel_event_listeners.py A discord.py cog containing channel event listeners for the Vyrtuous bot.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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
"""

from datetime import datetime, timezone

import discord
from discord.ext import commands

from vyrtuous.aliases import unvoice_mute_alias_service, voice_mute_alias_service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.utils.channels import automute_channel_service, video_channel_service
from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    server_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.users import moderator_service


class ChannelEventListeners(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    @commands.Cog.listener()
    async def on_voice_state_update(
        self,
        member: discord.Member,
        before: discord.VoiceState,
        after: discord.VoiceState,
    ) -> None:
        if member.bot:
            return None
        self.__bot.registry.get(MemberState).active[member.id] = (
            member.display_name,
            datetime.now(timezone.utc),
        )
        if not after.channel:
            return None
        if before.channel == after.channel:
            if before.mute == after.mute:
                if before.self_mute == after.self_mute:
                    return None
                return None
        await ban_service.is_banned_then_kick_and_reset_cooldown(
            channel=after.channel, member=member
        )
        if await server_mute_service.is_server_muted(
            channel_snowflake=after.channel.id,
            guild_snowflake=after.channel.guild.id,
            member_snowflake=member.id,
        ):
            return None
        if video_channel_service.is_active_video_channel(channel=after.channel):
            await video_channel_service.update_video_channel_tasks(
                after=after, before=before, member=member
            )
        duration_value = "1h"
        if (
            after.channel.guild.id,
            member.id,
        ) in self.__bot.registry.get(MemberState).invincible:
            await unvoice_mute_alias_service.unvoice_mute(
                channel_snowflake=after.channel.id,
                guild_snowflake=after.channel.guild.id,
                member_snowflake=member.id,
                target="user",
            )
            await unvoice_mute_alias_service.log_unvoice_mute(
                author_snowflake=None,
                channel_snowflake=after.channel.id,
                display=True,
                guild_snowflake=after.channel.guild.id,
                is_channel_scope=True,
                member_snowflake=member.id,
                message_snowflake=None,
                message_channel_snowflake=None,
                target="user",
            )
            return None
        elif await automute_channel_service.is_active_automute_channel(
            channel_snowflake=after.channel.id,
            guild_snowflake=after.channel.guild.id,
        ):
            await voice_mute_alias_service.voice_mute(
                channel_snowflake=after.channel.id,
                duration_value="1h",
                guild_snowflake=after.channel.guild.id,
                member_snowflake=member.id,
                reason="Right-click automute",
                target="auto",
            )
            await voice_mute_alias_service.log_voice_mute(
                author_snowflake=None,
                channel_snowflake=after.channel.id,
                display=True,
                duration_value="1h",
                guild_snowflake=after.channel.guild.id,
                is_channel_scope=True,
                member_snowflake=member.id,
                message_snowflake=None,
                message_channel_snowflake=None,
                reason="Right-click automute",
                target="auto",
            )
        elif before.channel != after.channel:
            if await voice_mute_service.is_voice_muted(
                channel_snowflake=after.channel.id,
                guild_snowflake=after.channel.guild.id,
                member_snowflake=member.id,
            ):
                await voice_mute_alias_service.voice_mute(
                    channel_snowflake=after.channel.id,
                    duration_value=duration_value,
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                    target="user",
                    reason="Right-click muted",
                )
                await voice_mute_alias_service.log_voice_mute(
                    author_snowflake=None,
                    channel_snowflake=after.channel.id,
                    display=True,
                    duration_value="1h",
                    guild_snowflake=after.channel.guild.id,
                    is_channel_scope=True,
                    member_snowflake=member.id,
                    message_snowflake=None,
                    message_channel_snowflake=None,
                    reason="Right-click muted",
                    target="user",
                )
            elif after.mute:
                await unvoice_mute_alias_service.unvoice_mute(
                    channel_snowflake=after.channel.id,
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                    target="user",
                )
                await unvoice_mute_alias_service.log_unvoice_mute(
                    author_snowflake=None,
                    channel_snowflake=after.channel.id,
                    display=True,
                    guild_snowflake=after.channel.guild.id,
                    is_channel_scope=True,
                    member_snowflake=member.id,
                    message_snowflake=None,
                    message_channel_snowflake=None,
                    target="user",
                )
            await flag_service.warn(channel=after.channel, member=member)
        elif before.channel == after.channel:
            if before.mute and not after.mute:
                database_factory: DatabaseFactory = DatabaseFactory(VoiceMute)
                await database_factory.delete(
                    channel_snowflake=after.channel.id,
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                    target="user",
                )
            else:
                await voice_mute_alias_service.voice_mute(
                    channel_snowflake=after.channel.id,
                    duration_value=duration_value,
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                    target="user",
                    reason="Right-click muted",
                )
                await voice_mute_alias_service.log_voice_mute(
                    author_snowflake=None,
                    channel_snowflake=after.channel.id,
                    display=True,
                    duration_value="1h",
                    guild_snowflake=after.channel.guild.id,
                    is_channel_scope=True,
                    member_snowflake=member.id,
                    message_snowflake=None,
                    message_channel_snowflake=None,
                    reason="Right-click muted",
                    target="user",
                )


async def setup(bot: DiscordBot):
    await bot.add_cog(ChannelEventListeners(bot))
