"""!/bin/python3
channel_event_listeners.py A discord.py cog containing channel event listeners for the Vyrtuous bot.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import time
from datetime import datetime, timezone

import discord
from discord.ext import commands

from vyrtuous.aliases import unvoice_mute_alias_service, voice_mute_alias_service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import ChannelState, MemberState
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.utils.channels import automute_channel_service, video_channel_service
from vyrtuous.utils.moderation import ban_service, flag_service, voice_mute_service
from vyrtuous.utils.users import vegan_service


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
        bot: DiscordBot = DiscordBot.get_instance()
        await ban_service.is_banned_then_kick_and_reset_cooldown(
            channel=after.channel, member=member
        )
        if video_channel_service.is_active_video_channel(channel=after.channel):
            await video_channel_service.update_video_channel_tasks(
                after=after, before=before, member=member
            )
        if invincible := self.__bot.registry.get(MemberState).invincible.get(
            after.channel.guild.id, None
        ):
            if member.id in invincible:
                database_factory: DatabaseFactory = DatabaseFactory(VoiceMute)
                await database_factory.delete(
                    channel_snowflake=after.channel.id,
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                )
                await member.edit(
                    mute=False,
                    reason="Channel event listeners unmute. Target: unknown.",
                )
                return None
        if before.channel != after.channel:
            if targets := await voice_mute_service.is_voice_muted(
                channel_snowflake=after.channel.id,
                guild_snowflake=after.channel.guild.id,
                member_snowflake=member.id,
                targets=["auto", "click", "command", "server"],
            ):
                for target in targets:
                    database_factory = DatabaseFactory(VoiceMute)
                    mute = await database_factory.select(
                        channel_snowflake=after.channel.id,
                        guild_snowflake=after.channel.guild.id,
                        target=target,
                        singular=True,
                    )
                    if mute:
                        duration_builder = DurationBuilder()
                        expires_in = duration_builder.from_timestamp(
                            mute.expires_in
                        ).build()
                        await member.edit(
                            mute=True,
                            reason="Channel event listeners remute. Target: command or server.",
                        )
                        await voice_mute_alias_service.log_voice_mute(
                            author_snowflake=None,
                            channel_snowflake=after.channel.id,
                            display=True,
                            duration=expires_in,
                            guild_snowflake=after.channel.guild.id,
                            is_channel_scope=True,
                            member_snowflake=member.id,
                            message_snowflake=None,
                            message_channel_snowflake=None,
                            reason=f"Channel event listeners mute. Target: {target}.",
                            target=target,
                        )
                        break
            elif await automute_channel_service.is_active_automute_channel(
                channel_snowflake=after.channel.id,
                guild_snowflake=after.channel.guild.id,
            ):
                database_factory = DatabaseFactory(AutoMute)
                automute = await database_factory.select(
                    channel_snowflake=after.channel.id,
                    guild_snowflake=after.channel.guild.id,
                    singular=True,
                )
                duration_builder = DurationBuilder()
                expires_in = duration_builder.from_timestamp(
                    automute.expires_in
                ).build()
                await voice_mute_alias_service.voice_mute(
                    channel_snowflake=after.channel.id,
                    duration=expires_in,
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                    target="auto",
                    reason=f"Automuted.",
                )
                await voice_mute_alias_service.log_voice_mute(
                    author_snowflake=None,
                    channel_snowflake=after.channel.id,
                    display=True,
                    duration=expires_in,
                    guild_snowflake=after.channel.guild.id,
                    is_channel_scope=True,
                    member_snowflake=member.id,
                    message_snowflake=None,
                    message_channel_snowflake=None,
                    reason=f"Channel event listeners automute. Target: auto.",
                    target="auto",
                )
            elif after.mute:
                await member.edit(
                    mute=False,
                    reason="Channel event listeners unmute. Target: unknown.",
                )
            await flag_service.warn(
                channel_snowflake=after.channel.id,
                guild_snowflake=after.channel.guild.id,
                member_snowflake=member.id,
            )
            await vegan_service.notify(channel=after.channel, member=member)
            joined_at = bot.registry.get(ChannelState).joined_at
            if before.channel:
                joined_at.pop(
                    (before.channel.id, before.channel.guild.id, member.id), None
                )
            joined_at[(after.channel.id, after.channel.guild.id, member.id)] = (
                time.time()
            )
        elif before.channel == after.channel:
            if before.mute and not after.mute:
                if await voice_mute_service.is_voice_muted(
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                    targets=["server"],
                ):
                    await member.edit(
                        mute=True,
                        reason="Channel event listeners remute. Target: server.",
                    )
                elif await voice_mute_service.is_voice_muted(
                    channel_snowflake=after.channel.id,
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                    targets=["command"],
                ):
                    await member.edit(
                        mute=True,
                        reason="Channel event listeners remute. Target: command.",
                    )
                    await voice_mute_service.alert_mute(
                        channel_snowflake=after.channel.id,
                        guild_snowflake=after.channel.guild.id,
                        member_snowflake=member.id,
                    )
                elif await voice_mute_service.is_voice_muted(
                    channel_snowflake=after.channel.id,
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                    targets=["auto"],
                ):
                    await unvoice_mute_alias_service.unvoice_mute(
                        channel_snowflake=after.channel.id,
                        guild_snowflake=after.channel.guild.id,
                        member_snowflake=member.id,
                        target="auto",
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
                        reason="Channel event listeners unmute. Target: auto.",
                        target="auto",
                    )
                elif await voice_mute_service.is_voice_muted(
                    channel_snowflake=after.channel.id,
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                    targets=["click"],
                ):
                    await unvoice_mute_alias_service.unvoice_mute(
                        channel_snowflake=after.channel.id,
                        guild_snowflake=after.channel.guild.id,
                        member_snowflake=member.id,
                        target="click",
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
                        reason="Channel event listeners unmute. Target: click.",
                        target="click",
                    )
            elif not before.mute and after.mute:
                if await automute_channel_service.is_active_automute_channel(
                    channel_snowflake=after.channel.id,
                    guild_snowflake=after.channel.guild.id,
                ):
                    database_factory = DatabaseFactory(AutoMute)
                    automute = await database_factory.select(
                        channel_snowflake=after.channel.id,
                        guild_snowflake=after.channel.guild.id,
                        singular=True,
                    )
                    duration_builder = DurationBuilder()
                    expires_in = duration_builder.from_timestamp(
                        automute.expires_in
                    ).build()
                    await voice_mute_alias_service.voice_mute(
                        channel_snowflake=after.channel.id,
                        duration=expires_in,
                        guild_snowflake=after.channel.guild.id,
                        member_snowflake=member.id,
                        target="auto",
                        reason="Automuted.",
                    )
                    await voice_mute_alias_service.log_voice_mute(
                        author_snowflake=None,
                        channel_snowflake=after.channel.id,
                        display=True,
                        duration=expires_in,
                        guild_snowflake=after.channel.guild.id,
                        is_channel_scope=True,
                        member_snowflake=member.id,
                        message_snowflake=None,
                        message_channel_snowflake=None,
                        reason="Channel event listeners mute. Target: auto.",
                        target="auto",
                    )
                else:
                    duration = DurationObject(number=1, prefix="", sign=1, unit="h")
                    await voice_mute_alias_service.voice_mute(
                        channel_snowflake=after.channel.id,
                        duration=duration,
                        guild_snowflake=after.channel.guild.id,
                        member_snowflake=member.id,
                        target="click",
                        reason="Right-click muted.",
                    )
                    await voice_mute_alias_service.log_voice_mute(
                        author_snowflake=None,
                        channel_snowflake=after.channel.id,
                        display=True,
                        duration=duration,
                        guild_snowflake=after.channel.guild.id,
                        is_channel_scope=True,
                        member_snowflake=member.id,
                        message_snowflake=None,
                        message_channel_snowflake=None,
                        reason="Channel event listeners mute. Target: click.",
                        target="click",
                    )


async def setup(bot: DiscordBot):
    await bot.add_cog(ChannelEventListeners(bot))
