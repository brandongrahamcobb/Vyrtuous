"""!/bin/python3
video_channels_service.py The purpose of this program is to extend Service to service the video channel class.

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

import asyncio
from datetime import timedelta

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import ChannelState, VideoChannelState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.video_channel import VideoChannel
from vyrtuous.utils.messaging import emojis

MODEL = VideoChannel
COOLDOWN = timedelta(minutes=30)


async def toggle_video_channel(channel_snowflake: int, guild_snowflake: int) -> str:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(str(channel_snowflake))
    video_channels = bot.registry.get(ChannelState).video
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    video_channel = await database_factory.select(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        singular=True,
    )
    if video_channel:
        action = "disabled"
        video_channels.remove(video_channel.channel_snowflake)
        await database_factory.delete(
            channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake
        )
    else:
        video_channel = MODEL(
            channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake
        )
        await database_factory.create(video_channel)
        video_channels.add(video_channel.channel_snowflake)
        action = "enabled"
    return f"Video-only channel {action} in {channel.mention}."


async def populate() -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    original_set = bot.registry.get(ChannelState).video
    video_channels = await database_factory.select(singular=False)
    for video_channel in video_channels:
        original_set.add(video_channel.channel_snowflake)


def is_active_video_channel(channel) -> bool:
    bot: DiscordBot = DiscordBot.get_instance()
    video_channels = bot.registry.get(ChannelState).video
    if channel.id in video_channels:
        return True
    return False


async def update_video_channel_tasks(before, after, member) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    key = (member.guild.id, member.id)
    video = bot.registry.get(VideoChannelState)

    if before.channel and not after.channel:
        video.cancel(key)
        return

    if not before.channel and after.channel:
        if not after.self_video:
            await _prompt_enable_camera(member, after.channel)
            video.schedule(
                member,
                after.channel,
                coro=_enforce_video(member, after.channel, 300),
            )
        return

    if after.channel and before.self_video and not after.self_video:
        await _prompt_enable_camera(member, after.channel)
        video.schedule(
            member,
            after.channel,
            coro=_enforce_video(member, after.channel, 300),
        )
        return

    if after.channel and not before.self_video and after.self_video:
        video.cancel(key)


async def _prompt_enable_camera(
    member: discord.Member, channel: discord.VoiceChannel
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    video = bot.registry.get(VideoChannelState)
    if video.is_on_cooldown(member.id, COOLDOWN):
        return
    video.set_cooldown(member.id)
    await channel.send(
        f"{emojis.get_random_emoji()} "
        f"Hi {member.mention}, {channel.mention} is a video-only channel. "
        f"You have 5 minutes to enable your camera."
    )


async def _enforce_video(
    member: discord.Member, channel: discord.VoiceChannel, delay: int
) -> None:
    await asyncio.sleep(delay)
    bot: DiscordBot = DiscordBot.get_instance()
    if not member.voice or member.voice.channel != channel or member.voice.self_video:
        return
    try:
        await member.move_to(None)
    except Exception as e:
        bot.logger.info(f"Video enforcement failed: {e}")
        return
    try:
        await member.send(
            f"{emojis.get_random_emoji()} "
            f"You were removed from {channel.mention} because your camera was off. "
            f"There is a 30-minute cooldown before you can rejoin."
        )
    except Exception:
        pass
    video = bot.registry.get(VideoChannelState)
    video.set_cooldown(member.id)
