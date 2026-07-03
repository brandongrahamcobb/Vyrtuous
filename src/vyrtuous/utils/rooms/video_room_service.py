"""!/bin/python3
video_rooms_service.py The purpose of this program is to extend Service to service the video room class.

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

import asyncio
from datetime import timedelta

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import ChannelState, VideoRoomState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.video_room import VideoRoom
from vyrtuous.utils.messaging import emojis

MODEL = VideoRoom
COOLDOWN = timedelta(minutes=30)


async def toggle_video_room(channel):
    bot = DiscordBot.get_instance()
    video_rooms = bot.registry.get(ChannelState).video
    database_factory = DatabaseFactory(MODEL)
    video_room = await database_factory.select(
        channel_snowflake=channel.id, singular=True
    )
    if video_room:
        action = "removed"
        video_rooms.remove(video_room.channel_snowflake)
        await database_factory.delete(channel_snowflake=channel.id)
    else:
        video_room = MODEL(
            channel_snowflake=channel.id, guild_snowflake=channel.guild.id
        )
        await database_factory.create(video_room)
        video_rooms.add(video_room.channel_snowflake)
        action = "created"
    return f"Video-only room {action} in {channel.mention}."


async def populate():
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    original_set = bot.registry.get(ChannelState).video
    video_rooms = await database_factory.select(singular=False)
    for video_room in video_rooms:
        original_set.add(video_room.channel_snowflake)


def is_active_video_room(channel):
    bot = DiscordBot.get_instance()
    video_rooms = bot.registry.get(ChannelState).video
    if channel.id in video_rooms:
        return True
    return False


async def update_video_room_tasks(before, after, member):
    bot = DiscordBot.get_instance()
    key = (member.guild.id, member.id)
    video = bot.registry.get(VideoRoomState)

    if before.channel and not after.channel:
        video.cancel(key)
        return

    if not before.channel and after.channel:
        if not after.self_video:
            await _prompt_enable_camera(member, after.channel)
            video.schedule(
                member,
                after.channel,
                delay=300,
                coro=_enforce_video(member, after.channel, 300),
            )
        return

    if after.channel and before.self_video and not after.self_video:
        await _prompt_enable_camera(member, after.channel)
        video.schedule(
            member,
            after.channel,
            delay=300,
            coro=_enforce_video(member, after.channel, 300),
        )
        return

    if after.channel and not before.self_video and after.self_video:
        video.cancel(key)


async def _prompt_enable_camera(
    member: discord.Member, channel: discord.VoiceChannel
) -> None:
    bot = DiscordBot.get_instance()
    video = bot.registry.get(VideoRoomState)
    if video.is_on_cooldown(member.id, COOLDOWN):
        return
    video.set_cooldown(member.id)
    await channel.send(
        f"{emojis.get_random_emoji()} "
        f"Hi {member.mention}, {channel.mention} is a video-only room. "
        f"You have 5 minutes to enable your camera."
    )


async def _enforce_video(
    member: discord.Member, channel: discord.VoiceChannel, delay: int
) -> None:
    await asyncio.sleep(delay)
    bot = DiscordBot.get_instance()
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
    video = bot.registry.get(VideoRoomState)
    video.set_cooldown(member.id)
