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
from datetime import datetime, timedelta, timezone

import discord

from vyrtuous.alias.alias import Alias
from vyrtuous.base.record_service import RecordService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_channels,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger
from vyrtuous.video_room.video_room import VideoRoom


class VideoRoomService(RecordService):
    COOLDOWN = timedelta(minutes=30)
    cooldowns = {}
    lines, pages = [], []
    model = VideoRoom
    video_rooms = []
    video_tasks = {}

    @classmethod
    async def enforce_video(cls, member, channel, delay):
        await asyncio.sleep(delay)
        if not member.voice:
            return
        if member.voice.channel != channel:
            return
        if member.voice.self_video:
            return
        try:
            await member.move_to(None)
        except Exception as e:
            logger.info(f"Unable to enforce video by kicking the user. {e}")
        try:
            await member.send(
                f"{get_random_emoji()} You were kicked from {channel.mention} because your video feed stopped. {channel.mention} is a video-only channel."
            )
        except Exception as e:
            logger.info(f"Unable to send a message to enforce video. {e}")

    @classmethod
    def cancel_task(cls, key):
        task = cls.video_tasks.pop(key, None)
        if task:
            task.cancel()

    @classmethod
    async def enforce_video_message(cls, channel_snowflake, member_snowflake, message):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(channel_snowflake)
        now = datetime.now(timezone.utc)
        last_trigger = cls.cooldowns.get(member_snowflake, None)
        if last_trigger and now - last_trigger < cls.COOLDOWN:
            return
        cls.cooldowns[member_snowflake] = now
        await channel.send(message)

        async def reset_cooldown():
            await asyncio.sleep(cls.COOLDOWN.total_seconds())
            if cls.cooldowns.get(member_snowflake) == now:
                del cls.cooldowns[member_snowflake]

        asyncio.create_task(reset_cooldown())

    @classmethod
    async def reinforce_video_room(cls, member, before, after):
        if not after.channel:
            VideoRoomService.cancel_task((member.guild.id, member.id))
            return
        for video_room in VideoRoomService.video_rooms:
            if after.channel.id != video_room.channel_snowflake:
                continue
            if not after.self_video:
                if after.channel != before.channel:
                    if after.channel.permissions_for(
                        after.channel.guild.me
                    ).send_messages:
                        await VideoRoomService.enforce_video_message(
                            channel_snowflake=after.channel.id,
                            member_snowflake=member.id,
                            message=f"{get_random_emoji()} "
                            f"Hi {member.mention}, "
                            f"{after.channel.mention} is a video "
                            f"only room. You have 5 minutes to turn "
                            f"on your camera!",
                        )
            key = (member.guild.id, member.id)
            if before.channel != after.channel:
                VideoRoomService.cancel_task(key)
                if not after.self_video:
                    task = asyncio.create_task(
                        VideoRoomService.enforce_video(member, after.channel, 300)
                    )
                    VideoRoomService.video_tasks[key] = task
                break
            if before.self_video and not after.self_video:
                VideoRoomService.cancel_task(key)
                task = asyncio.create_task(
                    VideoRoomService.enforce_video(member, after.channel, 60)
                )
                VideoRoomService.video_tasks[key] = task
                break
            if not before.self_video and after.self_video:
                VideoRoomService.cancel_task(key)
                break

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        aliases = await Alias.select(singular=False, **where_kwargs)
        video_rooms = await VideoRoom.select(singular=False, **where_kwargs)
        for video_room in video_rooms:
            dictionary.setdefault(video_room.guild_snowflake, {"channels": {}})
            dictionary[video_room.guild_snowflake]["channels"].setdefault(
                video_room.channel_snowflake, {}
            )
            if aliases:
                for alias in aliases:
                    if (
                        alias.guild_snowflake == video_room.guild_snowflake
                        and alias.channel_snowflake == video_room.channel_snowflake
                    ):
                        dictionary[video_room.guild_snowflake]["channels"][
                            video_room.channel_snowflake
                        ].setdefault(alias.category, [])
                        dictionary[video_room.guild_snowflake]["channels"][
                            video_room.channel_snowflake
                        ][alias.category].append(alias.alias_name)
        skipped_channels = generate_skipped_channels(dictionary)
        skipped_guilds = generate_skipped_guilds(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_channels:
                VideoRoomService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_channels, title="Skipped Channels in Server,"
                    )
                )
            if skipped_guilds:
                VideoRoomService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Video Rooms"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await VideoRoomService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        vr_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                VideoRoomService.lines.append(f"Channel: {channel.mention}")
                field_count += 1
                for category, alias_names in channel_data.items():
                    VideoRoomService.lines.append(f"{category}")
                    field_count += 1
                    vr_n += 1
                    for name in alias_names:
                        VideoRoomService.lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= CHUNK_SIZE:
                            embed.add_field(
                                name="Information",
                                value="\n".join(VideoRoomService.lines),
                                inline=False,
                            )
                            embed = flush_page(
                                embed, VideoRoomService.pages, title, guild.name
                            )
                            VideoRoomService.lines = []
                            field_count = 0
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(VideoRoomService.lines),
                        inline=False,
                    )
                    embed = flush_page(embed, VideoRoomService.pages, title, guild.name)
                    field_count = 0
                    VideoRoomService.lines = []
            if VideoRoomService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(VideoRoomService.lines),
                    inline=False,
                )
            VideoRoomService.pages.append(embed)
        if VideoRoomService.pages:
            VideoRoomService.pages[0].description = f"**({vr_n})**"
        return VideoRoomService.pages

    @classmethod
    async def toggle_video_room(cls, channel_dict):
        kwargs = channel_dict.get("columns", None)
        video_room = await VideoRoom.select(**kwargs, singular=True)
        if video_room:
            action = "removed"
            VideoRoomService.video_rooms = [
                vr
                for vr in VideoRoomService.video_rooms
                if vr.channel_snowflake != video_room.channel_snowflake
            ]
            await VideoRoom.delete(**kwargs)
        else:
            video_room = VideoRoom(**kwargs)
            await video_room.create()
            VideoRoomService.video_rooms.append(video_room)
            action = "created"
        return f"Video-only room {action} in {channel_dict.get('mention', None)}."
