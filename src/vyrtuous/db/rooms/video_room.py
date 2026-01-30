"""video_rooms.py A utility module for managing video rooms in the Vyrtuous Discord bot.

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

from datetime import datetime, timedelta, timezone
from typing import Optional
import asyncio

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.utils.logger import logger
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    clean_guild_dictionary,
    flush_page,
)


class VideoRoom(DatabaseFactory):

    COOLDOWN = timedelta(minutes=30)
    cooldowns = {}
    video_rooms = []
    video_tasks = {}

    ACT = "vr"
    CATEGORY = "vr"
    PLURAL = "Video Rooms"
    SCOPES = ["channels"]
    SINGULAR = "Video Rooms"
    UNDO = "vr"

    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "updated_at",
    ]

    TABLE_NAME = "video_rooms"

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_mention = f"<#{channel_snowflake}>"
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.is_video_room: Optional[bool] = True
        self.updated_at = updated_at

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
            VideoRoom.cancel_task((member.guild.id, member.id))
            return
        for video_room in VideoRoom.video_rooms:
            if after.channel.id != video_room.channel_snowflake:
                continue
            if not after.self_video:
                if after.channel != before.channel:
                    if after.channel.permissions_for(
                        after.channel.guild.me
                    ).send_messages:
                        await VideoRoom.enforce_video_message(
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
                VideoRoom.cancel_task(key)
                if not after.self_video:
                    task = asyncio.create_task(
                        VideoRoom.enforce_video(member, after.channel, 300)
                    )
                    VideoRoom.video_tasks[key] = task
                break
            if before.self_video and not after.self_video:
                VideoRoom.cancel_task(key)
                task = asyncio.create_task(
                    VideoRoom.enforce_video(member, after.channel, 60)
                )
                VideoRoom.video_tasks[key] = task
                break
            if not before.self_video and after.self_video:
                VideoRoom.cancel_task(key)
                break

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {VideoRoom.PLURAL}"

        kwargs = object_dict.get("columns", None)

        aliases = await Alias.select(**kwargs)
        video_rooms = await VideoRoom.select(**kwargs)

        for video_room in video_rooms:
            guild_dictionary.setdefault(video_room.guild_snowflake, {"channels": {}})
            guild_dictionary[video_room.guild_snowflake]["channels"].setdefault(
                video_room.channel_snowflake, {}
            )
            if aliases:
                for alias in aliases:
                    if (
                        alias.guild_snowflake == video_room.guild_snowflake
                        and alias.channel_snowflake == video_room.channel_snowflake
                    ):
                        guild_dictionary[video_room.guild_snowflake]["channels"][
                            video_room.channel_snowflake
                        ].setdefault(alias.category, [])
                        guild_dictionary[video_room.guild_snowflake]["channels"][
                            video_room.channel_snowflake
                        ][alias.category].append(alias.alias_name)

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                lines.append(f"Channel: {channel.mention}")
                field_count += 1
                for category, alias_names in channel_data.items():
                    lines.append(f"{category}")
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= chunk_size:
                            embed.add_field(
                                name="Information",
                                value="\n".join(lines),
                                inline=False,
                            )
                            embed, field_count = flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if is_at_home:
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )

    @classmethod
    async def toggle_video_room(cls, channel_dict):
        kwargs = channel_dict.get("columns", None)
        video_room = await VideoRoom.select(**kwargs, singular=True)
        if video_room:
            action = "removed"
            VideoRoom.video_rooms = [
                vr
                for vr in VideoRoom.video_rooms
                if vr.channel_snowflake != video_room.channel_snowflake
            ]
            await VideoRoom.delete(**kwargs)
        else:
            video_room = VideoRoom(**kwargs)
            await video_room.create()
            VideoRoom.video_rooms.append(video_room)
            action = "created"
        return f"Video-only room {action} in {channel_dict.get("mention", None)}."
