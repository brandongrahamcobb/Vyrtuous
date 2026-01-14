"""video_rooms.py A utility module for managing video rooms in the Vyrtuous Discord bot.

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
"""

from datetime import datetime, timedelta, timezone
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.rooms.room import Room
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.setup_logging import logger
import asyncio

class VideoRoom(Room):

    ACT = "vr"
    COOLDOWN = timedelta(minutes=30)
    cooldowns = {}
    PLURAL = "Video Rooms"
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
    video_rooms = []
    video_tasks = {}

    def __init__(
        self, channel_snowflake: Optional[int], guild_snowflake: Optional[int], created_at: Optional[datetime] = None, updated_at: Optional[datetime] = None
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
        last_trigger = cls.cooldowns.get(member_snowflake)
        if last_trigger and now - last_trigger < cls.COOLDOWN:
            return
        cls.cooldowns[member_snowflake] = now
        await channel.send(message)
        async def reset_cooldown():
            await asyncio.sleep(cls.COOLDOWN.total_seconds())
            if cls.cooldowns.get(member_snowflake) == now:
                del cls.cooldowns[member_snowflake]
        asyncio.create_task(reset_cooldown())
