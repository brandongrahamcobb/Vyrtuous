from copy import copy

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

from vyrtuous.video_room.video_room import VideoRoom


class VideoRoomService:
    __CHUNK_SIZE = 7
    __COOLDOWN = timedelta(minutes=30)
    MODEL = VideoRoom

    def __init__(
        self, *, bot=None, database_factory=None, dictionary_service=None, emoji=None
    ):
        self.cooldowns = {}
        self.video_rooms = []
        self.video_tasks = {}
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji

    async def build_clean_dictionary(self, is_at_home, where_kwargs):
        dictionary = {}
        video_rooms = await self.__database_factory.select(
            singular=False, **where_kwargs
        )
        for video_room in video_rooms:
            dictionary.setdefault(video_room.guild_snowflake, {"channels": {}})
            dictionary[video_room.guild_snowflake]["channels"].setdefault(
                video_room.channel_snowflake, {}
            )
        skipped_channels = self.__dictionary_service.generate_skipped_channels(
            dictionary
        )
        skipped_guilds = self.__dictionary_service.generate_skipped_guilds(dictionary)
        cleaned_dictionary = self.__dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        # if is_at_home:
        #     if skipped_channels:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_dict_pages(
        #                 skipped=skipped_channels, title="Skipped Channels in Server,"
        #             )
        #         )
        #     if skipped_guilds:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_set_pages(
        #                 skipped=skipped_guilds,
        #                 title="Skipped Servers",
        #             )
        #         )
        return cleaned_dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Video Rooms"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        vr_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                lines.append(f"Channel: {channel.mention}")
                field_count += 1
                for category, alias_names in channel_data.items():
                    lines.append(f"{category}")
                    field_count += 1
                    vr_n += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= self.__CHUNK_SIZE:
                            embed.add_field(
                                name="Information",
                                value="\n".join(lines),
                                inline=False,
                            )
                            embed = self.__dictionary_service.flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                            field_count = 0
                if field_count >= self.__CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = self.__dictionary_service.flush_page(
                        embed, pages, title, guild.name
                    )
                    field_count = 0
                    lines = []
            if lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)
        if pages:
            pages[0].description = f"**({vr_n})**"
        return pages

    async def toggle_video_room(self, channel_dict):
        kwargs = channel_dict.get("columns", None)
        video_room = await self.__database_factory.select(**kwargs, singular=True)
        if video_room:
            action = "removed"
            self.video_rooms = [
                vr
                for vr in self.video_rooms
                if vr.channel_snowflake != video_room.channel_snowflake
            ]
            await self.__database_factory.delete(**kwargs)
        else:
            video_room = self.MODEL(**kwargs)
            await self.__database_factory.create(video_room)
            self.video_rooms.append(video_room)
            action = "created"
        return f"Video-only room {action} in {channel_dict.get('mention', None)}."

    async def load_video_rooms_into_memory(self):
        self.video_rooms = self.__database_factory.select()

    def is_active_video_room(self, channel):
        for video_room in self.video_rooms:
            if video_room.channel_snowflake == channel.id:
                return True
        return False

    async def update_video_room_tasks(self, before, after, member):
        key = (member.guild.id, member.id)
        if before.channel and not after.channel:
            self.cancel_task(key)
            return
        if not before.channel and after.channel:
            if not after.self_video:
                await self.prompt_enable_camera(member, after.channel)
                self.schedule_enforcement(member, after.channel, delay=300)
            return
        if after.channel and before.self_video and not after.self_video:
            await self.prompt_enable_camera(member, after.channel)
            self.schedule_enforcement(member, after.channel, delay=300)
            return
        if after.channel and not before.self_video and after.self_video:
            self.cancel_task(key)

    def schedule_enforcement(self, member, channel, delay):
        key = (member.guild.id, member.id)
        self.cancel_task(key)
        task = asyncio.create_task(self.enforce_video(member, channel, delay))
        self.video_tasks[key] = task

    def cancel_task(self, key):
        task = self.video_tasks.pop(key, None)
        if task:
            task.cancel()

    async def prompt_enable_camera(self, member, channel):
        now = datetime.now(timezone.utc)
        last = self.cooldowns.get(member.id)
        if last and now - last < self.__COOLDOWN:
            return
        self.cooldowns[member.id] = now
        await channel.send(
            f"{self.__emoji.get_random_emoji()} "
            f"Hi {member.mention}, {channel.mention} is a video-only room. "
            f"You have 5 minutes to enable your camera."
        )

    async def enforce_video(self, member, channel, delay):
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
            self.__bot.logger.info(f"Video enforcement failed: {e}")
            return
        try:
            await member.send(
                f"{self.__emoji.get_random_emoji()} "
                f"You were removed from {channel.mention} because your camera was off. "
                f"There is a 30-minute cooldown before you can rejoin."
            )
        except Exception:
            pass
        self.cooldowns[member.id] = datetime.now(timezone.utc)
