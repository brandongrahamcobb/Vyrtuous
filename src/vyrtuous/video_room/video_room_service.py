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
            self.__bot.logger.info(f"Unable to enforce video by kicking the user. {e}")
        try:
            await member.send(
                f"{self.__emoji.get_random_emoji()} You were kicked from {channel.mention} because your video feed stopped. {channel.mention} is a video-only channel."
            )
        except Exception as e:
            self.__bot.logger.info(f"Unable to send a message to enforce video. {e}")

    def cancel_task(self, key):
        task = self.video_tasks.pop(key, None)
        if task:
            task.cancel()

    async def enforce_video_message(self, channel_snowflake, member_snowflake, message):
        channel = self.__bot.get_channel(channel_snowflake)
        now = datetime.now(timezone.utc)
        last_trigger = self.cooldowns.get(member_snowflake, None)
        if last_trigger and now - last_trigger < self.__COOLDOWN:
            return
        self.cooldowns[member_snowflake] = now
        await channel.send(message)

        async def reset_cooldown():
            await asyncio.sleep(self.__COOLDOWN.total_seconds())
            if self.cooldowns.get(member_snowflake) == now:
                del self.cooldowns[member_snowflake]

        asyncio.create_task(reset_cooldown())

    async def reinforce_video_room(self, member, before, after):
        if not after.channel:
            self.cancel_task((member.guild.id, member.id))
            return
        for video_room in self.video_rooms:
            if after.channel.id != video_room.channel_snowflake:
                continue
            if not after.self_video:
                if after.channel != before.channel:
                    if after.channel.permissions_for(
                        after.channel.guild.me
                    ).send_messages:
                        await self.enforce_video_message(
                            channel_snowflake=after.channel.id,
                            member_snowflake=member.id,
                            message=f"{self.__emoji.get_random_emoji()} "
                            f"Hi {member.mention}, "
                            f"{after.channel.mention} is a video "
                            f"only room. You have 5 minutes to turn "
                            f"on your camera!",
                        )
            key = (member.guild.id, member.id)
            if before.channel != after.channel:
                self.cancel_task(key)
                if not after.self_video:
                    task = asyncio.create_task(
                        self.enforce_video(member, after.channel, 300)
                    )
                    self.video_tasks[key] = task
                break
            if before.self_video and not after.self_video:
                self.cancel_task(key)
                task = asyncio.create_task(
                    self.enforce_video(member, after.channel, 60)
                )
                self.video_tasks[key] = task
                break
            if not before.self_video and after.self_video:
                self.cancel_task(key)
                break

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
