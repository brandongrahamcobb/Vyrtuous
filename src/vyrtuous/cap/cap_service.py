"""!/bin/python3
cap_service.py The purpose of this program is to extend Service to service the cap command class.

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

from copy import copy
from dataclasses import dataclass, field
from typing import Dict, List

import discord

from vyrtuous.cap.cap import Cap


@dataclass
class CapDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[str, int]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


class CapService:
    __CHUNK_SIZE = 7
    MODEL = Cap

    def __init__(
        self,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        duration_service=None,
        emoji=None,
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__duration_service = duration_service
        self.__emoji = emoji

    async def build_dictionary(self, where_kwargs):
        dictionary = {}
        caps = await self.__database_factory.select(singular=False, **where_kwargs)
        for cap in caps:
            dictionary.setdefault(cap.guild_snowflake, {"channels": {}})
            dictionary[cap.guild_snowflake]["channels"].setdefault(
                cap.channel_snowflake, {"caps": {}}
            )
            dictionary[cap.guild_snowflake]["channels"][cap.channel_snowflake]["caps"][
                cap.category
            ] = cap.duration_seconds
        return dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Caps"

        where_kwargs = object_dict.get("columns", None)

        dictionary = await self.build_dictionary(where_kwargs=where_kwargs)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=CapDictionary, dictionary=dictionary
        )
        if is_at_home:
            pages.extend(processed_dictionary.skipped_channels)
            pages.extend(processed_dictionary.skipped_guilds)

        cap_n = 0
        for guild_snowflake, guild_data in processed_dictionary.data.items():
            field_count = 0
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, cap_dictionary in guild_data.get("channels").items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                for moderation_type, duration_seconds in cap_dictionary.get(
                    "caps", {}
                ).items():
                    lines.append(
                        f"  ↳ {moderation_type} ({self.__duration_service.from_seconds(duration_seconds)})"
                    )
                    cap_n += 1
                    field_count += 1
                    if field_count >= self.__CHUNK_SIZE:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed = self.__dictionary_service.flush_page(
                            embed, pages, title, guild.name
                        )
                        lines = []
                        field_count = 0
                if lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
            pages.append(embed)
        if pages:
            pages[0].description = f"**({cap_n})**"
        return pages

    async def toggle_cap(self, category, channel_dict, hours):
        seconds = int(hours) * 3600
        where_kwargs = channel_dict.get("columns", None)
        where_kwargs.update({"category": category})
        cap = await self.__database_factory.select(singular=True, **where_kwargs)
        if cap and seconds:
            await self.__database_factory.update(
                set_kwargs={"duration_seconds": seconds}, where_kwargs=where_kwargs
            )
            return f"Cap `{category}` modified for {channel_dict.get('mention', None)}."
        elif cap:
            await self.__database_factory.delete(**where_kwargs)
            return (
                f"Cap of type {category} "
                f"and channel {channel_dict.get('mention', None)} deleted successfully."
            )
        else:
            where_kwargs.update({"duration_seconds": seconds})
            cap = self.MODEL(**where_kwargs)
            await self.__database_factory.create(cap)
            return (
                f"Cap `{category}` created for "
                f"{channel_dict.get('mention', None)} successfully."
            )

    async def assert_duration_exceeds_cap(self, category, duration, source_kwargs):
        exceeds_cap = False
        cap = await self.__database_factory.select(
            **source_kwargs, category=category, singular=True
        )
        duration_seconds = self.__duration_service.from_expires_in(
            duration.expires_in
        ).to_seconds()
        if cap:
            if duration_seconds > cap.duration_seconds:
                exceeds_cap = True
        else:
            self.__duration_service.duration = "8h"
            if duration_seconds > self.__duration_service.to_seconds():
                exceeds_cap = True
        return exceeds_cap

    async def migrate(self, kwargs):
        self.__database_factory.update(**kwargs)
