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

from dataclasses import dataclass, field
from typing import Dict, List

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.cap import Cap
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.messaging import emojis

MODEL = Cap


@dataclass
class CapDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[str, int]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


async def build_dictionary(obj):
    database_factory = DatabaseFactory(MODEL)
    caps = []
    dictionary = {}
    if isinstance(obj, discord.Guild):
        caps = await database_factory.select(guild_snowflake=obj.id, singular=False)
    elif isinstance(obj, discord.abc.GuildChannel):
        caps = await database_factory.select(channel_snowflake=obj.id, singular=False)
    else:
        caps = await database_factory.select(singular=False)
    if caps:
        for cap in caps:
            dictionary.setdefault(cap.guild_snowflake, {"channels": {}})
            dictionary[cap.guild_snowflake]["channels"].setdefault(
                cap.channel_snowflake, {"caps": {}}
            )
            dictionary[cap.guild_snowflake]["channels"][cap.channel_snowflake]["caps"][
                cap.category
            ] = cap.duration_seconds
    return dictionary


async def build_pages(is_at_home: bool, obj):
    bot = DiscordBot.get_instance()
    duration_builder = DurationBuilder()
    lines, pages = [], []

    obj_name = "All Servers"
    if obj:
        obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Caps for {obj_name}"

    dictionary = await build_dictionary(obj=obj)
    processed_dictionary = await list_service.process_dictionary(
        cls=CapDictionary, dictionary=dictionary
    )

    cap_n = 0
    for guild_snowflake, guild_data in processed_dictionary.data.items():
        field_count = 0
        lines = []
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for channel_snowflake, cap_dictionary in guild_data.get("channels").items():
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                continue
            for moderation_type, duration_seconds in cap_dictionary.get(
                "caps", {}
            ).items():
                lines.append(
                    f"  ↳ {moderation_type} ({duration_builder.from_seconds(duration_seconds).build(as_str=True)})"
                )
                cap_n += 1
                field_count += 1
                if field_count >= list_service.CHUNK_SIZE:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = list_service.flush_page(embed, pages, title, guild.name)
                    lines = []
                    field_count = 0
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({cap_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_channels)
        pages.extend(processed_dictionary.skipped_guilds)
    if not pages:
        return "No caps found."
    return pages
