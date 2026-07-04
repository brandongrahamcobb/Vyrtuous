"""!/bin/python3
streaming_service.py The purpose of this program is to extend Service service the stream command class.

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
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.stream import Stream
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = Stream


@dataclass
class StreamDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, List[int]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


async def build_dictionary(obj) -> dict:
    database_factory = DatabaseFactory(MODEL)
    streaming = []
    dictionary = {}
    if isinstance(obj, discord.Guild):
        streaming = await database_factory.select(
            guild_snowflake=obj.id, singular=False
        )
    elif isinstance(obj, discord.abc.GuildChannel):
        streaming = await database_factory.select(
            channel_snowflake=obj.id, singular=False
        )
    else:
        streaming = await database_factory.select(singular=False)
    if streaming:
        for stream in streaming:
            dictionary.setdefault(stream.guild_snowflake, {"channels": {}})
            dictionary[stream.guild_snowflake]["channels"].setdefault(
                stream.target_channel_snowflake, {"sources": []}
            )
            dictionary[stream.guild_snowflake]["channels"][
                stream.target_channel_snowflake
            ]["sources"].append(stream.source_channel_snowflake)
    return dictionary


async def build_pages(is_at_home: bool, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = "All Servers"
    if obj:
        obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Streaming Routes for {obj_name}"

    dictionary = await build_dictionary(obj=obj)
    processed_dictionary = await list_service.process_dictionary(
        cls=StreamDictionary, dictionary=dictionary
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        stream_n = 0
        field_count = 0
        lines = []
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for channel_snowflake, entry in guild_data.get("channels", {}).items():
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                continue
            sources = []
            for source_snowflake in entry["sources"]:
                source_channel = guild.get_channel(source_snowflake)
                if source_channel:
                    sources.append(source_channel.mention)
            if sources:
                source_lines = "\n".join(f"• {s}" for s in sources)
            else:
                source_lines = "• All sources"
            lines.append(f"**Target:** {channel.mention}\n**Sources:**\n{source_lines}")
            field_count += 1
            stream_n += 1
            if field_count >= list_service.CHUNK_SIZE:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
                embed = list_service.flush_page(embed, pages, title, guild.name)
                lines = []
                field_count = 0
        if lines:
            embed.add_field(
                name="Information",
                value="\n".join(lines),
                inline=False,
            )
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({stream_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_channels)
        pages.extend(processed_dictionary.skipped_guilds)
    if not pages:
        return "No streaming channels found."
    return pages
