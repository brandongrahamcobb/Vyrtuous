"""!/bin/python3
list_video_channels.py The purpose of this program is to list video-only channels.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.video_channel import VideoChannel
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = VideoChannel


async def build_dictionary(obj) -> dict[int, dict[str, dict[int, dict]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    video_channels = []
    dictionary: dict[int, dict[str, dict[int, dict]]] = {}
    if isinstance(obj, discord.Guild):
        video_channels = await database_factory.select(
            guild_snowflake=obj.id, singular=False
        )
    elif isinstance(obj, discord.abc.GuildChannel):
        video_channels = await database_factory.select(
            channel_snowflake=obj.id, singular=False
        )
    else:
        video_channels = await database_factory.select(singular=False)
    if video_channels:
        for video_channel in video_channels:
            dictionary.setdefault(video_channel.guild_snowflake, {"channels": {}})
            dictionary[video_channel.guild_snowflake]["channels"][
                video_channel.channel_snowflake
            ] = {}
    return dictionary


async def build_pages(obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Video Rooms in {obj_name}"

    dictionary = await build_dictionary(obj=obj)
    processed_dictionary: list_service.VideoChannelDictionary = (
        await list_service.process_dictionary(
            cls=list_service.VideoChannelDictionary, dictionary=dictionary
        )
    )

    vc_n = 0
    for guild_snowflake, guild_data in processed_dictionary.data.items():
        field_count = 0
        lines = []
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for channel_snowflake, channel_data in guild_data.get("channels", {}).items():
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                continue
            lines.append(f"Channel: {channel.mention}")
            vc_n += 1
            field_count += 1
            for category, alias_names in channel_data.items():
                lines.append(f"{category}")
                field_count += 1
                for name in alias_names:
                    lines.append(f"  ↳ {name}")
                    field_count += 1
                    if field_count >= list_service.CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed = list_service.flush_page(embed, pages, title, guild.name)
                        lines = []
                        field_count = 0
            if field_count >= list_service.CHUNK_SIZE:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
                embed = list_service.flush_page(embed, pages, title, guild.name)
                field_count = 0
                lines = []
        if lines:
            embed.add_field(
                name="Information",
                value="\n".join(lines),
                inline=False,
            )
        pages.append(embed)
        original_description = embed.description or ""
        embed.description = f"**{original_description}** **({vc_n})**"
    if not pages:
        return "No video channels found."
    return pages
