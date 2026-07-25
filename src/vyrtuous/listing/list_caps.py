"""!/bin/python3
list_caps.py The purpose of this program is to list caps.

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
from vyrtuous.db.cap import Cap
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.messaging import emojis

MODEL = Cap


async def build_dictionary(
    guild_snowflake: int,
    obj,
) -> dict[int, dict[str, dict[int, dict[str, dict[str, int]]]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    caps = []
    dictionary: dict[int, dict[str, dict[int, dict[str, dict[str, int]]]]] = {}
    if isinstance(obj, discord.Guild):
        caps = await database_factory.select(guild_snowflake=obj.id, singular=False)
        guild_snowflake = obj.id
    elif isinstance(obj, discord.abc.GuildChannel):
        caps = await database_factory.select(
            channel_snowflake=obj.id, guild_snowflake=guild_snowflake, singular=False
        )
    else:
        caps = await database_factory.select(
            guild_snowflake=guild_snowflake, singular=False
        )
    if caps:
        for cap in caps:
            dictionary.setdefault(guild_snowflake, {"channels": {}})
            dictionary[guild_snowflake]["channels"].setdefault(
                cap.channel_snowflake, {"caps": {}}
            )
            dictionary[guild_snowflake]["channels"][cap.channel_snowflake]["caps"][
                cap.category
            ] = cap.duration_seconds
    return dictionary


async def build_pages(guild_snowflake: int, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        return "No caps found."

    duration_builder = DurationBuilder()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = guild.name
    if not isinstance(obj, int):
        obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Caps for {obj_name}"

    dictionary = await build_dictionary(guild_snowflake=guild_snowflake, obj=obj)
    processed_dictionary: list_service.CapDictionary = (
        await list_service.process_dictionary(
            cls=list_service.CapDictionary, dictionary=dictionary
        )
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
        for channel_snowflake, cap_dictionary in guild_data.get("channels", {}).items():
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                continue
            for moderation_type, duration_seconds in cap_dictionary.get(
                "caps", {}
            ).items():
                lines.append(
                    f"  ↳ {moderation_type} ({duration_builder.from_seconds(duration_seconds).as_str()})"
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
    if not pages:
        return "No caps found."
    return pages
