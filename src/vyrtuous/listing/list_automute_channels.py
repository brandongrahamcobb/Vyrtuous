"""!/bin/python3stage"
automute_channel_service.py The purpose of this program is to extend Service to service the stage class.

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

from typing import Any

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = AutoMute


async def build_dictionary(
    guild_snowflake: int,
    obj,
) -> dict[int, dict[str, dict[int, dict[str, dict[str, Any]]]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    automutes = []
    dictionary: dict[int, dict[str, dict[int, dict[str, dict[str, Any]]]]] = {}
    if isinstance(obj, discord.Guild):
        automutes = await database_factory.select(
            guild_snowflake=obj.id, singular=False
        )
        guild_snowflake = obj.id
    elif isinstance(obj, discord.abc.GuildChannel):
        automutes = await database_factory.select(
            channel_snowflake=obj.id, guild_snowflake=guild_snowflake, singular=False
        )
    else:
        automutes = await database_factory.select(
            guild_snowflake=guild_snowflake, singular=False
        )
    if automutes:
        for automute in automutes:
            dictionary.setdefault(guild_snowflake, {"channels": {}})
            dictionary[guild_snowflake]["channels"].setdefault(
                automute.channel_snowflake, {}
            )
            dictionary[guild_snowflake]["channels"][
                automute.channel_snowflake
            ].setdefault("automutes", {})
            dictionary[guild_snowflake]["channels"][automute.channel_snowflake][
                "automutes"
            ].update({"expires_in": automute.expires_in})
    return dictionary


async def build_pages(guild_snowflake: int, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Automute Rooms for {obj_name}"

    dictionary = await build_dictionary(guild_snowflake=guild_snowflake, obj=obj)
    processed_dictionary: list_service.AutoMuteDictionary = (
        await list_service.process_dictionary(
            cls=list_service.AutoMuteDictionary, dictionary=dictionary
        )
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        automute_n = 0
        field_count = 0
        lines = []
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for channel_snowflake, automute_dictionary in guild_data.get(
            "channels", {}
        ).items():
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                continue
            lines.append(
                f"**Expires in:** {automute_dictionary.get('automutes', {}).get('expires_in', None)}"
            )
            automute_n += 1
            field_count += 1
            if field_count == list_service.CHUNK_SIZE:
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
        embed.description = f"**{original_description} ({automute_n})**"
        pages.append(embed)
    if not pages:
        return "No automute channels found."
    return pages
