"""cap.py The purpose of this program is to provide the Cap utility class.

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

from datetime import datetime
from typing import Optional

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.fields.duration import DurationObject
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    clean_guild_dictionary,
    flush_page,
)
from vyrtuous.utils.emojis import get_random_emoji


class Cap(DatabaseFactory):

    ACT = None
    CATEGORY = "cap"
    PLURAL = "Caps"
    SCOPES = ["channels"]
    SINGULAR = "Cap"
    UNDO = None

    REQUIRED_INSTANTIATION_ARGS = [
        "category",
        "channel_snowflake",
        "duration_seconds",
        "guild_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "updated_at",
    ]

    TABLE_NAME = "active_caps"

    def __init__(
        self,
        category: str,
        channel_snowflake: int,
        duration_seconds: int,
        guild_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
    ):
        self.category = category
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.duration_seconds = duration_seconds
        self.guild_snowflake = guild_snowflake
        self.updated_at = updated_at

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {Cap.PLURAL}"
        kwargs = object_dict.get("columns", None)

        caps = await Cap.select(**kwargs)
        for cap in caps:
            guild_dictionary.setdefault(cap.guild_snowflake, {"channels": {}})
            guild_dictionary[cap.guild_snowflake]["channels"].setdefault(
                cap.channel_snowflake, {"caps": {}}
            )
            guild_dictionary[cap.guild_snowflake]["channels"][cap.channel_snowflake][
                "caps"
            ][cap.category] = cap.duration_seconds

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
            for channel_snowflake, cap_dictionary in guild_data.get("channels").items():
                channel = guild.get_channel(channel_snowflake)
                for moderation_type, duration_seconds in cap_dictionary.get(
                    "caps", {}
                ).items():
                    lines.append(
                        f"  ↳ {moderation_type} ({DurationObject.from_seconds(duration_seconds)})"
                    )
                    field_count += 1
                    if field_count >= chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed, field_count = flush_page(embed, pages, title, guild.name)
                        lines = []
                if lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
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
        return pages
