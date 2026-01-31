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

from datetime import datetime, timezone


import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.fields.duration import DurationObject
from vyrtuous.utils.dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    clean_dictionary,
    flush_page,
)
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.inc.helpers import CHUNK_SIZE


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
    lines, pages = [], []

    def __init__(
        self,
        category: str,
        channel_snowflake: int,
        duration_seconds: int,
        guild_snowflake: int,
        created_at: datetime = datetime.now(timezone.utc),
        updated_at: datetime = datetime.now(timezone.utc),
    ):
        self.category = category
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.duration_seconds = duration_seconds
        self.guild_snowflake = guild_snowflake
        self.updated_at = updated_at

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        caps = await Cap.select(**where_kwargs)
        for cap in caps:
            dictionary.setdefault(cap.guild_snowflake, {"channels": {}})
            dictionary[cap.guild_snowflake]["channels"].setdefault(
                cap.channel_snowflake, {"caps": {}}
            )
            dictionary[cap.guild_snowflake]["channels"][cap.channel_snowflake]["caps"][
                cap.category
            ] = cap.duration_seconds
        skipped_channels = generate_skipped_channels(dictionary)
        skipped_guilds = generate_skipped_guilds(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_channels:
                Cap.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
            if skipped_guilds:
                Cap.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} {Cap.PLURAL}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await Cap.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
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
                    Cap.lines.append(
                        f"  ↳ {moderation_type} ({DurationObject.from_seconds(duration_seconds)})"
                    )
                    field_count += 1
                    if field_count >= CHUNK_SIZE:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(Cap.lines),
                            inline=False,
                        )
                        embed = flush_page(embed, Cap.pages, title, guild.name)
                        Cap.lines = []
                        field_count = 0
                if Cap.lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(Cap.lines),
                        inline=False,
                    )
            Cap.pages.append(embed)
        return Cap.pages

    @classmethod
    async def toggle_cap(cls, category, channel_dict, hours):
        seconds = int(hours) * 3600
        where_kwargs = channel_dict.get("columns", None)
        where_kwargs.update({"category": category})
        cap = await Cap.select(**where_kwargs, singular=True)
        if cap and seconds:
            await Cap.update(
                set_kwargs={"duration_seconds": seconds}, where_kwargs=where_kwargs
            )
            return f"Cap `{category}` modified for {channel_dict.get("mention", None)}."
        elif cap:
            await Cap.delete(**where_kwargs)
            return (
                f"Cap of type {category} "
                f"and channel {channel_dict.get("mention", None)} deleted successfully."
            )
        else:
            where_kwargs.update({"duration_seconds": seconds})
            cap = Cap(**where_kwargs)
            await cap.create()
            return (
                f"Cap `{category}` created for "
                f"{channel_dict.get("mention", None)} successfully."
            )
