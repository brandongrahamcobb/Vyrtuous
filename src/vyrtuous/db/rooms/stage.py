"""stage.py The purpose of this program is to provide the Stage utility class.

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
import time

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.fields.duration import DurationObject
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    clean_guild_dictionary,
    flush_page,
)


class Stage(DatabaseFactory):

    ACT = "stage"
    CATEGORY = "stage"
    PLURAL = "Stages"
    SCOPES = ["channels"]
    SINGULAR = "Stage"
    UNDO = "stage"

    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "expires_in",
        "guild_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "expired",
        "updated_at",
    ]

    TABLE_NAME = "active_stages"

    def __init__(
        self,
        channel_snowflake: int,
        expires_in: Optional[datetime],
        guild_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.expires_in = expires_in
        self.guild_snowflake = guild_snowflake
        self.updated_at = updated_at

    async def send_stage_ask_to_speak_message(
        self, join_log: dict[int, discord.Member], member: discord.Member
    ):
        bot = DiscordBot.get_instance()
        now = time.time()
        join_log[member.id] = [t for t in join_log[member.id] if now - t < 300]
        if len(join_log[member.id]) < 1:
            join_log[member.id].append(now)
            embed = discord.Embed(
                title=f"{get_random_emoji()} {self.channel_snowflake} — Stage Mode",
                description=f"Ends <t:{int(self.expires_in.timestamp())}:R>",
                color=discord.Color.green(),
            )
            embed.add_field(name="\u200b", value="**Ask to speak!**", inline=False)
            await bot.get_channel(self.channel_snowflake).send(embed=embed)

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} Stages"
        kwargs = object_dict.get("columns", None)

        stages = await Stage.select(**kwargs)

        for stage in stages:
            guild_dictionary.setdefault(stage.guild_snowflake, {"channels": {}})
            guild_dictionary[stage.guild_snowflake]["channels"].setdefault(
                stage.channel_snowflake, {}
            )
            guild_dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ].setdefault("stages", {})
            guild_dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ]["stages"].update(
                {"expires_in": DurationObject.from_expires_in(stage.expires_in)}
            )

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
            for channel_snowflake, stage_dictionary in guild_data.get(
                "channels"
            ).items():
                channel = guild.get_channel(channel_snowflake)
                lines.append(
                    f"**Expires in:** {stage_dictionary.get("expires_in", None)}"
                )
                field_count += 1
                if field_count == chunk_size:
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
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
        return pages
