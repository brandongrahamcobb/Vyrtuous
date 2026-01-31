"""vegan.py The purpose of this program is to inherit from the DatabaseFactory to provide the vegan role.

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
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    clean_guild_dictionary,
    flush_page,
)
from vyrtuous.utils.emojis import get_random_emoji


class Vegan(DatabaseFactory):

    ACT = "vegan"
    CATEGORY = "vegan"
    PLURAL = "Vegans"
    SCOPES = ["member"]
    SINGULAR = "Vegan"
    UNDO = "vegan"

    REQUIRED_INSTANTIATION_ARGS = ["guild_snowflake", "member_snowflake"]
    OPTIONAL_ARGS = ["created_at", "updated_at"]

    TABLE_NAME = "vegans"

    def __init__(
        self,
        guild_snowflake: int,
        member_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.updated_at = updated_at

    @classmethod
    async def act_embed(cls, action_information, source, **kwargs):
        author = resolve_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        embed = discord.Embed(
            title=f"\U0001f525\U0001f525 {member.display_name} "
            f"is going Vegan!!!\U0001f525\U0001f525",
            description=(f"**By:** {author.mention}\n" f"**User:** {member.mention}\n"),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, action_information, source, **kwargs):
        author = resolve_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        embed = discord.Embed(
            title=f"\U0001f44e\U0001f44e "
            f"{member.display_name} is a Carnist \U0001f44e\U0001f44e",
            description=(f"**By:** {author.mention}\n" f"**User:** {member.mention}\n"),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def build_dictionary(cls, where_kwargs):
        dictionary = {}
        vegans = await Vegan.select(**where_kwargs)
        for vegan in vegans:
            dictionary.setdefault(vegan.guild_snowflake, {"members": {}})
            dictionary[vegan.guild_snowflake]["members"].setdefault(
                vegan.member_snowflake, {"vegans": {}}
            )
            dictionary[vegan.guild_snowflake]["members"][vegan.member_snowflake][
                "vegans"
            ].setdefault({"placeholder": {}})
        return dictionary
    
    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines, pages = 7, 0, [], []
        thumbnail = False
        title = f"{get_random_emoji()} {Vegan.PLURAL} {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"
        where_kwargs = object_dict.get("columns", None)

        guild_dictionary = await Vegan.build_dictionary(where_kwargs=where_kwargs)

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, vegan_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                else:
                    if not thumbnail:
                        embed.set_thumbnail(
                            url=object_dict.get("object", None).display_avatar.url
                        )
                        thumbnail = True
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
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
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )
        return pages
