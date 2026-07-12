"""!/bin/python3
permission_service.py The purpose of this program is to provide the service for deciding whether a member has sufficient permissions.

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

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis


def build_dictionary() -> dict[int, dict[str, dict[int, bool]]]:
    bot: DiscordBot = DiscordBot.get_instance()
    dictionary: dict[int, dict[str, dict[int, bool]]] = {}
    for (
        guild_snowflake,
        member_snowflakes,
    ) in bot.registry.get(MemberState).invincible.items():
        dictionary.setdefault(guild_snowflake, {"members": {}})
        for member_snowflake in member_snowflakes:
            dictionary[guild_snowflake]["members"][member_snowflake] = True
    return dictionary


async def build_pages(is_at_home: bool, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = "All Servers"
    if obj is not None and not isinstance(obj, (int, str)):
        obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Heroes for {obj_name}"

    dictionary = build_dictionary()
    processed_dictionary: list_service.HeroDictionary = (
        await list_service.process_dictionary(
            cls=list_service.HeroDictionary, dictionary=dictionary
        )
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        hero_n = 0
        field_count = 0
        lines = []
        thumbnail_set = False
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for member_snowflake, _ in guild_data.get("members", {}).items():
            member = guild.get_member(member_snowflake)
            if member is None:
                continue
            hero_n += 1
            if not isinstance(obj, discord.Member):
                lines.append(f"**User:** {member.display_name} {member.mention}")
                field_count += 1
            elif not thumbnail_set:
                embed.set_thumbnail(url=member.display_avatar.url)
                thumbnail_set = True
            field_count += 1
            if field_count >= list_service.CHUNK_SIZE:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
                embed = list_service.flush_page(embed, pages, title, guild.name)
                lines = []
                field_count = 0
        if lines:
            embed.add_field(name="Information", value="\n".join(lines), inline=False)
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({hero_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_guilds)
        pages.extend(processed_dictionary.skipped_members)
    return pages
