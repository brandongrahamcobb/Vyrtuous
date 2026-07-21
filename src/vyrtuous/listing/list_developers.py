"""!/bin/python3
developer_service.py The purpose of this program is to extend Service to service the developer class.

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
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.developer import Developer
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = Developer


async def build_dictionary(obj) -> dict[str, dict[int, dict[str, dict[str, str]]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    developers = []
    dictionary: dict[str, dict[int, dict[str, dict[str, str]]]] = {}
    if isinstance(obj, discord.Member):
        developers = await database_factory.select(
            member_snowflake=obj.id, singular=False
        )
    else:
        developers = await database_factory.select(singular=False)
    if developers:
        for developer in developers:
            dictionary.setdefault("members", {})
            dictionary["members"].setdefault(
                developer.member_snowflake, {"developers": {}}
            )
            dictionary["members"][developer.member_snowflake]["developers"].update(
                {"placeholder": "placeholder"}
            )
    return dictionary


async def build_pages(obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = "All Members"
    if obj is not None and not isinstance(obj, int):
        obj_name = obj.name
    elif isinstance(obj, int):
        simplified_member = bot.registry.get(MemberState).active.get(obj, None)
        if simplified_member:
            obj_name = simplified_member[0]
        else:
            return "This command must target a valid member."
    title = f"{emojis.get_random_emoji()} Developers for {obj_name}"

    dictionary = await build_dictionary(obj=obj)

    embed = discord.Embed(
        title=title, description="Information", color=discord.Color.blue()
    )
    dev_n = 0
    for _, values in dictionary.items():
        field_count = 0
        thumbnail = False
        for member_snowflake, _ in values.items():
            user = bot.get_user(member_snowflake)
            if user:
                if not isinstance(obj, discord.Member):
                    lines.append(f"**User:** {user.display_name} {user.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(url=obj.display_avatar.url)
                    thumbnail = True
            else:
                simplified_member = bot.registry.get(MemberState).active.get(
                    member_snowflake, None
                )
                if simplified_member:
                    display_name = simplified_member[0]
                    lines.append(f"**User:** {display_name} {member_snowflake}")
                else:
                    continue
            dev_n += 1
            if field_count >= list_service.CHUNK_SIZE:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
                pages.append(embed)
                embed = discord.Embed(
                    title=title,
                    description=title,
                    color=discord.Color.blue(),
                )
                field_count = 0
                lines = []
        if lines:
            embed.add_field(
                name="Information",
                value="\n".join(lines),
                inline=False,
            )
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({dev_n})**"
        pages.append(embed)
    if not pages:
        return "No developers found."
    return pages
