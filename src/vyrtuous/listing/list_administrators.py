"""!/bin/python3
administrator_service.py The purpose of this program is to extend Service to service the administrator and administrator role classes.

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
from vyrtuous.db.administrator import Administrator
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = Administrator


async def build_dictionary(
    obj,
) -> dict[int, dict[str, dict[int, dict[str, dict[int, bool]]]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    administrators = []
    dictionary: dict[int, dict[str, dict[int, dict[str, dict[int, bool]]]]] = {}
    if isinstance(obj, discord.Guild):
        administrators = await database_factory.select(
            guild_snowflake=obj.id, singular=False
        )
    elif isinstance(obj, discord.Member):
        administrators = await database_factory.select(
            member_snowflake=obj.id, singular=False
        )
    else:
        administrators = await database_factory.select(singular=False)
    if administrators:
        for administrator in administrators:
            dictionary.setdefault(administrator.guild_snowflake, {"members": {}})
            dictionary[administrator.guild_snowflake]["members"].setdefault(
                administrator.member_snowflake, {"administrators": {}}
            )

            for role_snowflake in administrator.role_snowflakes:
                dictionary[administrator.guild_snowflake]["members"][
                    administrator.member_snowflake
                ]["administrators"].update({role_snowflake: True})
    return dictionary


async def build_pages(is_at_home: bool, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = "All Servers"
    if obj is not None and not isinstance(obj, (int, str)):
        obj_name = obj.name
    elif isinstance(obj, int):
        simplified_member = bot.registry.get(MemberState).active.get(obj, None)
        if simplified_member:
            obj_name = simplified_member[0]
        else:
            return "No administrators found."
    title = f"{emojis.get_random_emoji()} Administrators for {obj_name}"

    dictionary = await build_dictionary(obj=obj)
    processed_dictionary: list_service.AdministratorDictionary = (
        await list_service.process_dictionary(
            cls=list_service.AdministratorDictionary,
            dictionary=dictionary,
        )
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        admin_n = 0
        field_count = 0
        lines = []
        thumbnail = False
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for member_snowflake, member_dictionary in guild_data.get(
            "members", {}
        ).items():
            member = guild.get_member(member_snowflake)
            if member:
                if not thumbnail and isinstance(obj, discord.Member):
                    embed.set_thumbnail(url=obj.display_avatar.url)
                    thumbnail = True
                else:
                    lines.append(f"**User:** {member.display_name} {member.mention}")
            else:
                simplified_member = bot.registry.get(MemberState).active.get(
                    member_snowflake, None
                )
                if simplified_member:
                    display_name = simplified_member[0]
                    lines.append(f"**User:** {display_name} ({member_snowflake})")
            role_mentions = [
                role.mention
                for role_snowflake in member_dictionary.get("administrators", {})
                if (role := guild.get_role(role_snowflake))
            ]
            if isinstance(obj, discord.Member) and role_mentions:
                lines.append("\n**Roles:** " + "\n".join(role_mentions))
            admin_n += 1
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
        if lines:
            embed.add_field(
                name="Information",
                value="\n".join(lines),
                inline=False,
            )
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({admin_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_guilds)
        pages.extend(processed_dictionary.skipped_members)
    if not pages:
        return "No administrators found."
    return pages
