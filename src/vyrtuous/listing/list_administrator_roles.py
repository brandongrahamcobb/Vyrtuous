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

from dataclasses import dataclass, field
from typing import Dict, List

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.administrator import AdministratorRole
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = AdministratorRole


@dataclass(frozen=True)
class AdministratorRoleDictionary:
    data: Dict[int, Dict[str, Dict[int, dict]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_roles: List[discord.Embed] = field(default_factory=list)


async def build_dictionary(obj) -> dict:
    database_factory = DatabaseFactory(MODEL)
    administrator_roles = []
    dictionary = {}
    if isinstance(obj, discord.Guild):
        administrator_roles = await database_factory.select(
            guild_snowflake=obj.id, singular=False
        )
    else:
        administrator_roles = await database_factory.select(singular=False)
    if administrator_roles:
        for administrator_role in administrator_roles:
            dictionary.setdefault(administrator_role.guild_snowflake, {"roles": {}})
            dictionary[administrator_role.guild_snowflake]["roles"].setdefault(
                administrator_role.role_snowflake, {}
            )
    return dictionary


async def build_pages(is_at_home: bool, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    pages: list[discord.Embed] = []

    obj_name = "All Servers"
    if obj:
        obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Administrator Roles in {obj_name}"

    full_dictionary = await build_dictionary(obj=obj)
    processed_dictionary = await list_service.process_dictionary(
        cls=AdministratorRoleDictionary, dictionary=full_dictionary
    )

    admin_role_n = 0
    for guild_snowflake, guild_data in processed_dictionary.data.items():
        field_count = 0
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for role_snowflake, _ in guild_data.get("roles", {}).items():
            role = guild.get_role(role_snowflake)
            if not role:
                continue
            if field_count >= list_service.CHUNK_SIZE:
                embed = list_service.flush_page(embed, pages, title, guild.name)
            embed.add_field(name=role.name, value=role.mention, inline=False)
            field_count += 1
            admin_role_n += 1
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({admin_role_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_guilds)
        pages.extend(processed_dictionary.skipped_roles)
    if not pages:
        return "No administrator roles found."
    return pages
