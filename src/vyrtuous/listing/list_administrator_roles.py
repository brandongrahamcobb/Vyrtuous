"""!/bin/python3
list_administrator_role.py The purpose of this program is to list administrator roles.

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
from vyrtuous.db.administrator import AdministratorRole
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = AdministratorRole


async def build_dictionary(
    guild_snowflake: int,
) -> dict[int, dict[str, dict[int, dict]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    dictionary: dict[int, dict[str, dict[int, dict]]] = {}
    administrator_roles = await database_factory.select(
        guild_snowflake=guild_snowflake, singular=False
    )
    if administrator_roles:
        for administrator_role in administrator_roles:
            dictionary.setdefault(guild_snowflake, {"roles": {}})
            dictionary[guild_snowflake]["roles"].setdefault(
                administrator_role.role_snowflake, {}
            )
    return dictionary


async def build_pages(
    guild_snowflake: int,
) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        return "No administrator roles found."
    pages: list[discord.Embed] = []

    title = f"{emojis.get_random_emoji()} Administrator Roles in {guild.name}"

    full_dictionary = await build_dictionary(guild_snowflake=guild_snowflake)
    processed_dictionary: list_service.AdministratorRoleDictionary = (
        await list_service.process_dictionary(
            cls=list_service.AdministratorRoleDictionary, dictionary=full_dictionary
        )
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
    if not pages:
        return "No administrator roles found."
    return pages
