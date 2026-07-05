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

from vyrtuous.db.administrator import AdministratorRole
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.users import administrator_service

MODEL = AdministratorRole


async def is_added_role_administrator(
    guild_snowflake: int, role_snowflake: int
) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    where_kwargs = {
        "guild_snowflake": int(guild_snowflake),
        "role_snowflake": int(role_snowflake),
    }
    administrator_roles = await database_factory.select(
        singular=False,
        **where_kwargs,
    )
    if not administrator_roles:
        return False
    return True


async def toggle_administrator_role(role) -> list[discord.Embed]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    title = f"{emojis.get_random_emoji()} Administrators and Roles"
    administrators = await administrator_service.administrators_by_role(
        role_snowflake=role.id
    )
    administrator_roles = await is_added_role_administrator(
        guild_snowflake=role.guild.id, role_snowflake=role.id
    )
    if administrator_roles:
        action = "revoked"
        if administrator_roles:
            await database_factory.delete(role_snowflake=role.id)
        revoked_members: dict[int, dict[int, list[discord.Member]]] = {}
        for administrator in administrators:
            member = role.guild.get_member(administrator.member_snowflake)
            await administrator_service.removed_role(
                guild_snowflake=administrator.guild_snowflake,
                member_snowflake=administrator.member_snowflake,
                role_snowflake=role.id,
            )
            revoked_members.setdefault(role.guild.id, {}).setdefault(
                role.id, []
            ).append(member)
        members = revoked_members.get(role.guild.id, {}).get(role.id, [])
    else:
        action = "granted"
        granted_members: dict[int, dict[int, list[discord.Member]]] = {}
        granted_members.setdefault(role.guild.id, {})[role.id] = []
        administrator_role = AdministratorRole(
            guild_snowflake=role.guild.id, role_snowflake=role.id
        )
        await database_factory.create(administrator_role)
        for member in role.members:
            await administrator_service.added_role(
                guild_snowflake=role.guild.id,
                member_snowflake=member.id,
                role_snowflake=role.id,
            )
            granted_members[role.guild.id][role.id].append(member)
        members = granted_members.get(role.guild.id, {}).get(role.id, [])
    embed = discord.Embed(
        title=title,
        description=f"`{role.name}` was `{action}`.",
        color=discord.Color.red() if action == "revoked" else discord.Color.green(),
    )
    embed.add_field(name="Role ID", value=str(role.id), inline=False)
    embed.add_field(name="Guild", value=str(role.guild.name), inline=False)

    chunks = []
    chunk = []
    pages: list[discord.Embed] = []
    for member in members:
        chunk.append(member)
        if len(chunk) == list_service.CHUNK_SIZE:
            chunks.append(chunk)
            chunk = []
    if chunk:
        chunks.append(chunk)
    field_count = 1
    page_number = 1
    for chunk in chunks:
        embed = discord.Embed(
            title=f"Members {action.capitalize()}",
            color=(
                discord.Color.red() if action == "revoked" else discord.Color.green()
            ),
        )
        for member in chunk:
            embed.add_field(
                name=f"{field_count}. {member}",
                value=f"{member.mention} ({member.id})",
                inline=False,
            )
            field_count += 1
        embed.set_footer(text=f"Page {page_number}")
        pages.append(embed)
        page_number += 1
    return pages


async def get_administrator_roles(guild_snowflake=None) -> list[int]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    roles = await database_factory.select(
        guild_snowflake=guild_snowflake, singular=False
    )
    role_snowflakes = []
    for role in roles:
        role_snowflakes.append(role.role_snowflake)
    return role_snowflakes
