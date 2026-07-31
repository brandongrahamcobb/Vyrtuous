"""!/bin/python3
autoassign_role_service.py The purpose of this program is to service autoassign roles.

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
from vyrtuous.cache.permissions import PermissionGroup
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.db.autoassign import AutoAssignRole
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.permission_entry import PermissionEntry
from vyrtuous.listing import list_service
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import GuildNotFound, MemberNotFound, RoleNotFound
from vyrtuous.utils.messaging import emojis

MODEL = AutoAssignRole


async def toggle_autoassign_role(
    author_snowflake: int,
    group: PermissionGroup,
    guild_snowflake: int,
    message_snowflake: int,
    message_channel_snowflake: int,
    role_snowflake: int,
) -> list[discord.Embed]:
    title = f"{emojis.get_random_emoji()} Autoassign Roles"
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    role = guild.get_role(role_snowflake)
    if role is None:
        raise RoleNotFound(str(role_snowflake))
    autoassign_database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    autoassign_roles = await autoassign_database_factory.select(
        group_alias=group.alias,
        guild_snowflake=guild_snowflake,
        role_snowflake=role_snowflake,
        singular=True,
    )
    group_database_factory: DatabaseFactory = DatabaseFactory(PermissionEntry)
    group_members = await group_database_factory.select(
        group_alias=group.alias,
        guild_snowflake=role.guild.id,
        role_snowflake=role_snowflake,
        singular=False,
    )
    if autoassign_roles:
        action = "revoked"
        revoked_members: dict[int, dict[int, list[str]]] = {}
        for group_member in group_members:
            member = role.guild.get_member(group_member.member_snowflake)
            if member is None:
                simplified_member = bot.registry.get(MemberState).active.get(
                    group_member.member_snowflake, None
                )
                if simplified_member is None:
                    raise MemberNotFound(str(group_member.member_snowflake))
                else:
                    display_name = simplified_member[0]
            else:
                display_name = member.mention
                await removed_role(
                    guild_snowflake=role.guild.id,
                    member_snowflake=group_member.member_snowflake,
                    role_snowflake=role_snowflake,
                )
            await permission_service.log_xgroup(
                author_snowflake=author_snowflake,
                display=True,
                group_alias=group.alias,
                guild_snowflake=role.guild.id,
                member_snowflake=group_member.member_snowflake,
                message_snowflake=message_snowflake,
                message_channel_snowflake=message_channel_snowflake,
                role_snowflake=role_snowflake,
            )
            revoked_members.setdefault(guild_snowflake, {}).setdefault(
                role_snowflake, []
            ).append(display_name)
        members = revoked_members.get(guild_snowflake, {}).get(role_snowflake, [])
    else:
        action = "granted"
        granted_members: dict[int, dict[int, list[discord.Member]]] = {}
        granted_members.setdefault(role.guild.id, {})[role_snowflake] = []
        autoassign_role = AutoAssignRole(
            group_alias=group.alias,
            guild_snowflake=role.guild.id,
            role_snowflake=role_snowflake,
        )
        await autoassign_database_factory.create(autoassign_role)
        for member in role.members:
            await added_role(
                guild_snowflake=role.guild.id,
                member_snowflake=member.id,
                role_snowflake=role_snowflake,
            )
            await permission_service.log_group(
                author_snowflake=author_snowflake,
                display=True,
                group_alias=group.alias,
                guild_snowflake=role.guild.id,
                member_snowflake=member.id,
                message_snowflake=message_snowflake,
                message_channel_snowflake=message_channel_snowflake,
                role_snowflake=role_snowflake,
            )
            granted_members[role.guild.id][role_snowflake].append(member)
        members = granted_members.get(role.guild.id, {}).get(role_snowflake, [])
    embed = discord.Embed(
        title=title,
        description=f"`{role.name}` was `{action}` `{group.name}.",
        color=discord.Color.red() if action == "revoked" else discord.Color.green(),
    )
    embed.add_field(name="Role ID", value=str(role_snowflake), inline=False)
    embed.add_field(name="Guild", value=str(guild.name), inline=False)
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


async def get_autoassign_roles(guild_snowflake=None) -> list[int]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    roles = await database_factory.select(
        guild_snowflake=guild_snowflake, singular=False
    )
    role_snowflakes = []
    for role in roles:
        role_snowflakes.append(role.role_snowflake)
    return role_snowflakes


async def added_role(
    guild_snowflake: int, member_snowflake: int, role_snowflake: int
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    autoassign_database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    autoassign_role_snowflakes = []
    autoassign_role = await autoassign_database_factory.select(
        guild_snwoflake=guild_snowflake,
        role_snowflake=role_snowflake,
        singular=True,
    )
    group_database_factory: DatabaseFactory = DatabaseFactory(PermissionEntry)
    group_member = await group_database_factory.select(
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
        role_snowflakes=[role_snowflake],
        inside_fields=["role_snowflakes"],
        singular=True,
    )
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        return
    if not group_member:
        permission_entry = PermissionEntry(
            guild_snowflake=int(guild_snowflake),
            group_alias=group_member.group_alias,
            member_snowflake=int(member_snowflake),
            role_snowflakes=[int(role_snowflake)],
        )
        await group_database_factory.create(permission_entry)
        return
    permission_state = bot.registry.get(PermissionState)
    group = permission_state.groups.get(autoassign_role.group_alias, None)
    if group is None:
        name = "Unknown"
    else:
        name = group.name
    autoassign_role_snowflakes = group_member.role_snowflakes
    autoassign_role_snowflakes.append(role_snowflake)
    where_kwargs = {
        "guild_snowflake": int(guild_snowflake),
        "member_snowflake": member_snowflake,
    }
    set_kwargs = {"role_snowflakes": autoassign_role_snowflakes}
    await group_database_factory.update(
        set_kwargs=set_kwargs, where_kwargs=where_kwargs
    )
    bot.logger.debug(
        f"Granted permission group ({name}) to member ({member_snowflake}) in guild ({guild_snowflake})."
    )


async def removed_role(
    guild_snowflake: int, member_snowflake: int, role_snowflake: int
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    autoassign_database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    autoassign_role = await autoassign_database_factory.select(
        guild_snwoflake=guild_snowflake,
        role_snowflake=role_snowflake,
        singular=True,
    )
    group_database_factory: DatabaseFactory = DatabaseFactory(PermissionEntry)
    group_member = await group_database_factory.select(
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
        role_snowflakes=[role_snowflake],
        inside_fields=["role_snowflakes"],
        singular=True,
    )
    if not group_member:
        return
    permission_state = bot.registry.get(PermissionState)
    group = permission_state.groups.get(autoassign_role.group_alias, None)
    if group is None:
        name = "Unknown"
    else:
        name = group.name
    autoassign_role_snowflakes = group_member.role_snowflakes
    autoassign_role_snowflakes.remove(role_snowflake)
    if not autoassign_role_snowflakes:
        await group_database_factory.delete_by_cls(group_member)
    else:
        where_kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "member_snowflake": member_snowflake,
        }
        set_kwargs = {"role_snowflakes": autoassign_role_snowflakes}
        await group_database_factory.update(
            set_kwargs=set_kwargs, where_kwargs=where_kwargs
        )
        bot.logger.debug(
            f"Revoked permission group ({name}) to member ({member_snowflake}) in guild ({guild_snowflake})."
        )
