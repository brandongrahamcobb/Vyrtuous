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
from vyrtuous.models.duration import DurationObject
from vyrtuous.utils.errors.error import GuildNotFound, RoleNotFound
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.tracking import stream_service

MODEL = AutoAssignRole


async def toggle_autoassign_role(
    author_snowflake: int,
    group: PermissionGroup,
    guild_snowflake: int,
    role_snowflake: int,
) -> list[discord.Embed]:
    is_channel_scope = False
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    role = guild.get_role(role_snowflake)
    if role is None:
        raise RoleNotFound(str(role_snowflake))
    autoassign_database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    autoassign_role = await autoassign_database_factory.select(
        group_alias=group.alias,
        guild_snowflake=guild_snowflake,
        role_snowflake=role_snowflake,
        singular=True,
    )
    group_database_factory: DatabaseFactory = DatabaseFactory(PermissionEntry)
    group_members = await group_database_factory.select(
        group_alias=group.alias,
        guild_snowflake=role.guild.id,
        role_snowflakes=role_snowflake,
        inside_fields=["role_snowflakes"],
        singular=False,
    )
    if autoassign_role:
        members = []
        for group_member in group_members:
            member = guild.get_member(group_member.member_snowflake)
            if member is None:
                simplified_member = bot.registry.get(MemberState).active.get(
                    group_member.member_snowflake, None
                )
                if simplified_member is None:
                    continue
                else:
                    member_snowflake = group_member.member_snowflake
            else:
                if member.voice and member.voice.channel:
                    is_channel_scope = True
                member_snowflake = member.id
            await removed_role(
                guild_snowflake=role.guild.id,
                member_snowflake=member_snowflake,
                role_snowflake=role_snowflake,
            )
            await stream_service.log(
                author_snowflake=author_snowflake,
                channel_snowflake=None,
                display=True,
                duration=DurationObject(number=0, prefix="", sign=1, unit=""),
                identifier=f"x{group.alias}",
                is_channel_scope=is_channel_scope,
                guild_snowflake=role.guild.id,
                member_snowflake=group_member.member_snowflake,
                message_snowflake=None,
                message_channel_snowflake=None,
                reason="No reason provided.",
                role_snowflake=role_snowflake,
                target=None,
            )
            members.append(member)
        pages = await disable(
            group=group,
            guild_snowflake=guild_snowflake,
            members=members,
            role_snowflake=role_snowflake,
        )
        await autoassign_database_factory.delete_by_cls(autoassign_role)
    else:
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
            await stream_service.log(
                author_snowflake=author_snowflake,
                channel_snowflake=None,
                display=True,
                duration=DurationObject(number=0, prefix="", sign=1, unit=""),
                identifier=group.alias,
                is_channel_scope=is_channel_scope,
                guild_snowflake=role.guild.id,
                member_snowflake=member.id,
                message_snowflake=None,
                message_channel_snowflake=None,
                reason="No reason provided.",
                role_snowflake=role_snowflake,
                target=None,
            )
        pages = await enable(
            group=group,
            guild_snowflake=guild_snowflake,
            members=role.members,
            role_snowflake=role_snowflake,
        )
    return pages


async def enable(
    group: PermissionGroup,
    guild_snowflake: int,
    members: list[discord.Member],
    role_snowflake: int,
) -> list[discord.Embed]:
    pages: list[discord.Embed] = []
    chunks = []
    chunk = []
    action = "granted"
    title = f"{emojis.get_random_emoji()} Autoassign Roles"
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    role = guild.get_role(role_snowflake)
    if role is None:
        raise RoleNotFound(str(role_snowflake))
    embed = discord.Embed(
        title=title,
        description=f"`{role.name}` was {action} `{group.name}`.",
        color=discord.Color.red() if action == "revoked" else discord.Color.green(),
    )
    embed.add_field(name="Role ID", value=str(role_snowflake), inline=False)
    embed.add_field(name="Guild", value=str(guild.name), inline=False)
    pages.append(embed)
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
            title=f"Members",
            color=discord.Color.green(),
        )
        for member in chunk:
            embed.add_field(
                name=f"{field_count}",
                value=f"{member}",
                inline=True,
            )
            field_count += 1
        embed.set_footer(text=f"Page {page_number}")
        pages.append(embed)
        page_number += 1
    return pages


async def disable(
    group: PermissionGroup,
    guild_snowflake: int,
    members: list[discord.Member],
    role_snowflake: int,
) -> list[discord.Embed]:
    pages: list[discord.Embed] = []
    chunks = []
    chunk = []
    action = "revoked"
    title = f"{emojis.get_random_emoji()} Autoassign Roles"
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    role = guild.get_role(role_snowflake)
    if role is None:
        raise RoleNotFound(str(role_snowflake))
    embed = discord.Embed(
        title=title,
        description=f"`{role.name}` was {action} `{group.name}`.",
        color=discord.Color.red() if action == "revoked" else discord.Color.green(),
    )
    embed.add_field(name="Role ID", value=str(role_snowflake), inline=False)
    embed.add_field(name="Guild", value=str(guild.name), inline=False)
    pages.append(embed)
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
            title=f"Members",
            color=discord.Color.red(),
        )
        for member in chunk:
            embed.add_field(
                name=f"{field_count}",
                value=f"{member}",
                inline=True,
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
        role_snowflakes=role_snowflake,
        inside_fields=["role_snowflakes"],
        singular=True,
    )
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        return
    if not group_member:
        permission_entry = PermissionEntry(
            guild_snowflake=int(guild_snowflake),
            group_alias=autoassign_role.group_alias,
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
        role_snowflakes=role_snowflake,
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
