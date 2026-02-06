"""administrator.py The purpose of this program is to inherit from the DatabaseFactory to provide the administrator role.

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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.author import resolve_author
from vyrtuous.commands.errors import NotAdministrator
from vyrtuous.db.roles.admin.administrator import Administrator, AdministratorRole
from vyrtuous.db.roles.dev.developer_service import is_developer_wrapper
from vyrtuous.db.roles.owner.guild_owner_service import is_guild_owner_wrapper
from vyrtuous.db.roles.sysadmin.sysadmin_service import is_sysadmin_wrapper
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    generate_skipped_roles,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger


async def is_administrator_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    return await is_administrator(
        guild_snowflake=int(source.guild.id),
        member_snowflake=int(member_snowflake),
    )


async def is_administrator(guild_snowflake: int, member_snowflake: int) -> bool:
    administrator = await Administrator.select(
        guild_snowflake=int(guild_snowflake),
        member_snowflake=int(member_snowflake),
        singular=True,
    )
    if not administrator:
        raise NotAdministrator
    return True


async def is_administrator_at_all(
    member_snowflake: int,
):
    administrator = await Administrator.select(
        member_snowflake=int(member_snowflake),
    )
    if not administrator:
        raise NotAdministrator
    return True


def administrator_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (
            is_sysadmin_wrapper,
            is_developer_wrapper,
            is_guild_owner_wrapper,
            is_administrator_wrapper,
        ):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise NotAdministrator

    predicate._permission_level = "Administrator"
    return commands.check(predicate)


class AdministratorService:

    lines, pages = [], []

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        administrators = await Administrator.select(singular=False, **where_kwargs)
        for administrator in administrators:
            dictionary.setdefault(administrator.guild_snowflake, {"members": {}})
            dictionary[administrator.guild_snowflake]["members"].setdefault(
                administrator.member_snowflake, {"administrators": {}}
            )

            for role_snowflake in administrator.role_snowflakes:
                dictionary[administrator.guild_snowflake]["members"][
                    administrator.member_snowflake
                ]["administrators"].update({role_snowflake: True})
        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_members = generate_skipped_members(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        if is_at_home:
            if skipped_guilds:
                AdministratorService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                AdministratorService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_members,
                        title="Skipped Members in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Administrators {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await AdministratorService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            admin_n = 0
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, administrators_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                primary_dictionary = administrators_dictionary.get("administrators", {})
                role_mentions = [
                    guild.get_role(role_snowflake).mention
                    for role_snowflake in primary_dictionary
                    if guild.get_role(role_snowflake)
                ]
                if not thumbnail and isinstance(
                    object_dict.get("object", None), discord.Member
                ):
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                else:
                    member_line = f"**User:** {member.display_name} {member.mention}"
                    if role_mentions:
                        member_line += "\n**Roles:** " + "\n".join(role_mentions)
                admin_n += 1
                AdministratorService.lines.append(member_line)
                field_count += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(AdministratorService.lines),
                        inline=False,
                    )
                    embed = flush_page(
                        embed, AdministratorService.pages, title, guild.name
                    )
                    AdministratorService.lines = []
                    field_count = 0
            if AdministratorService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(AdministratorService.lines),
                    inline=False,
                )
            AdministratorService.pages.append(embed)
            AdministratorService.pages[0].description = f'{guild.name} **({admin_n})**'
        return AdministratorService.pages


class AdministratorRoleService:

    lines, pages = [], []

    @classmethod
    async def added_role(cls, updated_kwargs):
        administrator_role_snowflakes = []
        guild_snowflake = updated_kwargs.get("guild_snowflake", None)
        member_snowflake = updated_kwargs.get("member_snowflake", None)
        role_snowflake = updated_kwargs.get("role_snowflake", None)
        where_kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "role_snowflake": int(role_snowflake),
        }
        administrator_role = await AdministratorRole.select(
            **where_kwargs, singular=True
        )
        if not administrator_role:
            return
        administrator = await Administrator.select(
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
            singular=True,
        )
        if not administrator:
            administrator = Administrator(
                guild_snowflake=int(guild_snowflake),
                member_snowflake=int(member_snowflake),
                role_snowflakes=[int(role_snowflake)],
            )
            await administrator.create()
        else:
            administrator_role_snowflakes = administrator.role_snowflakes
            administrator_role_snowflakes.append(role_snowflake)
            where_kwargs = {
                "guild_snowflake": int(guild_snowflake),
                "member_snowflake": member_snowflake,
            }
            set_kwargs = {"role_snowflakes": administrator_role_snowflakes}
            await Administrator.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
            logger.info(
                f"Granted administrator to member ({member_snowflake}) in guild ({guild_snowflake})."
            )

    @classmethod
    async def removed_role(cls, updated_kwargs):
        administrator_role_snowflakes = []
        guild_snowflake = updated_kwargs.get("guild_snowflake", None)
        member_snowflake = updated_kwargs.get("member_snowflake", None)
        role_snowflake = updated_kwargs.get("role_snowflake", None)
        where_kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "role_snowflake": int(role_snowflake),
        }
        administrator_role = await AdministratorRole.select(
            **where_kwargs, singular=True
        )
        if not administrator_role:
            return
        administrator = await Administrator.select(
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
            role_snowflakes=[role_snowflake],
            singular=True,
            inside_fields=["role_snowflakes"],
        )
        if administrator:
            administrator_role_snowflakes = administrator.role_snowflakes
            administrator_role_snowflakes.remove(role_snowflake)
            where_kwargs = {
                "guild_snowflake": int(guild_snowflake),
                "member_snowflake": member_snowflake,
            }
            set_kwargs = {"role_snowflakes": administrator_role_snowflakes}
            await Administrator.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
            logger.info(
                f"Revoked administrator from member ({member_snowflake}) in guild ({guild_snowflake})."
            )

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        administrator_roles = await AdministratorRole.select(
            singular=False, **where_kwargs
        )
        for administrator_role in administrator_roles:
            dictionary.setdefault(administrator_role.guild_snowflake, {"roles": {}})
            dictionary[administrator_role.guild_snowflake]["roles"].setdefault(
                administrator_role.role_snowflake, {}
            )

        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_roles = generate_skipped_roles(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_roles=skipped_roles,
        )
        if is_at_home:
            if skipped_guilds:
                AdministratorRoleService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_roles:
                AdministratorRoleService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_roles,
                        title="Skipped Roles in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        dictionary = {}
        title = f"{get_random_emoji()} Administrator Roles"
        where_kwargs = object_dict.get("columns", None)

        dictionary = await AdministratorRoleService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            admin_role_n = 0
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for role_snowflake, entry in guild_data.get("roles", {}).items():
                role = guild.get_role(role_snowflake)
                if not role:
                    continue
                if field_count >= CHUNK_SIZE:
                    embed = flush_page(
                        embed, AdministratorRoleService.pages, title, guild.name
                    )
                embed.add_field(name=role.name, value=role.mention, inline=False)
                field_count += 1
                admin_role_n += 1
            AdministratorRoleService.pages.append(embed)
            AdministratorRoleService.pages[0].description = f'{guild.name} **({admin_role_n})**'
        return AdministratorRoleService.pages

    @classmethod
    async def toggle_administrator_role(cls, updated_kwargs, role_dict):
        bot = DiscordBot.get_instance()
        guild_snowflake = updated_kwargs.get("guild_snowflake", None)
        chunk_size, field_count, pages = 7, 0, []
        guild = bot.get_guild(guild_snowflake)
        title = f"{get_random_emoji()} Administrators and Roles"
        kwargs = role_dict.get("columns", None)

        administrators = await Administrator.select(
            role_snowflakes=role_dict.get("id", None), inside_fields=["role_snowflakes"]
        )
        administrator_roles = await AdministratorRole.select(
            role_snowflake=role_dict.get("id", None)
        )

        if administrator_roles:
            for administrator_role in administrator_roles:
                await AdministratorRole.delete(**kwargs)
            action = "revoked"
        else:
            action = "granted"
        if administrators:
            for administrator in administrators:
                revoked_members = {}
                member = guild.get_member(administrator.member_snowflake)
                await Administrator.delete(
                    guild_snowflake=int(guild_snowflake),
                    member_snowflake=administrator.member_snowflake,
                )
                revoked_members.setdefault(
                    administrator.guild_snowflake, {}
                ).setdefault(role_dict.get("id", None), []).append(member)
        else:
            revoked_members = {}

        if action != "revoked":
            granted_members = {}
            granted_members.setdefault(guild_snowflake, {})[
                role_dict.get("id", None)
            ] = []
            administrator_role = AdministratorRole(**kwargs)
            await administrator_role.create()
            role_snowflakes = [role_dict.get("id", None)]
            for member in role_dict.get("object", None).members:
                administrator = Administrator(
                    guild_snowflake=int(guild_snowflake),
                    member_snowflake=member.id,
                    role_snowflakes=role_snowflakes,
                )
                await administrator.create()
                granted_members[guild_snowflake][role_dict.get("id", None)].append(
                    member
                )

        embed = discord.Embed(
            title=title,
            description=f"`{role_dict.get('name', None)}` was `{action}`.",
            color=discord.Color.red() if action == "revoked" else discord.Color.green(),
        )
        embed.add_field(
            name="Role ID", value=str(role_dict.get("id", None)), inline=False
        )
        embed.add_field(name="Guild", value=guild.name, inline=False)
        pages.append(embed)
        if action == "revoked":
            members = revoked_members.get(guild_snowflake, {}).get(
                role_dict.get("id", None), []
            )
        else:
            members = granted_members.get(guild_snowflake, {}).get(
                role_dict.get("id", None), []
            )

        chunks = []
        chunk = []
        for member in members:
            chunk.append(member)
            if len(chunk) == chunk_size:
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
                    discord.Color.red()
                    if action == "revoked"
                    else discord.Color.green()
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
