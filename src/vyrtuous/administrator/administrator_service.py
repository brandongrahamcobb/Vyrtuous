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

from copy import copy
from dataclasses import dataclass, field
from typing import Dict, List

import discord
from discord.ext import commands

from vyrtuous.administrator.administrator import Administrator, AdministratorRole


class NotAdministrator(commands.CheckFailure):
    def __init__(
        self,
        message="Member is not an administrator in this server.",
    ):
        super().__init__(message)


@dataclass(frozen=True)
class AdministratorDictionary:
    data: Dict[int, Dict[int, Dict[int, bool]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


class AdministratorService:
    __CHUNK_SIZE = 12
    MODEL = Administrator

    def __init__(
        self,
        *,
        author_service=None,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        emoji=None,
    ):
        self.__author_service = author_service
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji

    async def is_administrator_wrapper(
        self,
        context,
    ):
        return await self.is_administrator(
            guild_snowflake=int(context.guild.id),
            member_snowflake=int(context.author.id),
        )

    async def is_administrator(
        self, guild_snowflake: int, member_snowflake: int
    ) -> bool:
        administrator = await self.__database_factory.select(
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
            singular=True,
        )
        if not administrator:
            raise NotAdministrator
        return True

    async def is_administrator_at_all(
        self,
        member_snowflake: int,
    ):
        administrator = await self.__database_factory.select(
            member_snowflake=int(member_snowflake),
        )
        if not administrator:
            raise NotAdministrator
        return True

    async def build_dictionary(self, kwargs):
        dictionary = {}
        administrators = await self.__database_factory.select(singular=False, **kwargs)
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

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []

        obj = object_dict.get("object")
        obj_name = "All Servers"
        if isinstance(obj, discord.Guild):
            obj_name = obj.name
        elif isinstance(obj, discord.abc.GuildChannel):
            obj_name = obj.name
        elif isinstance(obj, discord.Member):
            obj_name = object_dict.get("name", None)
        title = f"{self.__emoji.get_random_emoji()} Administrators for {obj_name}"

        dictionary = await self.build_dictionary(
            kwargs=object_dict.get("columns", None)
        )
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=AdministratorDictionary,
            dictionary=dictionary,
        )

        for guild_snowflake, guild_data in processed_dictionary.data.items():
            admin_n = 0
            field_count = 0
            lines = []
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, processed_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                role_mentions = [
                    guild.get_role(role_snowflake).mention
                    for role_snowflake in processed_dictionary.get("administrators", {})
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
                    lines.append(member_line)
                admin_n += 1
                field_count += 1
                if field_count >= self.__CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = self.__dictionary_service.flush_page(
                        embed, pages, title, guild.name
                    )
                    lines = []
                    field_count = 0
            if lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
            original_description = embed.description or ""
            embed.description = f"**{original_description}** **({admin_n})**"
            pages.append(embed)
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_members)
        return pages

    async def administrators_by_role(self, role_snowflake: int):
        administrators = await self.__database_factory.select(
            role_snowflakes=int(role_snowflake), inside_fields=["role_snowflakes"]
        )
        if administrators:
            return administrators
        return []

    async def administrator_existing(self, kwargs: dict[str, int]):
        guild_snowflake = kwargs.get("guild_snowflake", None)
        member_snowflake = kwargs.get("member_snowflake", None)
        role_snowflake = kwargs.get("role_snowflake", None)
        administrator = await self.__database_factory.select(
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
            role_snowflakes=int(role_snowflake),
            inside_fields=["role_snowflakes"],
            singular=True,
        )
        if administrator:
            return administrator
        return None

    async def added_role(self, kwargs):
        administrator_role_snowflakes = []
        guild_snowflake = kwargs.get("guild_snowflake", None)
        member_snowflake = kwargs.get("member_snowflake", None)
        role_snowflake = kwargs.get("role_snowflake", None)
        administrator_existing = await self.administrator_existing(kwargs=kwargs)
        if not administrator_existing:
            administrator = Administrator(
                guild_snowflake=int(guild_snowflake),
                member_snowflake=int(member_snowflake),
                role_snowflakes=[int(role_snowflake)],
            )
            await self.__database_factory.create(administrator)
            return
        administrator_role_snowflakes = administrator_existing.role_snowflakes
        administrator_role_snowflakes.append(role_snowflake)
        where_kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "member_snowflake": member_snowflake,
        }
        set_kwargs = {"role_snowflakes": administrator_role_snowflakes}
        await self.__database_factory.update(
            set_kwargs=set_kwargs, where_kwargs=where_kwargs
        )
        self.__bot.logger.info(
            f"Granted administrator to member ({member_snowflake}) in guild ({guild_snowflake})."
        )

    async def removed_role(self, kwargs):
        administrator_role_snowflakes = []
        guild_snowflake = kwargs.get("guild_snowflake", None)
        member_snowflake = kwargs.get("member_snowflake", None)
        role_snowflake = kwargs.get("role_snowflake", None)
        where_kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "role_snowflake": int(role_snowflake),
        }
        administrator = await self.__database_factory.select(
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
            role_snowflakes=role_snowflake,
            singular=True,
            inside_fields=["role_snowflakes"],
        )
        if not administrator:
            return
        administrator_role_snowflakes = administrator.role_snowflakes
        administrator_role_snowflakes.remove(role_snowflake)
        if administrator_role_snowflakes == []:
            await self.__database_factory.delete(
                guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
            )
        else:
            where_kwargs = {
                "guild_snowflake": int(guild_snowflake),
                "member_snowflake": member_snowflake,
            }
            set_kwargs = {"role_snowflakes": administrator_role_snowflakes}
            await self.__database_factory.update(
                set_kwargs=set_kwargs, where_kwargs=where_kwargs
            )
        self.__bot.logger.info(
            f"Revoked administrator from member ({member_snowflake}) in guild ({guild_snowflake})."
        )


@dataclass(frozen=True)
class AdministratorRoleDictionary:
    data: Dict[int, Dict[str, Dict[int, dict]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_roles: List[discord.Embed] = field(default_factory=list)


class AdministratorRoleService:
    __CHUNK_SIZE = 12
    MODEL = AdministratorRole

    def __init__(
        self,
        *,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        emoji=None,
    ):
        self.__administrator_service = AdministratorService(
            bot=bot, database_factory=database_factory
        )
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji

    async def is_added_role_administrator(self, guild_snowflake, role_snowflake):
        where_kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "role_snowflake": int(role_snowflake),
        }
        administrator_roles = await self.__database_factory.select(
            **where_kwargs,
        )
        if not administrator_roles:
            return False
        return True

    async def build_full_dictionary(self, where_kwargs):
        dictionary = {}
        administrator_roles = await self.__database_factory.select(
            singular=False, **where_kwargs
        )
        for administrator_role in administrator_roles:
            dictionary.setdefault(administrator_role.guild_snowflake, {"roles": {}})
            dictionary[administrator_role.guild_snowflake]["roles"].setdefault(
                administrator_role.role_snowflake, {}
            )
        return dictionary

    async def build_pages(self, object_dict, is_at_home):
        pages = []
        title = f"{self.__emoji.get_random_emoji()} Administrator Roles"

        where_kwargs = object_dict.get("columns", None)

        full_dictionary = await self.build_full_dictionary(where_kwargs=where_kwargs)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=AdministratorRoleDictionary, dictionary=full_dictionary
        )

        admin_role_n = 0
        for guild_snowflake, guild_data in processed_dictionary.data.items():
            field_count = 0
            lines = []
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for role_snowflake, entry in guild_data.get("roles", {}).items():
                role = guild.get_role(role_snowflake)
                if not role:
                    continue
                if field_count >= self.__CHUNK_SIZE:
                    embed = self.__dictionary_service.flush_page(
                        embed, pages, title, guild.name
                    )
                embed.add_field(name=role.name, value=role.mention, inline=False)
                field_count += 1
                admin_role_n += 1
            pages.append(embed)
        if pages:
            original_description = embed.description or ""
            embed.description = f"**{original_description}** **({admin_role_n})**"
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_roles)
        return pages

    async def toggle_administrator_role(self, role_dict):
        guild_snowflake = role_dict.get("columns", None).get("guild_snowflake")
        guild = self.__bot.get_guild(guild_snowflake)
        title = f"{self.__emoji.get_random_emoji()} Administrators and Roles"
        role_id = role_dict.get("id")
        administrators = await self.__administrator_service.administrators_by_role(
            role_snowflake=role_id
        )
        administrator_roles = await self.is_added_role_administrator(
            guild_snowflake=guild_snowflake, role_snowflake=role_id
        )
        if administrator_roles:
            action = "revoked"
            if administrator_roles:
                await self.__database_factory.delete(**role_dict.get("columns", None))
            revoked_members = {}
            for administrator in administrators:
                member = guild.get_member(administrator.member_snowflake)
                await self.__administrator_service.removed_role(
                    {
                        "guild_snowflake": guild_snowflake,
                        "member_snowflake": administrator.member_snowflake,
                        "role_snowflake": role_id,
                    }
                )
                revoked_members.setdefault(guild_snowflake, {}).setdefault(
                    role_id, []
                ).append(member)
            members = revoked_members.get(guild_snowflake, {}).get(role_id, [])
        else:
            action = "granted"
            granted_members = {}
            granted_members.setdefault(guild_snowflake, {})[role_id] = []
            administrator_role = AdministratorRole(**role_dict.get("columns", None))
            await self.__database_factory.create(administrator_role)
            for member in role_dict.get("object").members:
                await self.__administrator_service.added_role(
                    {
                        "guild_snowflake": guild_snowflake,
                        "member_snowflake": member.id,
                        "role_snowflake": role_id,
                    }
                )
                granted_members[guild_snowflake][role_id].append(member)
            members = granted_members.get(guild_snowflake, {}).get(role_id, [])

        embed = discord.Embed(
            title=title,
            description=f"`{role_dict.get('name')}` was `{action}`.",
            color=discord.Color.red() if action == "revoked" else discord.Color.green(),
        )
        embed.add_field(name="Role ID", value=str(role_id), inline=False)
        embed.add_field(name="Guild", value=guild.name, inline=False)

        chunks = []
        chunk = []
        pages = []
        for member in members:
            chunk.append(member)
            if len(chunk) == self.__CHUNK_SIZE:
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
