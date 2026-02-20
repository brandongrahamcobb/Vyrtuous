"""!/bin/python3
role_service.py The purpose of this program is to extend AliasService to service the role class.

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

from vyrtuous.role.role import Role


@dataclass
class RoleDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, int]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


class RoleService:
    __CHUNK_SIZE = 7
    MODEL = Role

    def __init__(
        self,
        *,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        emoji=None,
        stream_service=None,
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji
        self.__stream_service = stream_service

    async def act_embed(self, ctx):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        channel = guild.get_channel(ctx.target_channel_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        role = guild.get_role(ctx.target_role_snowflake)
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} {member.display_name} has been granted a role",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def undo_embed(self, ctx):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        channel = guild.get_channel(ctx.target_channel_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        role = guild.get_role(ctx.target_role_snowflake)
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} {member.display_name}'s role has been revoked",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def administer_role(self, guild_snowflake, member_snowflake, role_snowflake):
        guild = self.__bot.get_guild(guild_snowflake)
        member = guild.get_member(member_snowflake)
        role = guild.get_role(role_snowflake)
        try:
            await member.add_roles(role, reason="Granting role.")
        except discord.Forbidden as e:
            self.__bot.logger.error(str(e).capitalize())

    async def revoke_role(self, guild_snowflake, member_snowflake, role_snowflake):
        guild = self.__bot.get_guild(guild_snowflake)
        member = guild.get_member(member_snowflake)
        role = guild.get_role(role_snowflake)
        try:
            await member.remove_roles(role, reason="Revoking role.")
        except discord.Forbidden as e:
            self.__bot.logger.error(str(e).capitalize())

    async def added_role(
        self,
        category_class,
        category_role_class,
        guild_snowflake,
        member_snowflake,
        role_snowflake,
    ):
        kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "role_snowflake": str(role_snowflake),
        }
        role = await category_role_class.select(singular=True, **kwargs)
        if role:
            kwargs.update({"channel_snowflake": role.channel_snowflake})
            kwargs.update({"member_snowflake": role.member_snowflake})
            msg = f"Member ({member_snowflake}) was granted the role ({role_snowflake}) for category ({category_class.__name__()}) related to channel ({role.channel_snowflake}) in guild ({guild_snowflake})."
            category = category_class(**kwargs)
            await category.create()
            self.__bot.logger.info(msg)
        else:
            return

    async def removed_role(
        self,
        category_class,
        category_role_class,
        guild_snowflake,
        member_snowflake,
        role_snowflake,
    ):
        kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "role_snowflake": str(role_snowflake),
        }
        role = await category_role_class.select(singular=True, **kwargs)
        if role:
            kwargs.update({"channel_snowflake": role.channel_snowflake})
            kwargs.update({"member_snowflake": role.member_snowflake})
            msg = f"Member ({member_snowflake}) was revoked the role ({role_snowflake}) for category ({category_class.__name__()}) related to channel ({role.channel_snowflake}) in guild ({guild_snowflake})."
            await category_class.delete(**kwargs)
            self.__bot.logger.info(msg)
        else:
            return

    async def build_dictionary(self, where_kwargs):
        pages = []
        dictionary = {}
        roles = await Role.select(singular=False, **where_kwargs)
        for role in roles:
            dictionary.setdefault(role.guild_snowflake, {"members": {}})
            dictionary[role.guild_snowflake]["members"].setdefault(
                role.member_snowflake, {"roles": {}}
            )
            dictionary[role.guild_snowflake]["members"][role.member_snowflake][
                "roles"
            ].setdefault(role.channel_snowflake, {})
            dictionary[role.guild_snowflake]["members"][role.member_snowflake]["roles"][
                role.channel_snowflake
            ] = {"role": role.role_snowflake}
        return dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = []
        title = f"{self.__emoji.get_random_emoji()} Role {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_dictionary(where_kwargs=where_kwargs)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=RoleDictionary, dictionary=dictionary
        )
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_members)

        role_n = 0
        for guild_snowflake, guild_data in processed_dictionary.data.items():
            field_count = 0
            thumbnail = False
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, role_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in role_dictionary.get(
                    "roles"
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    role = channel_dictionary.get("role", None)
                    if not isinstance(
                        object_dict.get("object"), discord.abc.GuildChannel
                    ):
                        lines.append(f"**Channel:** {channel.mention}")
                    lines.append(f"**Role:** {role.mention}")
                    role_n += 1
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
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)
        if pages:
            pages[0].description = f"**({role_n})**"
        return pages

    async def enforce(self, ctx, source, state):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        added_role = self.MODEL(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
            role_snowflake=ctx.target_role_snowflake,
        )
        await self.__database_factory.create(added_role)
        role = guild.get_role(ctx.target_role_snowflake)
        if role:
            await self.administer_role(
                guild_snowflake=ctx.source_guild_snowflake,
                member_snowflake=ctx.target_member_snowflake,
                role_snowflake=ctx.target_role_snowflake,
            )
        await self.__stream_service.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            identifier="role",
            member=member,
            source=source,
        )
        embed = await self.act_embed(ctx=ctx)
        return await state.end(success=embed)

    async def undo(self, ctx, source, state):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        await self.__database_factory.delete(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
            role_snowflake=ctx.target_role_snowflake,
        )
        role = guild.get_role(ctx.target_role_snowflake)
        if role:
            await self.revoke_role(
                guild_snowflake=ctx.source_guild_snowflake,
                member_snowflake=ctx.target_member_snowflake,
                role_snowflake=ctx.target_role_snowflake,
            )
        await self.__stream_service.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            identifier="unrole",
            is_modification=True,
            member=member,
            source=source,
        )
        embed = await self.undo_embed(ctx=ctx)
        return await state.end(success=embed)
