"""!/bin/python3
moderator_service.py The purpose of this program is to extend Service to service the moderator class.

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

from vyrtuous.administrator.administrator_service import AdministratorService
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.moderator.moderator import Moderator
from vyrtuous.sysadmin.sysadmin_service import SysadminService


class NotModerator(commands.CommandError):
    def __init__(
        self,
        message="Member is not a moderator in this channel.",
    ):
        super().__init__(message)


class ModeratorService:
    __CHUNK_SIZE = 7
    MODEL = Moderator

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
        self.__database_factory = database_factory
        self.__dictionary_service = dictionary_service
        self.__dictionary_service.model = self.MODEL
        self.__emoji = emoji
        self.__sysadmin_service = SysadminService(
            author_service=author_service, bot=bot
        )
        self.__developer_service = DeveloperService(
            author_service=author_service, bot=bot, database_factory=database_factory
        )
        self.__guild_owner_service = GuildOwnerService(
            author_service=author_service, bot=bot, database_factory=database_factory
        )
        self.__administrator_service = AdministratorService(
            author_service=author_service, bot=bot, database_factory=database_factory
        )
        self.__coordinator_service = CoordinatorService(
            author_service=author_service,
            bot=bot,
            database_factory=database_factory,
            emoji=emoji,
        )

    async def is_moderator_wrapper(
        self,
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        member = self.__author_service.resolve_author(source=source)
        member_snowflake = member.id
        return await self.is_moderator(
            channel_snowflake=source.channel.id,
            guild_snowflake=source.guild.id,
            member_snowflake=int(member_snowflake),
        )

    async def is_moderator(
        self, channel_snowflake: int, guild_snowflake: int, member_snowflake: int
    ) -> bool:
        moderator = await self.__database_factory.select(
            channel_snowflake=int(channel_snowflake),
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
            singular=True,
        )
        if not moderator:
            raise NotModerator
        return True

    async def is_moderator_at_all_wrapper(
        self,
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ) -> bool:
        member = self.__author_service.resolve_author(source=source)
        member_snowflake = member.id
        return await self.is_moderator_at_all(member_snowflake=member_snowflake)

    async def is_moderator_at_all(
        self,
        member_snowflake: int,
    ) -> bool:
        moderator = await self.__database_factory.select(
            member_snowflake=int(member_snowflake)
        )
        if not moderator:
            raise NotModerator
        return True

    def moderator_predicator(self):
        async def predicate(
            source: Union[commands.Context, discord.Interaction, discord.Message],
        ):
            for verify in (
                self.__sysadmin_service.is_sysadmin_wrapper,
                self.__developer_service.is_developer_wrapper,
                self.__guild_owner_service.is_guild_owner_wrapper,
                self.__administrator_service.is_administrator_wrapper,
                self.__coordinator_service.is_coordinator_at_all_wrapper,
                self.is_moderator_at_all_wrapper,
            ):
                try:
                    if await verify(source):
                        return True
                except commands.CheckFailure:
                    continue
            raise NotModerator

        predicate._permission_level = "Moderator"
        return commands.check(predicate)

    async def build_clean_dictionary(self, is_at_home, where_kwargs):
        pages = []
        dictionary = {}
        moderators = await self.__database_factory.select(
            singular=False, **where_kwargs
        )
        for moderator in moderators:
            dictionary.setdefault(moderator.guild_snowflake, {"members": {}})
            dictionary[moderator.guild_snowflake]["members"].setdefault(
                moderator.member_snowflake, {"moderators": {}}
            )
            dictionary[moderator.guild_snowflake]["members"][
                moderator.member_snowflake
            ]["moderators"].setdefault(moderator.channel_snowflake, {})
            dictionary[moderator.guild_snowflake]["members"][
                moderator.member_snowflake
            ]["moderators"][moderator.channel_snowflake].update(
                {"placeholder": "placeholder"}
            )
        skipped_guilds = self.__dictionary_service.generate_skipped_guilds(dictionary)
        skipped_members = self.__dictionary_service.generate_skipped_members(dictionary)
        cleaned_dictionary = self.__dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        if is_at_home:
            if skipped_guilds:
                pages.extend(
                    self.__dictionary_service.generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                pages.extend(
                    self.__dictionary_service.generate_skipped_dict_pages(
                        skipped=skipped_members,
                        title="Skipped Members in Server",
                    )
                )
        return cleaned_dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = []
        title = f"{self.__emoji.get_random_emoji()} Moderator {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        mod_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, member_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                mod_n += 1
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in member_dictionary.get(
                    "moderators", {}
                ).items():
                    if not isinstance(
                        object_dict.get("object", None), discord.abc.GuildChannel
                    ):
                        channel = guild.get_channel(channel_snowflake)
                        if not channel:
                            continue
                        lines.append(f"**Channel:** {channel.mention}")
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
            pages.append(embed)
        if pages:
            pages[0].description = f"**({mod_n})**"
        return pages

    async def survey(self, channel_dict, guild_snowflake):

        chunk_size, pages = 7, []
        (
            sysadmins,
            developers,
            guild_owners,
            administrators,
            coordinators,
            moderators,
        ) = ([], [], [], [], [], [])

        for member in channel_dict.get("object", None).members:
            try:
                if await self.__sysadmin_service.is_sysadmin(member.id):
                    sysadmins.append(member)
            except commands.CheckFailure as e:
                self.__bot.logger.warning(str(e).capitalize())
            try:
                if await self.__developer_service.is_developer(member.id):
                    developers.append(member)
            except commands.CheckFailure as e:
                self.__bot.logger.warning(str(e).capitalize())
            try:
                if await self.__guild_owner_service.is_guild_owner(
                    guild_snowflake, member.id
                ):
                    guild_owners.append(member)
            except commands.CheckFailure as e:
                self.__bot.logger.warning(str(e).capitalize())
            try:
                if await self.__administrator_service.is_administrator(
                    guild_snowflake, member.id
                ):
                    administrators.append(member)
            except commands.CheckFailure as e:
                self.__bot.logger.warning(str(e).capitalize())
            try:
                if await self.__coordinator_service.is_coordinator(
                    channel_dict.get("id", None), guild_snowflake, member.id
                ):
                    coordinators.append(member)
            except commands.CheckFailure as e:
                self.__bot.logger.warning(str(e).capitalize())
            try:
                if await self.is_moderator(
                    channel_dict.get("id", None), guild_snowflake, member.id
                ):
                    moderators.append(member)
            except commands.CheckFailure as e:
                self.__bot.logger.warning(str(e).capitalize())
        sysadmins_chunks = [
            sysadmins[i : i + chunk_size] for i in range(0, len(sysadmins), chunk_size)
        ]
        guild_owners_chunks = [
            guild_owners[i : i + chunk_size]
            for i in range(0, len(guild_owners), chunk_size)
        ]
        developers_chunks = [
            developers[i : i + chunk_size]
            for i in range(0, len(developers), chunk_size)
        ]
        administrators_chunks = [
            administrators[i : i + chunk_size]
            for i in range(0, len(administrators), chunk_size)
        ]
        coordinators_chunks = [
            coordinators[i : i + chunk_size]
            for i in range(0, len(coordinators), chunk_size)
        ]
        moderators_chunks = [
            moderators[i : i + chunk_size]
            for i in range(0, len(moderators), chunk_size)
        ]
        roles_chunks = [
            ("Sysadmins", sysadmins, sysadmins_chunks),
            ("Developers", developers, developers_chunks),
            ("Guild Owners", guild_owners, guild_owners_chunks),
            ("Administrators", administrators, administrators_chunks),
            ("Coordinators", coordinators, coordinators_chunks),
            ("Moderators", moderators, moderators_chunks),
        ]
        max_pages = max(len(c[2]) for c in roles_chunks)
        for page in range(max_pages):
            embed = discord.Embed(
                title=f"{self.__emoji.get_random_emoji()} Survey results for {channel_dict.get('name', None)}",
                description=f"Total surveyed: {len(channel_dict.get('object', None).members)}",
                color=discord.Color.blurple(),
            )
            for role_name, role_list, chunks in roles_chunks:
                chunk = chunks[page] if page < len(chunks) else []
                embed.add_field(
                    name=f"{role_name} ({len(chunk)}/{len(role_list)})",
                    value=", ".join(u.mention for u in chunk) if chunk else "*None*",
                    inline=False,
                )
            pages.append(embed)
        return pages

    async def toggle_moderator(self, channel_dict, default_kwargs, member_dict):
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(channel_dict.get("columns", None))
        updated_kwargs.update(member_dict.get("columns", None))
        moderator = await self.__database_factory.select(
            singular=True, **updated_kwargs
        )
        if moderator:
            await self.__database_factory.delete(**updated_kwargs)
            action = "revoked"
        else:
            moderator = self.MODEL(**updated_kwargs)
            await self.__database_factory.create(moderator)
            action = "granted"
        return (
            f"Moderator access for {member_dict.get('mention', None)} has been "
            f"{action} in {channel_dict.get('mention', None)}."
        )
