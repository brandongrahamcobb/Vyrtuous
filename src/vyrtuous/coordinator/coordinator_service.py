from copy import copy

"""!/bin/python3
coordinator_service.py The purpose of this program is to extend Service to service the coordinator class.

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
from vyrtuous.coordinator.coordinator import Coordinator
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.sysadmin.sysadmin_service import SysadminService


class NotCoordinator(commands.CheckFailure):
    def __init__(
        self,
        message="Member is not a coordinator in this channel.",
    ):
        super().__init__(message)


class CoordinatorService:
    __CHUNK_SIZE = 7
    MODEL = Coordinator

    def __init__(
        self,
        *,
        author_service=None,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        emoji,
    ):
        self.__author_service = author_service
        self.__bot = bot
        self.__database_factory = database_factory
        self.__dictionary_service = dictionary_service
        self.__database_factory.model = self.MODEL
        self.__emoji = emoji

    async def is_coordinator(
        self, channel_snowflake: int, guild_snowflake: int, member_snowflake: int
    ) -> bool:
        coordinator = await self.__database_factory.select(
            channel_snowflake=int(channel_snowflake),
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
            singular=True,
        )
        if not coordinator:
            raise NotCoordinator
        return True

    async def is_coordinator_at_all(
        self,
        member_snowflake: int,
    ):
        coordinator = await self.__database_factory.select(
            member_snowflake=int(member_snowflake),
        )
        if not coordinator:
            raise NotCoordinator
        return True

    async def is_coordinator_at_all_wrapper(
        self,
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        member = self.__author_service.resolve_author(source=source)
        member_snowflake = member.id
        return await self.is_coordinator_at_all(member_snowflake=member_snowflake)

    async def is_coordinator_wrapper(
        self,
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        member = self.__author_service.resolve_author(source=source)
        member_snowflake = member.id
        return await self.is_coordinator(
            channel_snowflake=source.channel.id,
            guild_snowflake=source.guild.id,
            member_snowflake=int(member_snowflake),
        )

    async def build_clean_dictionary(self, is_at_home, where_kwargs):
        pages = []
        dictionary = {}
        coordinators = await self.__database_factory.select(
            singular=False, **where_kwargs
        )
        for coordinator in coordinators:
            dictionary.setdefault(coordinator.guild_snowflake, {"members": {}})
            dictionary[coordinator.guild_snowflake]["members"].setdefault(
                coordinator.member_snowflake, {"coordinators": {}}
            )
            dictionary[coordinator.guild_snowflake]["members"][
                coordinator.member_snowflake
            ]["coordinators"].setdefault(coordinator.channel_snowflake, {})
            dictionary[coordinator.guild_snowflake]["members"][
                coordinator.member_snowflake
            ]["coordinators"][coordinator.channel_snowflake].update(
                {"placeholder": "placeholder"}
            )
        skipped_guilds = self.__dictionary_service.generate_skipped_guilds(dictionary)
        skipped_members = self.__dictionary_service.generate_skipped_members(dictionary)
        cleaned_dictionary = self.__dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        # if is_at_home:
        #     if skipped_guilds:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_set_pages(
        #                 skipped=skipped_guilds,
        #                 title="Skipped Servers",
        #             )
        #         )
        #     if skipped_members:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_dict_pages(
        #                 skipped=skipped_members,
        #                 title="Skipped Members in Server",
        #             )
        #         )
        return cleaned_dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Coordinators {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )
        embed = discord.Embed(
            title=title, description="Default view", color=discord.Color.blue()
        )
        if not dictionary:
            pages = [embed]

        coord_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, coordinator_dictionary in guild_data.get(
                "members"
            ).items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                coord_n += 1
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in coordinator_dictionary.get(
                    "coordinators"
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
                        embed = flush_page(embed, pages, title, guild.name)
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
            pages[0].description = f"**({coord_n})**"
        return pages

    async def toggle_coordinator(self, channel_dict, default_kwargs, member_dict):
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(channel_dict.get("columns", None))
        updated_kwargs.update(member_dict.get("columns", None))
        coordinator = await self.__database_factory.select(
            singular=True, **updated_kwargs
        )
        if coordinator:
            await self.__database_factory.delete(**updated_kwargs)
            action = "revoked"
        else:
            coordinator = self.MODEL(**updated_kwargs)
            await self.__database_factory.create(coordinator)
            action = "granted"
        return (
            f"Coordinator access has been {action} for {member_dict.get('mention', None)} "
            f"in {channel_dict.get('mention', None)}."
        )
