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

from copy import copy
from dataclasses import dataclass, field
from typing import Dict, List, Union

import discord
from discord.ext import commands

from vyrtuous.active_members import active_member_service
from vyrtuous.coordinator.coordinator import Coordinator


class NotCoordinator(commands.CheckFailure):
    def __init__(
        self,
        message="Member is not a coordinator in this channel.",
    ):
        super().__init__(message)


@dataclass
class CoordinatorDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, str]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


class CoordinatorService:
    __CHUNK_SIZE = 12
    MODEL = Coordinator
    coordinators = {}

    def __init__(
        self,
        *,
        active_member_service=None,
        author_service=None,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        emoji,
    ):
        self.__active_member_service = active_member_service
        self.__author_service = author_service
        self.__bot = bot
        self.__database_factory = copy(database_factory)
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

    async def is_coordinator_at_all_wrapper(self, context):
        return await self.is_coordinator_at_all(member_snowflake=context.author.id)

    async def is_coordinator_wrapper(self, context):
        return await self.is_coordinator(
            channel_snowflake=int(context.channel.id),
            guild_snowflake=int(context.guild.id),
            member_snowflake=int(context.member.id),
        )

    async def populate(self):
        coordinators = await self.__database_factory.select()
        for coordinator in coordinators:
            guild = self.__bot.get_guild(coordinator.guild_snowflake)
            if not guild:
                continue
            self.coordinators[coordinator.member_snowflake] = {
                "last_active": None,
                "name": coordinator.display_name,
            }

    async def build_dictionary(self, obj):
        coordinators = []
        dictionary = {}
        if isinstance(obj, discord.Guild):
            coordinators = await self.__database_factory.select(guild_snowflake=obj.id)
        elif isinstance(obj, discord.abc.GuildChannel):
            coordinators = await self.__database_factory.select(
                channel_snowflake=obj.id
            )
        elif isinstance(obj, discord.Member):
            coordinators = await self.__database_factory.select(member_snowflake=obj.id)
        else:
            coordinators = await self.__database_factory.select()
        if coordinators:
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
        return dictionary

    async def build_pages(self, is_at_home, obj):
        lines, pages = [], []

        obj_name = "All Servers"
        if not isinstance(obj, int):
            obj_name = obj.name
        else:
            member = self.__active_member_service.active_member.get(obj, None)
            if member:
                obj_name = member.get("name", None)
            else:
                return "No coordinators found."
        title = f"{self.__emoji.get_random_emoji()} Coordinators for {obj_name}"

        dictionary = await self.build_dictionary(obj=obj)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=CoordinatorDictionary, dictionary=dictionary
        )

        for guild_snowflake, guild_data in processed_dictionary.data.items():
            coord_n = 0
            field_count = 0
            lines = []
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, coordinator_dictionary in guild_data.get(
                "members"
            ).items():
                member = guild.get_member(member_snowflake)
                if member:
                    if not isinstance(obj, discord.Member):
                        lines.append(
                            f"**User:** {member.display_name} {member.mention}"
                        )
                        field_count += 1
                    elif not thumbnail:
                        embed.set_thumbnail(url=obj.display_avatar.url)
                        thumbnail = True
                else:
                    display_name = self.__active_member_service.active_members.get(
                        member_snowflake, None
                    )
                    if member:
                        display_name = member.get("name", None)
                        lines.append(f"**User:** {display_name} ({member_snowflake})")
                coord_n += 1
                for channel_snowflake, channel_dictionary in coordinator_dictionary.get(
                    "coordinators"
                ).items():
                    if not isinstance(obj, discord.abc.GuildChannel):
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
            original_description = embed.description or ""
            embed.description = f"**{original_description} ({coord_n})**"
            pages.append(embed)
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_members)
        return pages

    async def toggle_coordinator(self, channel, member_snowflake, display_name=None):
        coordinator = await self.__database_factory.select(
            singular=True,
            channel_snowflake=channel.id,
            member_snowflake=member_snowflake,
        )
        if coordinator:
            await self.__database_factory.delete(
                channel_snowflake=channel.id, member_snowflake=member_snowflake
            )
            del self.coordinators[member_snowflake]
            action = "revoked"
        else:
            coordinator = self.MODEL(
                channel_snowflake=channel.id,
                display_name=display_name,
                guild_snowflake=channel.guild.id,
                member_snowflake=member_snowflake,
            )
            await self.__database_factory.create(coordinator)
            self.coordinators.update({member_snowflake: {"name": display_name}})
            action = "granted"
        member = channel.guild.get_member(member_snowflake)
        if member:
            member_str = member.mention
        else:
            member_str = self.__active_member_service.active_members.get(
                member_snowflake, None
            ).get("name", None)
        return (
            f"Coordinator access has been {action} for {member_str} "
            f"in {channel.mention}."
        )

    async def migrate(self, kwargs):
        self.__database_factory.update(**kwargs)
