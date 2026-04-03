"""!/bin/python3
guild_owner_service.py The purpose of this program is to extend Service to service the guild owner class.

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
from enum import member
from typing import Union

import discord
from discord.ext import commands

from vyrtuous.active_members import active_member_service
from vyrtuous.owner.guild_owner import GuildOwner


class NotGuildOwner(commands.CheckFailure):
    def __init__(
        self,
        message="Member is not a guild owner in this server.",
    ):
        super().__init__(message)


class GuildOwnerService:
    MODEL = GuildOwner
    guild_owners = {}

    def __init__(
        self,
        *,
        active_member_service=None,
        author_service=None,
        bot=None,
        database_factory=None,
    ):
        self.__active_member_service = active_member_service
        self.__author_service = author_service
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL

    async def populate(self):
        guild_owners = await self.__database_factory.select()
        for guild_owner in guild_owners:
            guild = self.__bot.get_guild(guild_owner.guild_snowflake)
            if not guild:
                continue
            self.guild_owners[guild_owner.member_snowflake] = {
                "last_active": None,
                "name": guild_owner.display_name,
            }

    async def update_guild_owners(self):
        for guild in self.__bot.guilds:
            guild_owner = await self.__database_factory.select(
                guild_snowflake=guild.id, singular=True
            )
            if guild_owner and guild_owner.member_snowflake == guild.owner_id:
                self.__bot.logger.info(
                    f"Guild owner ({guild_owner.member_snowflake}) already in the db."
                )
                continue
            else:
                guild_owner = self.MODEL(
                    guild_snowflake=guild.id, member_snowflake=guild.owner_id
                )
                await self.__database_factory.create(guild_owner)
                member = guild.get_member(guild.owner.id)
                self.guild_owners.update(
                    {guild.owner_id: {"name": member.display_name}}
                )
                self.__bot.logger.info(
                    f"Guild owner ({guild_owner.member_snowflake}) added to the db."
                )

    async def add_guild_owner(self, guild_snowflake, member_snowflake):
        guild_owner = self.MODEL(
            guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
        )
        await self.__database_factory.create(guild_owner)
        member = guild.get_member(member_snowflake)
        self.guild_owners.update({member_snowflake: {"name": member.display_name}})
        self.__bot.logger.info(f"Guild owner ({member_snowflake}) added.")

    async def remove_guild_owner(self, guild_snowflake, member_snowflake):
        await self.__database_factory.delete(
            guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
        )
        del self.guild_owners[member_snowflake]
        self.__bot.logger.info(f"Guild owner ({member_snowflake}) removed.")

    async def is_guild_owner_wrapper(
        self,
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        member = self.__author_service.resolve_author(source=source)
        member_snowflake = member.id
        return await self.is_guild_owner(
            guild_snowflake=source.guild.id, member_snowflake=int(member_snowflake)
        )

    async def is_guild_owner_at_all(
        self,
        member_snowflake: int,
    ):
        for guild in self.__bot.guilds:
            if guild and guild.owner_id == member_snowflake:
                return True
        raise NotGuildOwner

    async def is_guild_owner(self, guild_snowflake: int, member_snowflake: int) -> bool:
        guild = self.__bot.get_guild(guild_snowflake)
        if guild and guild.owner_id == member_snowflake:
            return True
        raise NotGuildOwner
