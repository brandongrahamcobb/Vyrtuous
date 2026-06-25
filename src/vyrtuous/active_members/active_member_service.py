"""!/bin/python3
active_member.py The purpose of this program is to service active members.

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
from datetime import datetime, timezone
import discord
from vyrtuous.active_members.active_member import ActiveMember


class ActiveMemberService:

    __MODEL = ActiveMember
    active_members = {}

    def __init__(self, *, bot=None, database_factory=None):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.__MODEL

    async def is_active(self, member_snowflake):
        if member_snowflake in self.active_members:
            return True
        return False

    async def populate(self):
        active_members = await self.__database_factory.select()
        for member in active_members:
            await self.update_active_member(
                last_active=member.last_active,
                guild_snowflake=member.guild_snowflake,
                member_snowflake=member.member_snowflake,
                name=member.display_name,
            )

    async def update_active_member(
        self, guild_snowflake, member_snowflake, name, last_active=None
    ):
        if last_active is None:
            last_active = datetime.now(timezone.utc)
        self.active_members.update(
            {
                member_snowflake: {
                    "last_active": last_active,
                    "name": name,
                    "guild_snowflake": guild_snowflake,
                    "id": member_snowflake,
                }
            }
        )

    async def save_active_members(self):
        active_members = await self.__database_factory.select()
        member_snowflakes = [
            active_member.member_snowflake for active_member in active_members
        ]
        for member_snowflake in self.active_members:
            if member_snowflake not in member_snowflakes:
                active_member = ActiveMember(
                    guild_snowflake=self.active_members[member_snowflake][
                        "guild_snowflake"
                    ],
                    last_active=self.active_members[member_snowflake]["last_active"],
                    member_snowflake=member_snowflake,
                    display_name=self.active_members[member_snowflake]["name"],
                )
                await self.__database_factory.create(active_member)

    async def remove_inactive_members(self, guild):
        for member in self.active_members:
            member_snowflakes = [
                active_member.member_snowflake for active_member in self.active_members
            ]
            for member_snowflake in member_snowflakes:
                if self.active_members[member_snowflake]["guild_snowflake"] == guild.id:
                    del self.active_members[member]
