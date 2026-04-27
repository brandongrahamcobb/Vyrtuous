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
