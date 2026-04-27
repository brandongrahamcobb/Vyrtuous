"""!/bin/python3
alias_context.py The purpose of this program is to classify Alias context.

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
from datetime import datetime

from discord.ext import commands
import discord

from vyrtuous.alias.alias import Alias
from vyrtuous.alias.alias_service import AliasService


class AliasContext:
    MODEL = Alias

    def __init__(
        self,
        *,
        active_member_service=None,
        bot=None,
        cap_service=None,
        content=None,
        database_factory=None,
        default_ctx=None,
        dictionary_service=None,
        duration_builder=None,
        emoji=None,
        moderator_service=None,
    ):
        self.alias = None
        self.__active_member_service = active_member_service
        self.__alias_name: str
        self.__bot = bot
        self.__cap_service = cap_service
        self.__content = content
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__duration_builder = duration_builder
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji
        self.__alias_service = AliasService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__d_ctx = default_ctx
        self.__moderator_service = moderator_service
        self.__args = []
        self.__kwargs: dict[str, tuple] = {}
        self.category: str | None = None
        self.channel: discord.abc.GuildChannel | None = None
        self.guild: discord.Guild | None = None
        self.member_snowflake: int | None = None
        self.role: discord.Role | None = None
        self.expires_in: datetime | None = None
        self.duration_value: str | None = None
        self.reason = "No reason provided"
        self.record = None

    async def setup(self):
        self._message_to_args()
        self._alias_name_from_args()
        self.alias = await self._populate_alias()
        if not self.alias:
            return False
        self._fill_map()
        await self._convert_args_to_values()
        return True

    def _message_to_args(self) -> None:
        self.__args = (
            self.__content[len(self.__bot.config["discord_command_prefix"]) :]
            .strip()
            .split()
        )

    def _fill_map(self) -> None:
        map = self.alias.ARGS_MAP
        sorted_args = sorted(map.items(), key=lambda x: x[1])
        for i, (key, pos) in enumerate(sorted_args):
            if i == len(sorted_args) - 1:
                value = (
                    " ".join(str(a) for a in self.__args[pos - 1 :])
                    if len(self.__args) >= pos
                    else ""
                )
            else:
                value = str(self.__args[pos - 1]) if len(self.__args) >= pos else ""
            self.__kwargs[key] = (pos, value)

    def _alias_name_from_args(self):
        if not self.__args:
            return
        self.__alias_name = self.__args[0]

    async def _populate_alias(self):
        alias_entry = await self.__database_factory.select(
            alias_name=self.__alias_name,
            guild_snowflake=self.__d_ctx.guild.id,
            singular=True,
        )
        if not alias_entry:
            return None
        self.guild = self.__bot.get_guild(int(alias_entry.guild_snowflake))
        self.channel = self.guild.get_channel(int(alias_entry.channel_snowflake))
        if getattr(alias_entry, "role_snowflake"):
            self.role = self.guild.get_role(int(alias_entry.role_snowflake))
        alias = self.__alias_service.alias_category_to_alias(
            alias_category=alias_entry.category
        )
        self.category = alias_entry.category
        self.record = alias.record
        return alias

    async def _convert_args_to_values(self):
        for field, tuple in self.__kwargs.items():
            value = tuple[1]
            if field == "duration":
                if not value:
                    self.duration_value = "8h"
                else:
                    self.duration_value = value
                if await self.__cap_service.assertion(
                    ctx=self,
                    default_ctx=self.__d_ctx,
                ):
                    await self.__moderator_service.check_minimum_role(
                        channel_snowflake=self.channel.id,
                        guild_snowflake=self.guild.id,
                        member_snowflake=self.__d_ctx.author.id,
                        lowest_role="Coordinator",
                    )
                self.expires_in = self.__duration_builder.parse(
                    self.duration_value
                ).to_expires_in()
            elif field == "member":
                self.member_snowflake = int(value.replace("<@", "").replace(">", ""))
                self.member = self.__d_ctx.guild.get_member(self.member_snowflake)
                if not self.member:
                    self.__bot.logger.info(
                        self.__active_member_service.active_members.get(
                            self.member_snowflake, None
                        )
                    )
                    member = self.__active_member_service.active_members.get(
                        self.member_snowflake, None
                    )
                    if member:
                        display_name = member.get("name", None)
                    else:
                        raise commands.MemberNotFound(self.member_snowflake)
                else:
                    display_name = self.member.display_name
                self.display_name = display_name
            elif field == "reason":
                if not value:
                    self.reason = "No reason provided."
                else:
                    self.reason = value
