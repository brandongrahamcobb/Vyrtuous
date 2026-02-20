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
from datetime import datetime, timezone
from typing import Dict, Tuple

import discord

from vyrtuous.alias.alias import Alias
from vyrtuous.alias.alias_service import AliasService


class AliasContext:
    MODEL = Alias

    def __init__(
        self,
        *,
        bot=None,
        cap_service=None,
        database_factory=None,
        dictionary_service=None,
        duration_service=None,
        emoji=None,
        message: discord.Message | None = None,
        moderator_service=None,
    ):
        self.alias = None
        self.alias_name: str
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji
        self.__alias_service = AliasService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.args = []
        self.expires_in = None
        self.kwargs: Dict[str, Tuple[int, str]] = {}
        self.source_kwargs: Dict[str, int] = {}
        self.message = message
        self.reason = "No reason provided"
        self.target_channel_snowflake: int
        self.target_member_snowflake: int
        self.target_role_snowflake: int
        self.source_guild_snowflake: int
        self.source_channel_snowflake = message.channel.id
        self.source_guild_snowflake = message.guild.id
        self.source_member_snowflake = message.author.id
        self.record = None
        self.__cap_service = cap_service
        self.__duration_service = duration_service
        self.__moderator_service = moderator_service

    async def setup(self):
        self.build_source_kwargs()
        self.message_to_args()
        self.alias_name_from_args()
        await self.populate_alias()
        self.fill_map()
        await self.convert_args_to_values()

    def build_source_kwargs(self):
        self.source_kwargs = {
            "channel_snowflake": self.source_channel_snowflake,
            "guild_snowflake": self.source_guild_snowflake,
            "member_snowflake": self.source_member_snowflake,
        }

    def message_to_args(self) -> None:
        self.args = (
            self.message.content[len(self.__bot.config["discord_command_prefix"]) :]
            .strip()
            .split()
        )

    def fill_map(self) -> None:
        map = self.alias.ARGS_MAP
        sorted_args = sorted(map.items(), key=lambda x: x[1])
        for i, (key, pos) in enumerate(sorted_args):
            if i == len(sorted_args) - 1:
                value = (
                    " ".join(str(a) for a in self.args[pos - 1 :])
                    if len(self.args) >= pos
                    else ""
                )
            else:
                value = str(self.args[pos - 1]) if len(self.args) >= pos else ""
            self.kwargs[key] = (pos, value)

    def alias_name_from_args(self):
        if not self.args:
            return
        self.alias_name = self.args[0]

    async def populate_alias(self):
        alias_entry = await self.__database_factory.select(
            alias_name=self.alias_name,
            guild_snowflake=self.source_guild_snowflake,
            singular=True,
        )
        if not alias_entry:
            return
        self.target_channel_snowflake = int(alias_entry.channel_snowflake)
        if getattr(alias_entry, "role_snowflake"):
            self.target_role_snowflake = int(alias_entry.role_snowflake)
        alias = self.__alias_service.alias_category_to_alias(
            alias_category=alias_entry.category
        )
        self.alias = alias
        self.record = alias.record

    async def convert_args_to_values(self):
        for field, tuple in self.kwargs.items():
            value = tuple[1]
            if field == "duration":
                if not value:
                    duration = self.__duration_service.parse("8h")
                else:
                    duration = self.__duration_service.parse(value)
                if await self.__cap_service.assert_duration_exceeds_cap(
                    duration=duration,
                    source_kwargs=self.source_kwargs,
                    category=self.alias.category,
                ):
                    await self.__moderator_service.check(
                        channel_snowflake=self.target_channel_snowflake,
                        guild_snowflake=self.source_guild_snowflake,
                        member_snowflake=self.target_member_snowflake,
                        lowest_role="Coordinator",
                    )
                self.expires_in = (
                    None
                    if duration.number == 0
                    else self.__duration_service.to_expires_in(duration)
                )
            elif field == "member":
                self.target_member_snowflake = int(value)
            elif field == "reason":
                self.reason = value
