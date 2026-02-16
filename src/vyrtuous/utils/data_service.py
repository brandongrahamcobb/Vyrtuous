"""!/bin/python3

data.py The purpose of this program is to manage statistics of Vyrtuous.

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

from dataclasses import replace
from datetime import datetime, timezone
from typing import Self

import discord

from vyrtuous.utils.data import Data


class DataService:
    def __init__(self, *, duration_service, moderator_service=None):
        self.__data = Data()
        self.__duration_service = duration_service
        self.__moderator_service = moderator_service

    def set_counts(
        self,
        *,
        current_channel_members=0,
        total_guild_members=0,
        online_members=0,
        total_voice_members=0,
    ) -> Self:
        replace(self.__data, current_channel_members=current_channel_members)
        replace(self.__data, total_guild_members=total_guild_members)
        replace(self.__data, online_members=online_members)
        replace(self.__data, total_voice_members=total_voice_members)
        return self

    def set_snowflakes(
        self,
        *,
        author_snowflake=None,
        channel_snowflake=None,
        guild_snowflake=None,
        target_snowflake=None,
    ) -> Self:
        if author_snowflake:
            replace(self.__data, author_snowflake=author_snowflake)
        if channel_snowflake:
            replace(self.__data, channel_snowflake=channel_snowflake)
        if guild_snowflake:
            replace(self.__data, guild_snowflake=guild_snowflake)
        if target_snowflake:
            replace(self.__data, target_snowflake=target_snowflake)
        return self

    def set_expires_at(self, *, expires_at=None) -> Self:
        replace(
            self.__data,
            expires_at=expires_at if expires_at else datetime.now(timezone.utc),
        )
        return self

    def set_identifier(self, *, identifier=None) -> Self:
        replace(self.__data, identiifer=identifier)
        return self

    def set_reason(self, *, reason="No reason provided") -> Self:
        replace(self.__data, reason=reason)
        return self

    def set_is_modification(self, *, is_modification=False) -> Self:
        replace(self.__data, is_modification=is_modification)
        return self

    async def set_highest_roles(
        self,
        *,
        author: discord.Member | None = None,
        target: discord.Member | None = None,
    ) -> Self:
        executor_role = await self.__moderator_service.resolve_highest_role_at_all(
            member_snowflake=int(author.id),
        )
        target_role = await self.__moderator_service.resolve_highest_role_at_all(
            member_snowflake=int(target.id),
        )
        replace(self.__data, executor_role=executor_role)
        replace(self.__data, target_role=target_role)
        return self

    async def save_data(
        self,
        author: discord.Member,
        channel: discord.abc.GuildChannel,
        identifier: str,
        target: discord.Member,
        duration: str | None = None,
        is_modification: bool = False,
        reason: str = "No reason provided",
    ):
        duration = self.__duration_service.parse(duration)
        if duration is not None:
            expires_at = datetime.now(timezone.utc) + duration.to_timedelta()
        else:
            expires_at = None
        current_channel_members = len(channel.members)
        total_guild_members = sum(
            1 for member in channel.guild.members if not member.bot
        )
        online_members = sum(
            1
            for member in channel.guild.members
            if not member.bot and member.status != discord.Status.offline
        )
        total_voice_members = sum(
            len([member for member in channel.members if not member.bot])
            for channel in channel.guild.voice_channels
        )
        # fmt: off
        self.set_counts(current_channel_members=current_channel_members) \
            .set_counts(total_guild_members=total_guild_members) \
            .set_counts(online_members=online_members) \
            .set_counts(total_voice_members=total_voice_members) \
            .set_snowflakes(channel_snowflake=int(channel.id)) \
            .set_snowflakes(guild_snowflake=int(channel.guild.id)) \
            .set_snowflakes(author_snowflake=int(author.id)) \
            .set_snowflakes(target_snowflake=int(target.id)) \
            .set_expires_at(expires_at=expires_at) \
            .set_identifier(identifier=identifier) \
            .set_is_modification(is_modification=is_modification) \
            .set_reason(reason=reason)
        # fmt: on
        await self.__database_factory.create(self.data)
