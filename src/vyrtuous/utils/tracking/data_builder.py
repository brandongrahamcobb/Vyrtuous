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
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.data import Data
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.utils.users import moderator_service

MODEL = Data


class DataBuilder:
    def __init__(self, target_snowflake: int):
        self.__data = Data(target_snowflake=target_snowflake)

    def set_counts(
        self,
        *,
        current_channel_members=0,
        total_guild_members=0,
        online_members=0,
        total_voice_members=0,
    ) -> Self:
        self.__data = replace(
            self.__data, current_channel_members=current_channel_members
        )
        self.__data = replace(self.__data, total_guild_members=total_guild_members)
        self.__data = replace(self.__data, online_members=online_members)
        self.__data = replace(self.__data, total_voice_members=total_voice_members)
        return self

    def set_snowflakes(
        self,
        *,
        author_snowflake=None,
        channel_snowflake=None,
        guild_snowflake=None,
        role_snowflake=None,
    ) -> Self:
        if author_snowflake:
            self.__data = replace(self.__data, author_snowflake=author_snowflake)
        if channel_snowflake:
            self.__data = replace(self.__data, channel_snowflake=channel_snowflake)
        if guild_snowflake:
            self.__data = replace(self.__data, guild_snowflake=guild_snowflake)
        if role_snowflake:
            self.__data = replace(self.__data, role_snowflake=role_snowflake)
        return self

    def set_expires_at(self, *, expires_at=None) -> Self:
        self.__data = replace(
            self.__data,
            expires_at=expires_at if expires_at else datetime.now(timezone.utc),
        )
        return self

    def set_identifier(self, *, identifier=None) -> Self:
        self.__data = replace(self.__data, identifier=identifier)
        return self

    def set_reason(self, *, reason="No reason provided") -> Self:
        self.__data = replace(self.__data, reason=reason)
        return self

    async def set_highest_roles(
        self,
        *,
        author_snowflake=None,
        channel_snowflake=None,
        guild_snowflake=None,
        member_snowflake=None,
    ):
        if author_snowflake and channel_snowflake and guild_snowflake:
            executor_highest_role = await moderator_service.resolve_highest_role(
                channel_snowflake=int(channel_snowflake),
                guild_snowflake=int(guild_snowflake),
                member_snowflake=int(author_snowflake),
            )
            self.__data = replace(
                self.__data, executor_highest_role=executor_highest_role
            )
            if member_snowflake:
                target_highest_role = await moderator_service.resolve_highest_role(
                    channel_snowflake=int(channel_snowflake),
                    guild_snowflake=int(guild_snowflake),
                    member_snowflake=int(member_snowflake),
                )
                self.__data = replace(
                    self.__data, target_highest_role=target_highest_role
                )
        return self

    def set_target(self, *, target: str | None = "user") -> Self:
        if target:
            self.__data = replace(self.__data, target=target)
        return self

    async def create(self):
        database_factory = DatabaseFactory(MODEL)
        await database_factory.create(self.__data)


async def save_data(
    author_snowflake: int | None,
    channel_snowflake: int | None,
    duration: DurationObject | None,
    guild_snowflake: int | None,
    identifier: str,
    member_snowflake: int,
    reason: str,
    role_snowflake: int | None,
    target: str | None,
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    data = DataBuilder(target_snowflake=member_snowflake)
    data.set_snowflakes(
        author_snowflake=author_snowflake,
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        role_snowflake=role_snowflake,
    )
    data.set_identifier(identifier=identifier)
    data.set_target(target=target)
    data.set_reason(reason=reason)
    await data.set_highest_roles(
        author_snowflake=author_snowflake,
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
    )
    if guild_snowflake:
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(guild_snowflake))
        total_guild_members = sum(1 for member in guild.members if not member.bot)
        data.set_counts(total_guild_members=total_guild_members)
        online_members = sum(
            1
            for member in guild.members
            if not member.bot and member.status != discord.Status.offline
        )
        data.set_counts(online_members=online_members)
        total_voice_members = sum(
            len([member for member in channel.members if not member.bot])
            for channel in guild.voice_channels
        )
        data.set_counts(
            total_guild_members=total_guild_members,
            online_members=online_members,
            total_voice_members=total_voice_members,
        )
        if channel_snowflake:
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                raise commands.ChannelNotFound(str(channel_snowflake))
            if isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
                current_channel_members = len(channel.members)
                data.set_counts(current_channel_members=current_channel_members)
        else:
            data.set_counts(current_channel_members=0)
    duration_builder = DurationBuilder()
    if duration is not None:
        expires_at = duration_builder.load(duration=duration).to_expires_in()
        data.set_expires_at(expires_at=expires_at)
    else:
        expires_at = None
        data.set_expires_at(expires_at=expires_at)
    await data.create()
