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
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.data import Data
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.users import moderator_service

MODEL = Data


class DataBuilder:
    def __init__(self, guild_snowflake: int, target_snowflake: int):
        self.__data = Data(
            guild_snowflake=guild_snowflake, target_snowflake=target_snowflake
        )

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
        target_snowflake=None,
    ) -> Self:
        if author_snowflake:
            self.__data = replace(self.__data, author_snowflake=author_snowflake)
        if channel_snowflake:
            self.__data = replace(self.__data, channel_snowflake=channel_snowflake)
        if guild_snowflake:
            self.__data = replace(self.__data, guild_snowflake=guild_snowflake)
        if role_snowflake:
            self.__data = replace(self.__data, role_snowflake=role_snowflake)
        if target_snowflake:
            self.__data = replace(self.__data, target_snowflake=target_snowflake)
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

    async def set_highest_roles(self, author_snowflake: int, member_snowflake: int):
        executor_highest_role = await moderator_service.resolve_highest_role_at_all(
            member_snowflake=int(author_snowflake),
        )
        target_highest_role = await moderator_service.resolve_highest_role_at_all(
            member_snowflake=int(member_snowflake),
        )
        self.__data = replace(self.__data, executor_highest_role=executor_highest_role)
        self.__data = replace(self.__data, target_highest_role=target_highest_role)
        return self

    def set_target(self, *, target: str | None = "user") -> Self:
        self.__data = replace(self.__data, target=target)
        return self

    async def create(self):
        database_factory = DatabaseFactory(MODEL)
        await database_factory.create(self.__data)


async def save_data(
    author_snowflake: int | None,
    channel_snowflake: int | None,
    duration_value: str | None,
    guild_snowflake: int,
    identifier: str,
    member_snowflake: int,
    reason: str,
    role_snowflake: int | None,
    target: str | None,
):
    data = DataBuilder(
        guild_snowflake=guild_snowflake, target_snowflake=member_snowflake
    )
    data.set_identifier(identifier=identifier)
    if target:
        data.set_target(target=target)
    if reason:
        data.set_reason(reason=reason)
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    data.set_snowflakes(guild_snowflake=guild_snowflake)
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
    data.set_counts(total_voice_members=total_voice_members)
    if channel_snowflake:
        channel = bot.get_channel(channel_snowflake)
        if channel is None:
            raise commands.ChannelNotFound(str(channel_snowflake))
        if isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
            current_channel_members = len(channel.members)
            data.set_counts(current_channel_members=current_channel_members)
    else:
        data.set_counts(current_channel_members=0)
    member = guild.get_member(member_snowflake)
    if member is None:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member is None:
            raise commands.MemberNotFound(str(member_snowflake))
        data.set_snowflakes(target_snowflake=member_snowflake)
    else:
        data.set_snowflakes(target_snowflake=member_snowflake)
    if author_snowflake:
        await data.set_highest_roles(
            author_snowflake=author_snowflake, member_snowflake=member_snowflake
        )
        data.set_snowflakes(author_snowflake=int(author_snowflake))
    if role_snowflake:
        role = guild.get_role(role_snowflake)
        if role is None:
            raise commands.RoleNotFound(str(role_snowflake))
        else:
            data.set_snowflakes(role_snowflake=role_snowflake)
    else:
        data.set_snowflakes(role_snowflake=None)
    duration_builder = DurationBuilder()
    if duration_value is not None:
        expires_at = duration_builder.parse(value=duration_value).to_expires_in()
        data.set_expires_at(expires_at=expires_at)
    else:
        expires_at = None
        data.set_expires_at(expires_at=expires_at)
    await data.create()
