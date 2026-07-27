"""!/bin/python3
coordinator_service.py The purpose of this program is to extend Service to service the coordinator class.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.coordinator import Coordinator, NotCoordinator
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationObject
from vyrtuous.utils.tracking import data_builder, stream_service

MODEL = Coordinator


async def is_coordinator(
    guild_snowflake: int, member_snowflake: int, *, channel_snowflake: int | None
) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    if channel_snowflake is None:
        coordinator = await database_factory.select(
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
            singular=True,
        )
        if not coordinator:
            raise NotCoordinator(
                channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake
            )
        return True
    else:
        coordinator = await database_factory.select(
            channel_snowflake=int(channel_snowflake),
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
            singular=True,
        )
        if not coordinator:
            raise NotCoordinator(
                channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake
            )
        return True


async def toggle_coordinator(
    author_snowflake: int,
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
    message_snowflake: int,
    message_channel_snowflake: int,
) -> str:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(Coordinator)
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(str(channel_snowflake))
    member = channel.guild.get_member(member_snowflake)
    if member:
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member is None:
            raise commands.MemberNotFound(str(member_snowflake))
        member_str = simplified_member[0]
    coordinator = await database_factory.select(
        singular=True,
        channel_snowflake=channel_snowflake,
        member_snowflake=member_snowflake,
    )
    if coordinator:
        await database_factory.delete(
            channel_snowflake=channel_snowflake, member_snowflake=member_snowflake
        )
        await log_xcoord(
            author_snowflake=author_snowflake,
            channel_snowflake=channel_snowflake,
            display=True,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake,
            message_channel_snowflake=message_channel_snowflake,
        )
        action = "revoked"
    else:
        coordinator = MODEL(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        await database_factory.create(coordinator)
        await log_coord(
            author_snowflake=author_snowflake,
            channel_snowflake=channel_snowflake,
            display=True,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake,
            message_channel_snowflake=message_channel_snowflake,
        )
        action = "granted"
    return (
        f"Coordinator access has been {action} for {member_str} "
        f"in {channel.mention}."
    )


async def log_coord(
    author_snowflake: int | None,
    channel_snowflake: int,
    display: bool,
    guild_snowflake: int,
    member_snowflake: int,
    message_snowflake: int | None,
    message_channel_snowflake: int | None,
):
    duration = DurationObject(number=0, prefix="", sign=1, unit="")
    is_channel_scope = None
    reason = None
    role_snowflake = None
    target = None
    await data_builder.save_data(
        author_snowflake=author_snowflake or None,
        channel_snowflake=channel_snowflake,
        duration=duration,
        guild_snowflake=guild_snowflake,
        identifier="coord",
        member_snowflake=member_snowflake,
        reason=reason or "No reason provided.",
        role_snowflake=role_snowflake or None,
        target=target or None,
    )
    if display:
        await stream_service.send_log(
            author_snowflake=author_snowflake or None,
            channel_snowflake=channel_snowflake,
            identifier="coord",
            duration=duration,
            guild_snowflake=guild_snowflake,
            is_channel_scope=is_channel_scope,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake or None,
            message_channel_snowflake=message_channel_snowflake or None,
            reason=reason or "No reason provided.",
            role_snowflake=role_snowflake or None,
            target=target or None,
        )


async def log_xcoord(
    author_snowflake: int | None,
    channel_snowflake: int,
    display: bool,
    guild_snowflake: int,
    member_snowflake: int,
    message_snowflake: int | None,
    message_channel_snowflake: int | None,
):
    duration = DurationObject(number=0, prefix="", sign=1, unit="")
    is_channel_scope = None
    reason = None
    role_snowflake = None
    target = None
    await data_builder.save_data(
        author_snowflake=author_snowflake or None,
        channel_snowflake=channel_snowflake,
        duration=duration,
        guild_snowflake=guild_snowflake,
        identifier="xcoord",
        member_snowflake=member_snowflake,
        reason=reason or "No reason provided.",
        role_snowflake=role_snowflake or None,
        target=target or None,
    )
    if display:
        await stream_service.send_log(
            author_snowflake=author_snowflake or None,
            channel_snowflake=channel_snowflake,
            identifier="xcoord",
            duration=duration,
            guild_snowflake=guild_snowflake,
            is_channel_scope=is_channel_scope,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake or None,
            message_channel_snowflake=message_channel_snowflake or None,
            reason=reason or "No reason provided.",
            role_snowflake=role_snowflake or None,
            target=target or None,
        )
