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

from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.coordinator import Coordinator, NotCoordinator
from vyrtuous.db.database_factory import DatabaseFactory

MODEL = Coordinator


async def is_coordinator(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> bool:
    database_factory = DatabaseFactory(MODEL)
    coordinator = await database_factory.select(
        channel_snowflake=int(channel_snowflake),
        guild_snowflake=int(guild_snowflake),
        member_snowflake=int(member_snowflake),
        singular=True,
    )
    if not coordinator:
        raise NotCoordinator
    return True


async def is_coordinator_at_all(
    member_snowflake: int,
):
    database_factory = DatabaseFactory(MODEL)
    coordinator = await database_factory.select(
        member_snowflake=int(member_snowflake),
        singular=True,
    )
    if not coordinator:
        raise NotCoordinator
    return True


async def is_coordinator_at_all_wrapper(context):
    return await is_coordinator_at_all(member_snowflake=context.member_snowflake)


async def is_coordinator_wrapper(context):
    return await is_coordinator(
        channel_snowflake=int(context.channel_snowflake),
        guild_snowflake=int(context.guild_snowflake),
        member_snowflake=int(context.member_snowflake),
    )


async def toggle_coordinator(channel, member_snowflake):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(Coordinator)
    coordinator = await database_factory.select(
        singular=True,
        channel_snowflake=channel.id,
        member_snowflake=member_snowflake,
    )
    if coordinator:
        await database_factory.delete(
            channel_snowflake=channel.id, member_snowflake=member_snowflake
        )
        action = "revoked"
    else:
        coordinator = MODEL(
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
            member_snowflake=member_snowflake,
        )
        await database_factory.create(coordinator)
        action = "granted"
    member = channel.guild.get_member(member_snowflake)
    if member:
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member is None:
            raise commands.MemberNotFound(member_snowflake)
        member_str = simplified_member[0]
    return (
        f"Coordinator access has been {action} for {member_str} "
        f"in {channel.mention}."
    )
