"""!/bin/python3
moderator_service.py The purpose of this program is to extend Service to service the moderator class.

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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.guild_owner import NotGuildOwner
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.administrator import NotAdministrator
from vyrtuous.db.coordinator import NotCoordinator
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.developer import NotDeveloper
from vyrtuous.db.moderator import Moderator, NotModerator
from vyrtuous.inc.helpers import resolve_author
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.users import (
    administrator_service,
    coordinator_service,
    developer_service,
    guild_owner_service,
    sysadmin_service,
)

MODEL = Moderator
PERMISSION_TYPES = [
    "Everyone",
    "Moderator",
    "Coordinator",
    "Administrator",
    "Guild Owner",
    "Developer",
    "Sysadmin",
]


class HasEqualOrLowerRole(commands.CheckFailure):
    def __init__(self, target_rank=str):
        super().__init__(
            message=f"You may not execute this command on this `{target_rank}` because they have equal or higher role than you in this channel/server."
        )


async def is_moderator_wrapper(context) -> bool:
    return await is_moderator(
        channel_snowflake=int(context.channel_snowflake),
        guild_snowflake=int(context.guild_snowflake),
        member_snowflake=int(context.member_snowflake),
    )


async def is_moderator(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    moderator = await database_factory.select(
        channel_snowflake=int(channel_snowflake),
        guild_snowflake=int(guild_snowflake),
        member_snowflake=int(member_snowflake),
        singular=True,
    )
    if not moderator:
        raise NotModerator
    return True


async def is_moderator_at_all_wrapper(context) -> bool:
    return await is_moderator_at_all(member_snowflake=context.author.id)


async def is_moderator_at_all(
    member_snowflake: int,
) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    moderator = await database_factory.select(
        member_snowflake=int(member_snowflake), singular=True
    )
    if not moderator:
        raise NotModerator
    return True


async def survey(channel) -> list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    chunk_size, pages = 7, []
    (
        sysadmins,
        developers,
        guild_owners,
        administrators,
        coordinators,
        moderators,
    ) = ([], [], [], [], [], [])

    for member in channel.members:
        try:
            if await sysadmin_service.is_sysadmin(member_snowflake=member.id):
                sysadmins.append(member)
        except commands.CheckFailure as e:
            bot.logger.warning(str(e).capitalize())
        try:
            if await developer_service.is_developer(member_snowflake=member.id):
                developers.append(member)
        except commands.CheckFailure as e:
            bot.logger.warning(str(e).capitalize())
        try:
            if await guild_owner_service.is_guild_owner(
                guild_snowflake=channel.guild.id, member_snowflake=member.id
            ):
                guild_owners.append(member)
        except commands.CheckFailure as e:
            bot.logger.warning(str(e).capitalize())
        try:
            if await administrator_service.is_administrator(
                guild_snowflake=channel.guild.id, member_snowflake=member.id
            ):
                administrators.append(member)
        except commands.CheckFailure as e:
            bot.logger.warning(str(e).capitalize())
        try:
            if await coordinator_service.is_coordinator(
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                member_snowflake=member.id,
            ):
                coordinators.append(member)
        except commands.CheckFailure as e:
            bot.logger.warning(str(e).capitalize())
        try:
            if await is_moderator(
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                member_snowflake=member.id,
            ):
                moderators.append(member)
        except commands.CheckFailure as e:
            bot.logger.warning(str(e).capitalize())
    sysadmins_chunks = [
        sysadmins[i : i + chunk_size] for i in range(0, len(sysadmins), chunk_size)
    ]
    guild_owners_chunks = [
        guild_owners[i : i + chunk_size]
        for i in range(0, len(guild_owners), chunk_size)
    ]
    developers_chunks = [
        developers[i : i + chunk_size] for i in range(0, len(developers), chunk_size)
    ]
    administrators_chunks = [
        administrators[i : i + chunk_size]
        for i in range(0, len(administrators), chunk_size)
    ]
    coordinators_chunks = [
        coordinators[i : i + chunk_size]
        for i in range(0, len(coordinators), chunk_size)
    ]
    moderators_chunks = [
        moderators[i : i + chunk_size] for i in range(0, len(moderators), chunk_size)
    ]
    roles_chunks = [
        ("Sysadmins", sysadmins, sysadmins_chunks),
        ("Developers", developers, developers_chunks),
        ("Guild Owners", guild_owners, guild_owners_chunks),
        ("Administrators", administrators, administrators_chunks),
        ("Coordinators", coordinators, coordinators_chunks),
        ("Moderators", moderators, moderators_chunks),
    ]
    max_pages = max(len(c[2]) for c in roles_chunks)
    for page in range(max_pages):
        embed = discord.Embed(
            title=f"{emojis.get_random_emoji()} Survey results for {channel.name}",
            description=f"Total surveyed: {len(channel.members)}",
            color=discord.Color.blurple(),
        )
        for role_name, role_list, chunks in roles_chunks:
            chunk = chunks[page] if page < len(chunks) else []
            embed.add_field(
                name=f"{role_name} ({len(chunk)}/{len(role_list)})",
                value=", ".join(u.mention for u in chunk) if chunk else "*None*",
                inline=False,
            )
        pages.append(embed)
    return pages


async def toggle_moderator(channel, member_snowflake: int) -> str:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    moderator = await database_factory.select(
        channel_snowflake=int(channel.id),
        member_snowflake=int(member_snowflake),
        singular=True,
    )
    if moderator:
        await database_factory.delete(
            channel_snowflake=int(channel.id),
            member_snowflake=int(member_snowflake),
        )
        action = "revoked"
    else:
        moderator = MODEL(
            channel_snowflake=int(channel.id),
            guild_snowflake=int(channel.guild.id),
            member_snowflake=int(member_snowflake),
        )
        await database_factory.create(moderator)
        action = "granted"
    member = channel.guild.get_member(member_snowflake)
    if member:
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member:
            member_str = simplified_member[0]
        else:
            raise commands.MemberNotFound(str(member_snowflake))
    return (
        f"Moderator access for {member_str} has been " f"{action} in {channel.mention}."
    )


async def check_minimum_role(
    channel_snowflake,
    guild_snowflake,
    member_snowflake,
    lowest_role: str,
) -> str:
    verifications = (
        ("Sysadmin", sysadmin_service.is_sysadmin),
        ("Developer", developer_service.is_developer),
        ("Guild Owner", guild_owner_service.is_guild_owner),
        ("Administrator", administrator_service.is_administrator),
        ("Coordinator", coordinator_service.is_coordinator),
        ("Moderator", is_moderator),
    )
    passed_lowest = False
    for role_name, verify in verifications:
        if role_name == lowest_role:
            passed_lowest = True
        try:
            if role_name in ("Sysadmin", "Developer"):
                if await verify(member_snowflake=int(member_snowflake)):
                    return role_name
            elif role_name in ("Guild Owner", "Administrator"):
                if await verify(
                    guild_snowflake=int(guild_snowflake),
                    member_snowflake=int(member_snowflake),
                ):
                    return role_name
            else:
                if await verify(
                    channel_snowflake=int(channel_snowflake),
                    guild_snowflake=int(guild_snowflake),
                    member_snowflake=int(member_snowflake),
                ):
                    return role_name
        except commands.CheckFailure:
            if lowest_role is not None and passed_lowest:
                raise
    return "Everyone"


async def check_minimum_role_at_all(
    guild_snowflake,
    member_snowflake,
    lowest_role: str,
) -> str:
    verifications = (
        ("Sysadmin", sysadmin_service.is_sysadmin),
        ("Developer", developer_service.is_developer),
        ("Guild Owner", guild_owner_service.is_guild_owner),
        ("Administrator", administrator_service.is_administrator),
        ("Coordinator", coordinator_service.is_coordinator_at_all),
        ("Moderator", is_moderator_at_all),
    )
    passed_lowest = False
    for role_name, verify in verifications:
        if role_name == lowest_role:
            passed_lowest = True
        try:

            if role_name in ("Sysadmin", "Developer"):
                if await verify(member_snowflake=int(member_snowflake)):
                    return role_name
            else:
                if await verify(
                    guild_snowflake=int(guild_snowflake),
                    member_snowflake=int(member_snowflake),
                ):
                    return role_name
        except commands.CheckFailure:
            if lowest_role is not None and passed_lowest:
                raise
        if role_name == lowest_role:
            passed_lowest = True
    return "Everyone"


async def has_equal_or_lower_role(
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
    target_member_snowflake: int,
) -> bool:
    sender_name = await resolve_highest_role(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
    )
    sender_rank = PERMISSION_TYPES.index(sender_name)
    target_name = await resolve_highest_role(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=target_member_snowflake,
    )
    target_rank = PERMISSION_TYPES.index(target_name)
    compare_ranks(sender_rank=sender_rank, target_rank=target_rank)
    return sender_name


async def resolve_highest_role(
    channel_snowflake: int,
    member_snowflake: int,
    guild_snowflake: int,
) -> str:
    bot: DiscordBot = DiscordBot.get_instance()
    try:
        if await sysadmin_service.is_sysadmin(member_snowflake=int(member_snowflake)):
            return "Sysadmin"
    except sysadmin_service.NotSysadmin as e:
        bot.logger.warning(str(e).capitalize())
    try:
        if await developer_service.is_developer(member_snowflake=int(member_snowflake)):
            return "Developer"
    except NotDeveloper as e:
        bot.logger.warning(str(e).capitalize())
    try:
        if await guild_owner_service.is_guild_owner(
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
        ):
            return "Guild Owner"
    except NotGuildOwner as e:
        bot.logger.warning(str(e).capitalize())
    try:
        if await administrator_service.is_administrator(
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
        ):
            return "Administrator"
    except NotAdministrator as e:
        bot.logger.warning(str(e).capitalize())
    if channel_snowflake:
        try:
            if await coordinator_service.is_coordinator(
                channel_snowflake=int(channel_snowflake),
                guild_snowflake=int(guild_snowflake),
                member_snowflake=int(member_snowflake),
            ):
                return "Coordinator"
        except NotCoordinator as e:
            bot.logger.warning(str(e).capitalize())
        try:
            if await is_moderator(
                channel_snowflake=int(channel_snowflake),
                guild_snowflake=int(guild_snowflake),
                member_snowflake=int(member_snowflake),
            ):
                return "Moderator"
        except NotModerator as e:
            bot.logger.warning(str(e).capitalize())
    return "Everyone"


async def resolve_highest_role_at_all(
    member_snowflake: int,
) -> str:
    bot: DiscordBot = DiscordBot.get_instance()
    try:
        if await sysadmin_service.is_sysadmin(member_snowflake=int(member_snowflake)):
            return "Sysadmin"
    except sysadmin_service.NotSysadmin as e:
        bot.logger.warning(str(e).capitalize())
    try:
        if await developer_service.is_developer(member_snowflake=int(member_snowflake)):
            return "Developer"
    except NotDeveloper as e:
        bot.logger.warning(str(e).capitalize())
    try:
        if await guild_owner_service.is_guild_owner_at_all(
            member_snowflake=int(member_snowflake),
        ):
            return "Guild Owner"
    except NotGuildOwner as e:
        bot.logger.warning(str(e).capitalize())
    try:
        if await administrator_service.is_administrator_at_all(
            member_snowflake=int(member_snowflake),
        ):
            return "Administrator"
    except NotAdministrator as e:
        bot.logger.warning(str(e).capitalize())
    try:
        if await coordinator_service.is_coordinator_at_all(
            member_snowflake=int(member_snowflake),
        ):
            return "Coordinator"
    except NotCoordinator as e:
        bot.logger.warning(str(e).capitalize())
    try:
        if await is_moderator_at_all(
            member_snowflake=int(member_snowflake),
        ):
            return "Moderator"
    except NotModerator as e:
        bot.logger.warning(str(e).capitalize())
    return "Everyone"


def compare_ranks(sender_rank, target_rank) -> bool:
    try:
        if sender_rank <= target_rank:
            raise HasEqualOrLowerRole(PERMISSION_TYPES[target_rank])
    except HasEqualOrLowerRole as e:
        raise e
    return True


# async def can_list(
#     source=Union[commands.Context, discord.Interaction, discord.Message]
# ) -> tuple[list[discord.], list[int]]:
#     bot: DiscordBot = DiscordBot.get_instance()
#     available_channels = {}
#     available_guilds = {}
#     member_snowflake = resolve_author(source=source).id
#     verifications = (
#         ("all", sysadmin_service.is_sysadmin),
#         ("all", developer_service.is_developer),
#         ("guild", guild_owner_service.is_guild_owner),
#         ("guild", administrator_service.is_administrator),
#         ("channel", coordinator_service.is_coordinator),
#         ("channel", is_moderator),
#     )
#     for role_scope, verify in verifications:
#         if role_scope == "all":
#             try:
#                 if await verify(member_snowflake=int(member_snowflake)):
#                     available_guilds["all"] = bot.guilds
#                     available_channels["all"] = []
#                     for guild in bot.guilds:
#                         available_guilds[guild.id] = guild
#                         available_channels.setdefault(guild.id, [])
#                         for channel in guild.channels:
#                             if isinstance(channel, discord.VoiceChannel):
#                                 available_channels[guild.id].append(channel)
#                                 available_channels["all"].append(channel)
#             except commands.CheckFailure:
#                 pass
#         elif role_scope == "guild":
#             try:
#                 for guild in bot.guilds:
#                     if await verify(
#                         guild_snowflake=int(guild.id),
#                         member_snowflake=int(member_snowflake),
#                     ):
#                         available_guilds[guild.id] = guild
#                         available_channels.setdefault(guild.id, [])
#                         for channel in guild.channels:
#                             if isinstance(channel, discord.VoiceChannel):
#                                 available_channels[guild.id].append(channel)
#             except commands.CheckFailure:
#                 pass
#         elif role_scope == "channel":
#             try:
#                 for guild in bot.guilds:
#                     for channel in guild.channels:
#                         if await verify(
#                             channel_snowflake=int(channel.id),
#                             guild_snowflake=int(guild.id),
#                             member_snowflake=int(member_snowflake),
#                         ):
#                             available_guilds[guild.id] = guild
#                             available_channels.setdefault(guild.id, [])
#                             if isinstance(channel, discord.VoiceChannel):
#                                 available_channels[guild.id].append(channel)
#             except commands.CheckFailure:
#                 pass
#     for gid in list(available_channels):
#         available_channels[gid] = list(
#             {c.id: c for c in available_channels[gid]}.values()
#         )
#     return available_channels, available_guilds
