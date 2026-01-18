"""check_service.py The purpose of this program is to provide the check_service module.

Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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

from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.roles.administrator import Administrator
from vyrtuous.database.roles.coordinator import Coordinator
from vyrtuous.database.roles.developer import Developer
from vyrtuous.database.roles.moderator import Moderator
from vyrtuous.service.logging_service import logger
from vyrtuous.utils.permission import PERMISSION_TYPES


class NotCoordinator(commands.CheckFailure):
    def __init__(
        self, message="You are not a coordinator in this channel and cannot do this."
    ):
        super().__init__(message)


class NotAdministrator(commands.CheckFailure):
    def __init__(self, message="You are not an administrator and cannot do this."):
        super().__init__(message)


class NotModerator(commands.CheckFailure):
    def __init__(
        self, message="You are not a moderator in this channel and cannot do this."
    ):
        super().__init__(message)


class NotGuildOwner(commands.CheckFailure):
    def __init__(self, message="You are not the guild owner and cannot do this."):
        super().__init__(message)


class NotSystemOwner(commands.CheckFailure):
    def __init__(self, message="You are not the system owner and cannot do this."):
        super().__init__(message)


class NotDeveloper(commands.CheckFailure):
    def __init__(self, message="You are not a developer and cannot do this."):
        super().__init__(message)


class NotAtHome(commands.CheckFailure):
    def __init__(self, message="You are not in the home server and cannot do this."):
        super().__init__(message)


def at_home(
    source: Union[commands.Context, discord.Interaction, discord.Message],
) -> bool:
    bot = DiscordBot.get_instance()
    if source.guild.id == int(bot.config["discord_testing_guild_snowflake"]):
        return True
    raise NotAtHome()


async def is_moderator(
    source: Union[commands.Context, discord.Interaction, discord.Message],
    member_snowflake: int,
):
    if isinstance(source, discord.Interaction):
        member_snowflake = source.user.id
    elif isinstance(source, commands.Context):
        member_snowflake = source.author.id
    elif isinstance(source, discord.Message):
        member_snowflake = source.author.id
    else:
        member_snowflake = None
    moderator = await Moderator.select(
        channel_snowflake=source.channel.id,
        guild_snowflake=source.guild.id,
        member_snowflake=member_snowflake,
    )
    if not moderator:
        raise NotModerator()
    return True


async def is_coordinator(
    source: Union[commands.Context, discord.Interaction, discord.Message],
    member_snowflake: int,
):
    if isinstance(source, discord.Interaction):
        member_snowflake = source.user.id
    elif isinstance(source, commands.Context):
        member_snowflake = source.author.id
    elif isinstance(source, discord.Message):
        member_snowflake = source.author.id
    else:
        member_snowflake = None
    coordinator = await Coordinator.select(
        channel_snowflake=source.channel.id,
        guild_snowflake=source.guild.id,
        member_snowflake=member_snowflake,
    )
    if not coordinator:
        raise NotCoordinator()
    return True


async def is_administrator(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    if isinstance(source, discord.Interaction):
        member_snowflake = source.user.id
    elif isinstance(source, commands.Context):
        member_snowflake = source.author.id
    elif isinstance(source, discord.Message):
        member_snowflake = source.author.id
    else:
        member_snowflake = None
    administrator = await Administrator.select(
        guild_snowflake=source.guild.id,
        member_snowflake=member_snowflake,
    )
    if not administrator:
        raise NotAdministrator()
    return True


async def is_developer(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    if isinstance(source, discord.Interaction):
        member_snowflake = source.user.id
    elif isinstance(source, commands.Context):
        member_snowflake = source.author.id
    elif isinstance(source, discord.Message):
        member_snowflake = source.author.id
    else:
        member_snowflake = None
    developer = await Developer.select(
        member_snowflake=member_snowflake,
    )
    if developer:
        return True
    else:
        raise NotDeveloper()


async def is_guild_owner(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    if isinstance(source, discord.Interaction):
        member_snowflake = source.user.id
    elif isinstance(source, commands.Context):
        member_snowflake = source.author.id
    elif isinstance(source, discord.Message):
        member_snowflake = source.author.id
    else:
        member_snowflake = None
    if source.guild.owner_id != member_snowflake:
        raise NotGuildOwner()
    return True


async def is_system_owner(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    bot = DiscordBot.get_instance()
    if isinstance(source, discord.Interaction):
        user = source.user
    elif isinstance(source, commands.Context):
        user = source.author
    elif isinstance(source, discord.Message):
        user = source.author
    system_owner_id = int(bot.config["discord_owner_id"])
    if system_owner_id != user.id:
        raise NotSystemOwner()
    return True


def developer_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (is_system_owner, is_developer):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure("You are not a system owner or developer.")

    predicate._permission_level = "Developer"
    return commands.check(predicate)


def guild_owner_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (is_system_owner, is_developer, is_guild_owner):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure(
            "You are not a system owner, developer or guild owner in this server."
        )

    predicate._permission_level = "Guild Owner"
    return commands.check(predicate)


def sys_owner_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        if await is_system_owner(source):
            return True
        raise commands.CheckFailure("You are not a system owner.")

    predicate._permission_level = "System Owner"
    return commands.check(predicate)


def administrator_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (is_system_owner, is_developer, is_guild_owner, is_administrator):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure(
            "You are not a system owner, developer, guild owner or administrator in this server."
        )

    predicate._permission_level = "Administrator"
    return commands.check(predicate)


def coordinator_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (
            is_system_owner,
            is_developer,
            is_guild_owner,
            is_administrator,
            is_coordinator,
        ):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure(
            "You are not a system owner, developer, guild owner, administrator or coordinator in this channel."
        )

    predicate._permission_level = "Coordinator"
    return commands.check(predicate)


def moderator_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (
            is_system_owner,
            is_developer,
            is_guild_owner,
            is_administrator,
            is_coordinator,
            is_moderator,
        ):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure(
            "You are not a system owner, developer, guild owner, administrator, coordinator or moderator in this channel."
        )

    predicate._permission_level = "Moderator"
    return commands.check(predicate)


async def check(
    source: Union[commands.Context, discord.Interaction, discord.Message],
    lowest_role: str = None,
    member_snowflake: int = None,
) -> str:
    verifications = (
        ("System Owner", is_system_owner),
        ("Developer", is_developer),
        ("Guild Owner", is_guild_owner),
        ("Administrator", is_administrator),
        ("Coordinator", lambda s: is_coordinator(s, member_snowflake)),
        ("Moderator", lambda s: is_moderator(s, member_snowflake)),
    )
    passed_lowest = False
    for role_name, verify in verifications:
        try:
            if await verify(source):
                return role_name
        except commands.CheckFailure:
            if lowest_role is not None and passed_lowest:
                raise
        if role_name == lowest_role:
            passed_lowest = True
    return "Everyone"


async def member_is_guild_owner(guild_snowflake: int, member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild and guild.owner_id == member_snowflake:
        return True
    raise NotGuildOwner


async def member_is_system_owner(member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    if int(bot.config["discord_owner_id"]) == member_snowflake:
        return True
    raise NotSystemOwner


async def member_is_developer(member_snowflake: int) -> bool:
    developer = await Developer.select(member_snowflake=member_snowflake)
    if not developer:
        raise NotDeveloper
    return True


async def member_is_administrator(guild_snowflake: int, member_snowflake: int) -> bool:
    administrator = await Administrator.select(
        guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
    )
    if not administrator:
        raise NotAdministrator
    return True


async def member_is_coordinator(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> bool:
    coordinator = await Coordinator.select(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
    )
    if not coordinator:
        raise NotCoordinator
    return True


async def member_is_moderator(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> bool:
    moderator = await Moderator.select(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
    )
    if not moderator:
        raise NotModerator
    return True


async def resolve_highest_permission_role(
    member_snowflake: int, channel_snowflake=None, guild_snowflake=None
):
    try:
        if await member_is_system_owner(member_snowflake=member_snowflake):
            return PERMISSION_TYPES.index("System Owner")
    except NotSystemOwner as e:
        logger.warning(str(e).capitalize())
    try:
        if await member_is_developer(member_snowflake=member_snowflake):
            return PERMISSION_TYPES.index("Developer")
    except NotDeveloper as e:
        logger.warning(str(e).capitalize())
    try:
        if await member_is_guild_owner(
            guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
        ):
            return PERMISSION_TYPES.index("Guild Owner")
    except NotGuildOwner as e:
        logger.warning(str(e).capitalize())
    try:
        if await member_is_administrator(
            guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
        ):
            return PERMISSION_TYPES.index("Administrator")
    except NotAdministrator as e:
        logger.warning(str(e).capitalize())
    if channel_snowflake:
        try:
            if await member_is_coordinator(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                member_snowflake=member_snowflake,
            ):
                return PERMISSION_TYPES.index("Coordinator")
        except NotCoordinator as e:
            logger.warning(str(e).capitalize())
        try:
            if await member_is_moderator(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                member_snowflake=member_snowflake,
            ):
                return PERMISSION_TYPES.index("Moderator")
        except NotModerator as e:
            logger.warning(str(e).capitalize())
    return PERMISSION_TYPES.index("Everyone")
