"""guild_owner.py The purpose of this program is to inherit from the DatabaseFactory to provide the guild owner role.

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

from vyrtuous.base.service import Service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.author import resolve_author
from vyrtuous.db.roles.dev.developer_service import is_developer_wrapper
from vyrtuous.db.roles.sysadmin.sysadmin_service import is_sysadmin_wrapper


class NotGuildOwner(commands.CheckFailure):
    def __init__(
        self,
        message="Member is not a sysadmin, developer or guild owner in this server.",
    ):
        super().__init__(message)


async def is_guild_owner_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    return await is_guild_owner(
        guild_snowflake=source.guild.id, member_snowflake=int(member_snowflake)
    )


def guild_owner_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (
            is_sysadmin_wrapper,
            is_developer_wrapper,
            is_guild_owner_wrapper,
        ):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure(
            "Member is not a sysadmin, developer or guild owner in this server."
        )

    predicate._permission_level = "Guild Owner"
    return commands.check(predicate)


async def is_guild_owner(guild_snowflake: int, member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild and guild.owner_id == member_snowflake:
        return True
    raise NotGuildOwner


class GuildOwnerService(Service):

    pass
