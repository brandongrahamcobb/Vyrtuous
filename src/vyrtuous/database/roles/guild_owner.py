"""check_service.py The purpose of this program is to provide the check_service module.

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

from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.roles.developer import is_developer_wrapper
from vyrtuous.database.roles.sysadmin import is_sysadmin_wrapper
from vyrtuous.service.member_snowflake import get_member_snowflake

class NotGuildOwner(commands.CheckFailure):
    def __init__(self, message="You are not the guild owner and cannot do this."):
        super().__init__(message)


async def is_guild_owner_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member_snowflake = get_member_snowflake(source=source)
    return is_guild_owner(guild_snowflake=source.guild.id, member_snowflake=member_snowflake)

def guild_owner_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (is_sysadmin_wrapper, is_developer_wrapper, is_guild_owner_wrapper):
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

async def is_guild_owner(guild_snowflake: int, member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild and guild.owner_id == member_snowflake:
        return True
    raise NotGuildOwner