"""!/bin/python3
guild_event_listeners.py A discord.py cog containing guild event listeners for the Vyrtuous bot.

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

import discord
from discord.ext import commands

from vyrtuous.administrator.administrator_service import AdministratorRoleService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.owner.guild_owner import GuildOwner
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.logger import logger


class GuildEventListeners(commands.Cog):
    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__database_factory.model = GuildOwner
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__administrator_role_service = AdministratorRoleService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )

    @commands.Cog.listener()
    async def on_guild_update(self, before: discord.Guild, after: discord.Guild):
        if before.owner_id != after.owner_id:
            where_kwargs = {
                "guild_snowflake": int(before.id),
                "member_snowflake": before.owner_id,
            }
            set_kwargs = {
                "guild_snowflake": int(after.id),
                "member_snowflake": after.owner_id,
            }
            await self.__database_factory.update(
                set_kwargs=set_kwargs, where_kwargs=where_kwargs
            )

    @commands.Cog.listener()
    async def on_member_update(self, before: discord.Member, after: discord.Member):
        if before.roles == after.roles:
            return
        guild_snowflake = before.guild.id
        before_role_snowflakes = {str(r.id) for r in before.roles}
        after_role_snowflakes = {str(r.id) for r in after.roles}
        added_roles = after_role_snowflakes - before_role_snowflakes
        removed_roles = before_role_snowflakes - after_role_snowflakes
        kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "member_snowflake": int(before.id),
        }
        if added_roles:
            for added_role in added_roles:
                kwargs.update({"role_snowflake": int(added_role)})
                await self.__administrator_role_service.added_role(kwargs=kwargs)
                logger.info(f"Added roles: {', '.join(added_roles)}")
        elif removed_roles:
            for removed_role in removed_roles:
                kwargs.update({"role_snowflake": int(removed_role)})
            await self.__administrator_role_service.removed_role(kwargs=kwargs)
            logger.info(f"Removed roles: {', '.join(removed_roles)}")

    @commands.Cog.listener()
    async def on_guild_role_delete(self, role: discord.Role):
        guild_snowflake = role.guild.id
        kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "role_snowflake": str(role.id),
        }
        for member in role.members:
            await self.__administrator_role_service.removed_role(kwargs=kwargs)
            logger.info(f"Removed role ({role.id}) from server ({role.guild.name}).")


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildEventListeners(bot=bot))
