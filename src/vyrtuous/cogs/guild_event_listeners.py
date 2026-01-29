"""event_listeners.py A discord.py cog containing event listeners for the Vyrtuous bot.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.actions.ban import Ban
from vyrtuous.db.actions.role import Role
from vyrtuous.db.actions.text_mute import TextMute
from vyrtuous.db.roles.administrator import AdministratorRole
from vyrtuous.db.roles.guild_owner import GuildOwner

from vyrtuous.utils.logger import logger


class GuildEventListeners(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot

    @commands.Cog.listener()
    async def on_guild_update(self, before: discord.Guild, after: discord.Guild):
        if before.owner_id != after.owner_id:
            where_kwargs = {
                "guild_snowflake": before.guild.id,
                "member_snowflake": before.owner_id,
            }
            set_kwargs = {
                "guild_snowflake": after.guild.id,
                "member_snowflake": after.owner_id,
            }
            await GuildOwner.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)

    @commands.Cog.listener()
    async def on_member_update(self, before: discord.Member, after: discord.Member):
        if before.roles == after.roles:
            return
        guild_snowflake = before.guild.id
        before_role_snowflakes = {r.id for r in before.roles}
        after_role_snowflakes = {r.id for r in after.roles}
        added_roles = after_role_snowflakes - before_role_snowflakes
        removed_roles = before_role_snowflakes - after_role_snowflakes
        if added_roles:
            await AdministratorRole.added_role(
                guild_snowflake=guild_snowflake,
                member_snowflake=before.id,
                role_snowflake=added_roles[0],
            )
            logger.info(f"Added roles: {', '.join(added_roles)}")
        elif removed_roles:
            await AdministratorRole.removed_role(
                guild_snowflake=guild_snowflake,
                member_snowflake=before.id,
                role_snowflake=removed_roles[0],
            )
            logger.info(f"Removed roles: {', '.join(removed_roles)}")

    @commands.Cog.listener()
    async def on_guild_role_delete(self, role: discord.Role):
        guild_snowflake = role.guild.id
        for member in role.members:
            await AdministratorRole.removed_role(
                guild_snowflake=guild_snowflake,
                member_snowflake=member.id,
                role_snowflake=role.id,
            )
            logger.info(f"Removed role ({role.id}) from server ({role.guild.name}).")


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildEventListeners(bot))
