"""!/bin/python3
guild_event_listeners.py A discord.py cog containing guild event listeners for the Vyrtuous bot.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.utils.users import autoassign_role_service


class GuildEventListeners(commands.Cog):
    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    @commands.Cog.listener()
    async def on_member_update(
        self, before: discord.Member, after: discord.Member
    ) -> None:
        bot: DiscordBot = DiscordBot.get_instance()
        if before.roles == after.roles:
            return
        guild_snowflake = before.guild.id
        autoassign_role_snowflakes = await autoassign_role_service.get_autoassign_roles(
            guild_snowflake
        )
        before_role_snowflakes = {str(r.id) for r in before.roles}
        after_role_snowflakes = {str(r.id) for r in after.roles}
        added_roles = after_role_snowflakes - before_role_snowflakes
        removed_roles = before_role_snowflakes - after_role_snowflakes
        invincible_member_set = bot.registry.get(MemberState).invincible.get(
            after.guild.id, None
        )
        if added_roles:
            for added_role in added_roles:
                if int(added_role) in autoassign_role_snowflakes:
                    await autoassign_role_service.added_role(
                        guild_snowflake=guild_snowflake,
                        member_snowflake=before.id,
                        role_snowflake=int(added_role),
                    )
                    self.__bot.logger.debug(f"Added roles: {', '.join(added_roles)}")
            if invincible_member_set and after.id in invincible_member_set:
                try:
                    roles = [
                        role
                        for role_snowflake in added_roles
                        if (role := after.guild.get_role(int(role_snowflake)))
                        is not None
                    ]
                    await after.remove_roles(
                        *roles, reason="Hero restricts role addition"
                    )
                except discord.HTTPException:
                    bot.logger.debug(
                        f"Unable to remove roles added to hero {after.display_name}."
                    )
        elif removed_roles:
            for removed_role in removed_roles:
                await autoassign_role_service.removed_role(
                    guild_snowflake=guild_snowflake,
                    member_snowflake=before.id,
                    role_snowflake=int(removed_role),
                )
                self.__bot.logger.debug(f"Removed roles: {', '.join(removed_roles)}")

    @commands.Cog.listener()
    async def on_guild_role_delete(self, role: discord.Role) -> None:
        guild_snowflake = role.guild.id
        role_snowflake = role.id
        for member in role.members:
            await autoassign_role_service.removed_role(
                guild_snowflake=guild_snowflake,
                member_snowflake=member.id,
                role_snowflake=role_snowflake,
            )
            self.__bot.logger.debug(
                f"Removed role ({role.id}) from server ({role.guild.name})."
            )


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildEventListeners(bot=bot))
