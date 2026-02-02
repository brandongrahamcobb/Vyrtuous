"""guild_owner_commands.py A discord.py cog containing guild owner commands for the Vyrtuous bot.

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
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cogs.help_command import skip_help_discovery
from vyrtuous.fields.snowflake import AppMemberSnowflake, AppRoleSnowflake
from vyrtuous.service.discord_object_service import DiscordObject
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.roles.administrator_service import AdministratorRoleService
from vyrtuous.service.roles.guild_owner_service import guild_owner_predicator
from vyrtuous.service.state_service import StateService
from vyrtuous.utils.invincibility import Invincibility


class GuildOwnerAppCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    @app_commands.command(name="admin", description="Role -> Administrator.")
    @guild_owner_predicator()
    async def toggle_administrator_by_role_app_command(
        self, interaction: discord.Interaction, role: AppRoleSnowflake
    ):
        state = StateService(interaction=interaction)
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        do = DiscordObject(interaction=interaction)
        role_dict = await do.determine_from_target(target=role)
        pages = await AdministratorRoleService.toggle_administrator_role(
            role_dict=role_dict, snowflake_kwargs=snowflake_kwargs
        )
        await StateService.send_pages(
            title="Administrator Role", pages=pages, state=state
        )

    @app_commands.command(name="hero", description="Grant/revoke invincibility.")
    @app_commands.describe(member="Tag a member or include their ID")
    @guild_owner_predicator()
    @skip_help_discovery()
    async def invincibility_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        where_kwargs = member_dict.get("columns", None)
        enabled = Invincibility.toggle_enabled()
        if enabled:
            Invincibility.add_invincible_member(**where_kwargs)
            await Invincibility.unrestrict(**where_kwargs)
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member_dict.get("mention", None)}."
            )
        else:
            Invincibility.remove_invincible_member(**where_kwargs)
            msg = f"Invincibility has been disabled for {member_dict.get("mention", None)}"
        return await state.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildOwnerAppCommands(bot))
