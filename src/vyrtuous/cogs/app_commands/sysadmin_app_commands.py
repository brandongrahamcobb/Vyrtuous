"""sysadmin_commands.py A discord.py cog containing sysadmin commands for the Vyrtuous bot.

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

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.bug import Bug
from vyrtuous.db.roles.developer import Developer
from vyrtuous.db.roles.sysadmin import sysadmin_predicator
from vyrtuous.fields.snowflake import (
    AppMemberSnowflake,
)
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import DiscordObject


class SysadminAppCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    @app_commands.command(name="assign", description="Assign developer.")
    @app_commands.describe(
        reference="Include an issue reference ID",
        member="Tag a member or include their ID",
    )
    @sysadmin_predicator()
    async def assign_bug_to_developer_app_command(
        self,
        interaction: discord.Interaction,
        reference: str,
        member: AppMemberSnowflake,
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        embed = await Bug.assign_bug_to_developer(
            reference=reference, member_dict=member_dict
        )
        return await state.end(success=embed)

    @app_commands.command(name="dev", description="Grant/revoke devs.")
    @app_commands.describe(member="Tag a member or include their ID")
    @sysadmin_predicator()
    async def toggle_developer_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = StateService(interaction=interaction)
        snowflake_kwargs = {
            "channel_snowflake": interaction.channel.id,
            "guild_snowflake": interaction.guild.id,
            "member_snowflake": interaction.user.id,
        }
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        msg = await Developer.toggle_developer(
            member_dict=member_dict, snowflake_kwargs=snowflake_kwargs
        )
        return await state.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(SysadminAppCommands(bot))
