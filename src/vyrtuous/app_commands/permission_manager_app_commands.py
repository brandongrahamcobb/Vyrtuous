"""!/bin/python3
permission_manager_app_commands.py A discord.py cog containing permission manager commands for the Vyrtuous bot.

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
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.permissions import permission_service
from vyrtuous.view.grant_view import GrantView
from vyrtuous.view.view_context import ViewContext


class PermissionManagerAppCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    @app_commands.command(name="grant", description="Grant permission levels.")
    @app_commands.describe(
        member="Specify a member ID/mention.", guild="Specify a server ID."
    )
    async def grant_group_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = interaction.channel.id
        ctx = ViewContext(
            interaction=interaction,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = GrantView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the group", view=view, ephemeral=True
        )

    async def cog_app_command_error(self, interaction, error):
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))


async def setup(bot: DiscordBot):
    await bot.add_cog(PermissionManagerAppCommands(bot=bot))
