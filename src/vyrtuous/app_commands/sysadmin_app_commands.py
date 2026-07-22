"""!/bin/python3
sysadmin_app_commands.py A discord.py cog containing sysadmin commands for the Vyrtuous bot.

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
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import developer_service, moderator_service


class SysadminAppCommands(commands.Cog):
    PERMISSION_LEVEL = "Sysadmin"

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def interaction_check(self, interaction: discord.Interaction):
        if interaction.guild is None:
            raise commands.CheckFailure(
                "This command must be executed inside a server."
            )
        if interaction.channel is None:
            raise commands.CheckFailure(
                "This command must be executed in a valid channel."
            )
        await moderator_service.check_minimum_role(
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @app_commands.command(name="dev", description="Grant/revoke devs.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def toggle_developer_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
    ) -> discord.Message:
        await interaction.response.defer()
        tick = Tick(bot=self.__bot, interaction=interaction)
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        message = await interaction.original_response()
        message_snowflake = message.id
        msg = await developer_service.toggle_developer(
            author_snowflake=interaction.user.id,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake,
            message_channel_snowflake=message.channel.id,
        )
        return await tick.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(SysadminAppCommands(bot=bot))
