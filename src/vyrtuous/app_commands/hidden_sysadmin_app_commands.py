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

from vyrtuous.app_commands.help_app_command import skip_app_command_help_discovery
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.upload import upload_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import moderator_service


class HiddenSysadminAppCommands(commands.Cog):
    PERMISSION_LEVEL = "Sysadmin"

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def interaction_check(self, interaction: discord.Interaction):
        if interaction.guild is None:
            raise commands.CheckFailure("This command must be executed inside a server.")
        if interaction.channel is None:
            raise commands.CheckFailure("This command must be executed in a valid channel.")
        await moderator_service.check_minimum_role(
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @app_commands.command(name="upload", description="Create the upload document.")
    @skip_app_command_help_discovery()
    async def uploads_app_command(self, interaction: discord.Interaction) -> None:
        tick = Tick(bot=self.__bot, interaction=interaction)
        return await upload_service.build_latex_document()


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenSysadminAppCommands(bot=bot))
