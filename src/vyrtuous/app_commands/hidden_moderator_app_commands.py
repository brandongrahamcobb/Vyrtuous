"""!/bin/python3
moderator_app_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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
from vyrtuous.inc.helpers import at_home
from vyrtuous.listing import list_administrators, list_vegans
from vyrtuous.models import multi_converter
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import moderator_service, vegan_service

# from vyrtuous.view.data_view import DataView
# from vyrtuous.view.modify_infraction_view import ModifyInfractionView


class HiddenModeratorAppCommands(commands.Cog):

    PERMISSION_LEVEL = "Moderator"

    def __init__(
        self,
        *,
        bot: DiscordBot,
    ):
        self.__bot = bot

    async def interaction_check(self, interaction: discord.Interaction):
        if interaction.guild is None:
            raise commands.CheckFailure("This command must be used inside a server.")
        if interaction.channel is None:
            raise commands.CheckFailure("This command must be used in a valid channel.")
        await moderator_service.check_minimum_role(
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    # @app_commands.command(name="data", description="Create a chart.")
    # async def create_data_app_command(self, interaction: discord.Interaction):
    #     tick = Tick(bot=self.__bot, interaction=interaction)
    #     view = DataView(
    #         bot=self.__bot,
    #         ban_service=self.__ban_service,
    #         duration_builder=self.__duration_builder,
    #         flag_service=self.__flag_service,
    #         interaction=interaction,
    #         moderator_service=self.__moderator_service,
    #         text_mute_service=self.__text_mute_service,
    #         voice_mute_service=self.__voice_mute_service,
    #         tick=tick,
    #     )
    #     await view.setup()
    #     await interaction.response.send_message(
    #         content="Select a channel, duration and infraction",
    #         view=view,
    #         ephemeral=True,
    #     )
    #
    @app_commands.command(name="admins", description="Lists admins.")
    @app_commands.describe(
        target="Specify one of: `all`, a member ID/mention or server ID.",
    )
    @skip_app_command_help_discovery()
    async def list_administrators_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_administrators.build_pages(is_at_home=is_at_home, obj=obj)
        return await tick.end(success=pages)

    @app_commands.command(name="survey", description="Survey stage members.")
    @app_commands.describe(
        channel="Specify channel ID/mention.",
    )
    @skip_app_command_help_discovery()
    async def stage_survey_app_command(
        self, interaction: discord.Interaction, channel: discord.abc.GuildChannel | None
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        obj = channel or interaction.channel
        pages = await moderator_service.survey(
            channel=obj,
        )
        return await tick.end(success=pages)

    @app_commands.command(name="vcow", description="Toggle vegan.")
    @app_commands.describe(member="Specify member ID/mention.", notes="Include notes.")
    async def toggle_vegan_app_command(
        self,
        interaction: discord.Interaction,
        member: str,
        notes: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server")
        if not vegan_service.is_vegan(
            guild_snowflake=interaction.guild.id, member_snowflake=interaction.user.id
        ):
            return await tick.end(warning="Author is not a vegan.")
        target = multi_converter.transform(interaction=interaction, argument=member)
        if isinstance(target, int):
            member_snowflake = target
        elif isinstance(target, discord.Member):
            member_snowflake = target.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        embed = await vegan_service.toggle_vegan(
            guild_snowflake=interaction.guild.id,
            member_snowflake=member_snowflake,
            notes=notes,
        )
        return await tick.end(success=embed)

    @app_commands.command(name="vegans", description="List new vegans.")
    @app_commands.describe(
        target="Specify one of: `all`, member ID/mention or server ID.",
    )
    @skip_app_command_help_discovery()
    async def list_new_vegans_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.guild
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_vegans.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenModeratorAppCommands(bot=bot))
