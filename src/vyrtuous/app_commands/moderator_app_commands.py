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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.ban import Ban
from vyrtuous.db.moderator import NotAppModerator
from vyrtuous.db.text_mute import TextMute
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.modal.duration_modal import DurationModal
from vyrtuous.modal.reason_modal import ReasonModal
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import (
    administrator_service,
    coordinator_service,
    developer_service,
    guild_owner_service,
    moderator_service,
    sysadmin_service,
)

# from vyrtuous.view.data_view import DataView
# from vyrtuous.view.infraction_view import InfractionView
# from vyrtuous.view.modify_infraction_view import ModifyInfractionView
# from vyrtuous.view.view_context import ViewContext


class ModeratorAppCommands(commands.Cog):

    PERMISSION_LEVEL = "Moderator"

    def __init__(
        self,
        *,
        bot: DiscordBot,
    ):
        self.__bot = bot
        self.__duration_builder = DurationBuilder()

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
    # @app_commands.command(name="duration", description="Modify a duration.")
    # @app_commands.describe(member="The ID or mention of the member.")
    # async def change_moderation_duration_app_command(
    #     self, interaction: discord.Interaction, member: discord.Member
    # ):
    #     tick = Tick(bot=self.__bot, interaction=interaction)
    #     ctx = ViewContext(
    #         ban_service=self.__ban_service,
    #         flag_service=self.__flag_service,
    #         interaction=interaction,
    #         moderator_service=self.__moderator_service,
    #         text_mute_service=self.__text_mute_service,
    #         voice_mute_service=self.__voice_mute_service,
    #     )
    #     await ctx.setup(target_member_snowflake=member.id)
    #     view = ModifyInfractionView(
    #         ctx=ctx,
    #         modal=DurationModal,
    #         tick=tick,
    #     )
    #     await view.setup()
    #     await interaction.response.send_message(
    #         content="Select a channel and an infraction", view=view, ephemeral=True
    #     )
    #
    # @app_commands.command(name="reason", description="Modify a reason.")
    # @app_commands.describe(member="The ID or mention of the member.")
    # async def change_moderation_reason_app_command(
    #     self, interaction: discord.Interaction, member: discord.Member
    # ):
    #     tick = Tick(bot=self.__bot, interaction=interaction)
    #     ctx = ViewContext(
    #         ban_service=self.__ban_service,
    #         flag_service=self.__flag_service,
    #         interaction=interaction,
    #         moderator_service=self.__moderator_service,
    #         text_mute_service=self.__text_mute_service,
    #         voice_mute_service=self.__voice_mute_service,
    #     )
    #     await ctx.setup(target_member_snowflake=member.id)
    #     view = ModifyInfractionView(
    #         ctx=ctx,
    #         modal=ReasonModal,
    #         tick=tick,
    #         ban_service=self.__ban_service,
    #         flag_service=self.__flag_service,
    #         text_mute_service=self.__text_mute_service,
    #         voice_mute_service=self.__voice_mute_service,
    #     )
    #     await view.setup()
    #     await interaction.response.send_message(
    #         content="Select a channel and an infraction", view=view, ephemeral=True
    #     )
    #
    # @app_commands.command(name="vban", description="Create a ban.")
    # @app_commands.describe(member="The ID or mention of the member.")
    # async def create_ban_app_command(
    #     self, interaction: discord.Interaction, member: discord.Member
    # ):
    #     tick = Tick(bot=self.__bot, interaction=interaction)
    #     ctx = ViewContext(
    #         ban_service=self.__ban_service,
    #         flag_service=self.__flag_service,
    #         interaction=interaction,
    #         moderator_service=self.__moderator_service,
    #         text_mute_service=self.__text_mute_service,
    #         voice_mute_service=self.__voice_mute_service,
    #     )
    #     ctx.infraction = Ban
    #     await ctx.setup(target_member_snowflake=member.id)
    #     view = InfractionView(
    #         cap_service=self.__cap_service,
    #         ctx=ctx,
    #         duration_builder=self.__duration_builder,
    #         modal=ReasonModal,
    #         tick=tick,
    #     )
    #     await view.setup()
    #     await interaction.response.send_message(
    #         content="Select a channel and a duration", view=view, ephemeral=True
    #     )
    #
    # @app_commands.command(name="vmute", description="Create a mute.")
    # @app_commands.describe(member="The ID or mention of the member.")
    # async def create_voice_mute_app_command(
    #     self, interaction: discord.Interaction, member: discord.Member
    # ):
    #     tick = Tick(bot=self.__bot, interaction=interaction)
    #     ctx = ViewContext(
    #         ban_service=self.__ban_service,
    #         flag_service=self.__flag_service,
    #         interaction=interaction,
    #         moderator_service=self.__moderator_service,
    #         text_mute_service=self.__text_mute_service,
    #         voice_mute_service=self.__voice_mute_service,
    #     )
    #     ctx.infraction = VoiceMute
    #     await ctx.setup(target_member_snowflake=member.id)
    #     view = InfractionView(
    #         cap_service=self.__cap_service,
    #         ctx=ctx,
    #         duration_builder=self.__duration_builder,
    #         modal=ReasonModal,
    #         tick=tick,
    #     )
    #     await view.setup()
    #     await interaction.response.send_message(
    #         content="Select a channel and a duration", view=view, ephemeral=True
    #     )
    #
    # @app_commands.command(name="vtmute", description="Create a text-mute.")
    # @app_commands.describe(member="The ID or mention of the member.")
    # async def create_text_mute_app_command(
    #     self, interaction: discord.Interaction, member: discord.Member
    # ):
    #     tick = Tick(bot=self.__bot, interaction=interaction)
    #     ctx = ViewContext(
    #         ban_service=self.__ban_service,
    #         flag_service=self.__flag_service,
    #         interaction=interaction,
    #         moderator_service=self.__moderator_service,
    #         text_mute_service=self.__text_mute_service,
    #         voice_mute_service=self.__voice_mute_service,
    #     )
    #     ctx.infraction = TextMute
    #     await ctx.setup(target_member_snowflake=member.id)
    #     view = InfractionView(
    #         cap_service=self.__cap_service,
    #         ctx=ctx,
    #         duration_builder=self.__duration_builder,
    #         modal=ReasonModal,
    #         tick=tick,
    #     )
    #     await view.setup()
    #     await interaction.response.send_message(
    #         content="Select a channel and a duration", view=view, ephemeral=True
    #     )
    #


async def setup(bot: DiscordBot):
    await bot.add_cog(ModeratorAppCommands(bot=bot))
