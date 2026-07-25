"""!/bin/python3
hidden_moderator_app_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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

from vyrtuous.app_commands.help_app_command import skip_app_command_help_discovery
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.listing import list_administrators, list_vegans
from vyrtuous.models.target import AppTarget, TargetObject
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
        target="Specify one of: a member ID/mention, role ID/mention or server ID.",
    )
    @skip_app_command_help_discovery()
    async def list_administrators_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None:
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if interaction.guild.id != guild_snowflake:
            await moderator_service.check_minimum_role(
                member_snowflake=interaction.user.id,
                lowest_role="Developer",
            )
        if target is None:
            obj = interaction.guild
        else:
            obj = target.target
        pages = await list_administrators.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @app_commands.command(name="survey", description="Survey stage members.")
    @app_commands.describe(
        channel="Specify channel ID/mention.",
    )
    @skip_app_command_help_discovery()
    async def survey_app_command(
        self,
        interaction: discord.Interaction,
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if channel is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        pages = await moderator_service.survey(
            channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake
        )
        return await tick.end(success=pages)

    @app_commands.command(name="vcow", description="Toggle vegan.")
    @app_commands.describe(member="Specify member ID/mention.", notes="Include notes.")
    async def toggle_vegan_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        notes: str | None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if not vegan_service.is_vegan(
            guild_snowflake=guild_snowflake, member_snowflake=interaction.user.id
        ):
            return await tick.end(warning="Author is not a vegan.")
        if isinstance(member, int):
            member_snowflake = member
        elif isinstance(member, discord.Member):
            member_snowflake = member.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        embed = await vegan_service.toggle_vegan(
            guild_snowflake=guild_snowflake,
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
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        pages = await list_vegans.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    async def cog_app_command_error(self, interaction, error):
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenModeratorAppCommands(bot=bot))
