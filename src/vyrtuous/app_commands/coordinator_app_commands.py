"""!/bin/python3

coordinator_app_commands.py A discord.py cog containing coordinator commands for the Vyrtuous bot.

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
from vyrtuous.models.duration import AppDuration, DurationObject, DurationWrapper
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.utils.channels import automute_channel_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import moderator_service


class CoordinatorAppCommands(commands.Cog):

    PERMISSION_LEVEL = "Coordinator"

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

    @app_commands.command(name="mod", description="Grant/revoke mods.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
        channel="Specify a channel ID/mention.",
    )
    async def toggle_moderator_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        await interaction.response.defer()
        tick = Tick(bot=self.__bot, interaction=interaction)
        if channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = interaction.channel.id
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        await moderator_service.has_equal_or_lower_role(
            target_member_snowflake=member_snowflake,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
        )
        message = await interaction.original_response()
        message_snowflake = message.id
        msg = await moderator_service.toggle_moderator(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake,
            message_channel_snowflake=message.channel.id,
        )
        return await tick.end(success=msg)

    @app_commands.command(name="roles", description="List role members.")
    @app_commands.describe(
        role="Specify a role ID/mention.", guild="Specify a server ID."
    )
    async def list_roles_app_command(
        self,
        interaction: discord.Interaction,
        role: app_commands.Transform[TargetObject, AppTarget],
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            members = interaction.guild.members
        else:
            if isinstance(guild.target, discord.Guild):
                members = guild.target.members
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if isinstance(role.target, discord.Role):
            role_name = role.target.name
            color = (
                role.target.color
                if role.target.color.value
                else discord.Color.blurple()
            )
        else:
            return await tick.end(warning="This command must target a valid role.")
        embeds = []
        members = [member for member in members if role in member.roles]
        chunk_size = 12
        for index in range(0, len(members), chunk_size):
            chunk = members[index : index + chunk_size]
            description = "\n".join(
                f"{position + 1}. {member.mention} ({member.id})"
                for position, member in enumerate(chunk, start=index)
            )
            embed = discord.Embed(
                title=f"{role_name} Members",
                description=description or "No members found.",
                color=color,
            )
            embeds.append(embed)
        if not embeds:
            embed = discord.Embed(
                title=f"{role_name} Members",
                description="No members found.",
                color=color,
            )
            embeds.append(embed)
        return await tick.end(success=embeds)

    @app_commands.command(name="automute", description="Start/stop automute")
    @app_commands.describe(
        channel="Specify a channel ID/mention.",
        duration="Specify duration as m/h/d.",
    )
    async def toggle_automute_app_command(
        self,
        interaction: discord.Interaction,
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
        duration: app_commands.Transform[DurationWrapper | None, AppDuration] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = interaction.channel.id
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        if duration is None:
            duration_obj = DurationObject(number=2, prefix="", sign=1, unit="h")
        else:
            duration_obj = duration.duration

        pages = await automute_channel_service.toggle_automute(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            duration=duration_obj,
        )
        return await tick.end(success=pages)

    async def cog_app_command_error(self, interaction, error):
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorAppCommands(bot=bot))
