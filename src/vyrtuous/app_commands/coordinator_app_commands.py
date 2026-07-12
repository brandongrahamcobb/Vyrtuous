"""!/bin/python3

coordinator_app_commands.py A discord.py cog containing coordinator commands for the Vyrtuous bot.

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
from vyrtuous.models import multi_converter
from vyrtuous.utils.channels import automute_channel_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import moderator_service


class CoordinatorAppCommands(commands.Cog):

    PERMISSION_LEVEL = "Coordinator"

    def __init__(self, bot: DiscordBot):
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

    @app_commands.command(name="mod", description="Grant/revoke mods.")
    @app_commands.describe(
        member="Specify a member ID/mention.", channel="Specify a channel ID/mention."
    )
    async def toggle_moderator_app_command(
        self,
        interaction: discord.Interaction,
        member: str,
        channel: discord.abc.GuildChannel,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        target = multi_converter.transform(interaction=interaction, argument=member)
        if isinstance(target, int):
            member_snowflake = target
        elif isinstance(target, discord.Member):
            member_snowflake = target.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        await moderator_service.has_equal_or_lower_role(
            target_member_snowflake=member_snowflake,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
        )
        message = await interaction.original_response()
        message_snowflake = message.id
        msg = await moderator_service.toggle_moderator(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake,
            message_channel_snowflake=message.channel.id,
        )
        return await tick.end(success=msg)

    @app_commands.command(name="roles", description="List role members.")
    @app_commands.describe(role="Specify a role ID/mention.")
    async def list_roles_app_command(
        self, interaction: discord.Interaction, role: discord.Role
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        embeds = []
        members = [
            member for member in interaction.guild.members if role in member.roles
        ]
        chunk_size = 12
        for index in range(0, len(members), chunk_size):
            chunk = members[index : index + chunk_size]
            description = "\n".join(
                f"{position + 1}. {member.mention} ({member.id})"
                for position, member in enumerate(chunk, start=index)
            )
            embed = discord.Embed(
                title=f"{role.name} Members",
                description=description or "No members found.",
                color=role.color if role.color.value else discord.Color.blurple(),
            )

            embeds.append(embed)
        if not embeds:
            embed = discord.Embed(
                title=f"{role.name} Members",
                description="No members found.",
                color=role.color if role.color.value else discord.Color.blurple(),
            )
            embeds.append(embed)
        return await tick.end(success=embeds)

    @app_commands.command(name="automute", description="Start/stop automute")
    @app_commands.describe(
        channel="Specify a channel ID/mention.", duration="Specify duration as m/h/d"
    )
    async def toggle_automute_app_command(
        self,
        interaction: discord.Interaction,
        channel: discord.abc.GuildChannel | None,
        duration: str | None = "1h",
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        resolved_channel = channel or interaction.channel
        if resolved_channel is None:
            return await tick.end(
                warning="This command must be executed in a valid channel."
            )
        await moderator_service.check_minimum_role(
            channel_snowflake=resolved_channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
            lowest_role="Coordinator",
        )
        pages = await automute_channel_service.toggle_automute(
            author_snowflake=interaction.user.id,
            channel_snowflake=resolved_channel.id,
            guild_snowflake=interaction.guild.id,
            duration_value=duration or "1h",
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorAppCommands(bot=bot))
