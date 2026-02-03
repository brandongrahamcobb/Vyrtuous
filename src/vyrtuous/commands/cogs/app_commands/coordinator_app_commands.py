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
from vyrtuous.commands.discord_object_service import DiscordObject
from vyrtuous.commands.fields.duration import AppDuration
from vyrtuous.commands.fields.snowflake import AppChannelSnowflake, AppMemberSnowflake
from vyrtuous.commands.messaging.message_service import MessageService
from vyrtuous.commands.messaging.state_service import StateService
from vyrtuous.db.roles.coord.coordinator_service import coordinator_predicator
from vyrtuous.db.roles.mod.moderator_service import ModeratorService
from vyrtuous.db.rooms.stage.stage_service import StageService


class CoordinatorAppCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    @app_commands.command(name="mod", description="Grant/revoke mods.")
    @app_commands.describe(
        member="Tag a member or include their ID",
        channel="Tag a channel or include its ID.",
    )
    @coordinator_predicator()
    async def toggle_moderator_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake,
    ):
        state = StateService(interaction=interaction)
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        msg = await ModeratorService.toggle_moderator(
            channel_dict=channel_dict,
            member_dict=member_dict,
            snowflake_kwargs=snowflake_kwargs,
        )
        return await state.end(success=msg)

    @app_commands.command(name="stage", description="Start/stop stage.")
    @app_commands.describe(
        channel="Tag a voice/stage channel",
        duration="Duration of the stage (e.g., 1h, 30m)",
    )
    @coordinator_predicator()
    async def toggle_stage_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        duration: AppDuration,
    ):
        state = StateService(interaction=interaction)
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        do = DiscordObject(interaction=interaction)
        channel = channel or int(interaction.channel.id)
        channel_dict = await do.determine_from_target(target=channel)
        pages = await StageService.toggle_stage(
            channel_dict=channel_dict,
            duration=duration,
            snowflake_kwargs=snowflake_kwargs,
        )
        await StateService.send_pages(title="Stages", pages=pages, state=state)


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorAppCommands(bot))
