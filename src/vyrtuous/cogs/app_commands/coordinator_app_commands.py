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

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.roles.coordinator_service import coordinator_predicator
from vyrtuous.db.roles.moderator import Moderator
from vyrtuous.db.infractions.voice_mute import VoiceMute
from vyrtuous.fields.snowflake import (
    AppChannelSnowflake,
    AppMemberSnowflake,
)
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import DiscordObject


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
        msg = await Moderator.toggle_moderator(
            channel_dict=channel_dict,
            member_dict=member_dict,
            snowflake_kwargs=snowflake_kwargs,
        )
        return await state.end(success=msg)

    @app_commands.command(name="rmute", description="Room mute (except yourself).")
    @app_commands.describe(channel="Tag a channel or include its ID.")
    @coordinator_predicator()
    async def room_mute_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        reason: str = "No reason provided.",
    ):
        state = StateService(interaction=interaction)
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        pages = VoiceMute.room_mute(
            channel_dict=channel_dict,
            guild_snowflake=interaction.guild.id,
            reason=reason,
            snowflake_kwargs=snowflake_kwargs,
        )
        await StateService.send_pages(plural=VoiceMute.PLURAL, pages=pages, state=state)

    @app_commands.command(name="xrmute", description="Unmute all.")
    @app_commands.describe(channel="Tag a channel or include its ID.")
    @coordinator_predicator()
    async def room_unmute_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        pages = VoiceMute.room_unmute(
            channel_dict=channel_dict, guild_snowflake=interaction.guild.id
        )
        await StateService.send_pages(plural=VoiceMute.PLURAL, pages=pages, state=state)


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorAppCommands(bot))
