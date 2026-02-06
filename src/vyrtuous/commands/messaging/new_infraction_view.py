"""!/bin/python3
new_infraction_view.py The purpose of this program is to provide the view for creating an infraction.

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

from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.roles.admin.administrator import Administrator
from vyrtuous.db.roles.coord.coordinator import Coordinator
from vyrtuous.db.roles.dev.developer import Developer
from vyrtuous.db.roles.mod.moderator import Moderator
from vyrtuous.db.roles.owner.guild_owner_service import (
    NotGuildOwner,
    is_guild_owner_wrapper,
)
from vyrtuous.db.roles.sysadmin.sysadmin_service import NotSysadmin, is_sysadmin_wrapper


class NewInfractionView(discord.ui.View):

    def __init__(
        self, interaction: discord.Interaction, infraction, member_snowflake: int, modal
    ):
        super().__init__(timeout=120)
        self.information = {}
        self.author_snowflake = int(interaction.user.id)
        self.infraction = infraction
        self.interaction = interaction
        self.member_snowflake = int(member_snowflake)
        self.modal = modal

    async def interaction_check(self, interaction):
        return interaction.user.id == self.author_snowflake

    async def setup(self):
        self.information["category"] = self.infraction.identifier
        self.information["infraction"] = self.infraction
        self.information["guild_snowflake"] = int(self.interaction.guild.id)
        channel_options = await self._build_channel_options()
        duration_options = self._build_duration_options()
        self.channel_select.options = channel_options
        self.duration_select.options = duration_options

    async def _build_channel_options(self):
        available_channel_snowflakes = set()
        available_channels = set()
        default_kwargs = {
            "guild_snowflake": int(self.interaction.guild.id),
            "member_snowflake": self.author_snowflake,
        }
        moderators = await Moderator.select(**default_kwargs)
        if moderators:
            for moderator in moderators:
                available_channel_snowflakes.add(int(moderator.channel_snowflake))
        coordinators = await Coordinator.select(**default_kwargs)
        if coordinators:
            for coordinator in coordinators:
                available_channel_snowflakes.add(int(coordinator.channel_snowflake))
        administrator = await Administrator.select(**default_kwargs)
        guild_owner = None
        try:
            guild_owner = await is_guild_owner_wrapper(source=self.interaction)
        except NotGuildOwner:
            pass
        del default_kwargs["guild_snowflake"]
        developer = await Developer.select(**default_kwargs)
        sysadmin = None
        try:
            sysadmin = await is_sysadmin_wrapper(source=self.interaction)
        except NotSysadmin:
            pass
        if administrator or developer or guild_owner or sysadmin:
            for channel in self.interaction.guild.channels:
                if isinstance(channel, discord.VoiceChannel):
                    available_channels.add(channel)
        for channel_snowflake in available_channel_snowflakes:
            channel = self.interaction.guild.get_channel(channel_snowflake)
            if channel:
                available_channels.add(channel)
        return [
            discord.SelectOption(label=ch.name, value=str(ch.id))
            for ch in available_channels
        ]

    def _build_duration_options(self):
        durations = ["0", "1h", "8h", "1d", "1w"]
        return [
            discord.SelectOption(
                label=duration if duration != "0" else "Permanent", value=duration
            )
            for duration in durations
        ]

    @discord.ui.select(
        placeholder="Select channel",
        options=[],
    )
    async def channel_select(self, interaction, select):
        channel = interaction.guild.get_channel(int(select.values[0]))
        self.information["channel_snowflake"] = channel.id
        existing = await self.infraction.select(
            channel_snowflake=self.information.get("channel_snowflake", None),
            member_snowflake=self.member_snowflake,
            singular=True,
        )
        if existing:
            await interaction.response.send_message(
                content=f"An {self.infraction.category} already exists in this channel for this member.",
                ephemeral=True,
            )
            self.stop()
            return
        self.channel_select.placeholder = channel.name
        executor_role = await PermissionService.resolve_highest_role(
            channel_snowflake=channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
        )
        self.information["executor_role"] = executor_role
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select duration",
        options=[],
    )
    async def duration_select(self, interaction, select):
        duration_name = select.values[0]
        self.duration_select.placeholder = self.information["duration"] = duration_name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        self.information["member_snowflake"] = self.member_snowflake
        modal = self.modal(information=self.information)
        if (
            "duration" not in self.information
            or "channel_snowflake" not in self.information
        ):
            await interaction.response.send_message(
                content="You must select both a channel and a duration.", ephemeral=True
            )
            return
        await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        await interaction.message.delete()
        self.stop()
