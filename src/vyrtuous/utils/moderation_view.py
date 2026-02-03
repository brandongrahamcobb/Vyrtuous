"""moderation_view.py The purpose of this program is to provide the moderation view utility class.

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

from pathlib import Path

import discord

from vyrtuous.fields.duration import DurationObject
from vyrtuous.utils.check import has_equal_or_lower_role_wrapper
from vyrtuous.utils.dir_to_classes import dir_to_classes


class ModerationView(discord.ui.View):

    infraction_paths = []
    infraction_paths.append(Path(__file__).resolve().parents[1] / "db/infractions")

    def __init__(self, interaction: discord.Interaction, member_snowflake: int, modal):
        super().__init__(timeout=120)
        self.information = {}
        self.author_snowflake = interaction.user.id
        self.infractions = dir_to_classes(dir_paths=self.infraction_paths)
        self.interaction = interaction
        self.member_snowflake = member_snowflake
        self.modal = modal

    async def interaction_check(self, interaction):
        return interaction.user.id == self.author_snowflake

    async def setup(self):
        channel_options = await self._build_channel_options()
        identifier_options = await self._build_identifier_options()
        if not identifier_options:
            await self.interaction.response.send_message(
                content="No moderation actions exist for this member.", ephemeral=True
            )
            self.stop()
            return
        self.channel_select.options = channel_options
        self.category_select.options = identifier_options

    async def _build_channel_options(self):
        channels = []
        for infraction in self.infractions:
            action = await infraction.select(
                guild_snowflake=self.interaction.guild.id,
                member_snowflake=self.member_snowflake,
                singular=True,
            )
            if action:
                channel = self.interaction.guild.get_channel(action.channel_snowflake)
                channels.append(channel)
        return [
            discord.SelectOption(label=ch.name, value=str(ch.id))
            for ch in channels
            if isinstance(ch, discord.VoiceChannel)
        ]

    async def _build_identifier_options(self):
        actions = []
        for infraction in self.infractions:
            action = await infraction.select(
                guild_snowflake=self.interaction.guild.id,
                member_snowflake=self.member_snowflake,
                singular=True,
            )
            if action:
                actions.append(infraction)
        return [
            discord.SelectOption(label=infraction.identifier, value=infraction.__name__)
            for infraction in actions
            if infraction.identifier != "smute"
        ]

    @discord.ui.select(
        placeholder="Select channel",
        options=[],
    )
    async def channel_select(self, interaction, select):
        channel = interaction.guild.get_channel(int(select.values[0]))
        self.information["channel_snowflake"] = channel.id
        self.channel_select.placeholder = channel.name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select category",
        options=[],
    )
    async def category_select(self, interaction, select):
        category_name = select.values[0]
        infraction = next(c for c in self.infractions if c.__name__ == category_name)
        self.category_select.placeholder = self.information["category"] = (
            infraction.identifier
        )
        self.information["infraction"] = infraction
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        executor_role = await has_equal_or_lower_role_wrapper(
            source=interaction,
            member_snowflake=self.member_snowflake,
            sender_snowflake=interaction.user.id,
        )
        existing = await self.information.get("infraction", None).select(
            channel_snowflake=self.information.get("channel_snowflake", None),
            member_snowflake=self.member_snowflake,
            singular=True,
        )
        self.information["executor_role"] = executor_role
        self.information["member_snowflake"] = self.member_snowflake
        self.information["existing"] = existing
        if hasattr(existing, "expires_in"):
            if DurationObject.from_expires_in_to_str(existing.expires_in) == 0:
                await interaction.response.send_message(
                    content="This moderation is permanent and can only be undone, not modified.",
                    ephemeral=True,
                )
                await interaction.message.delete()
                self.stop()
        modal = self.modal(information=self.information)
        await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        await interaction.message.delete()
        self.stop()
