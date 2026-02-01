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

    category_paths = []
    category_paths.append(Path(__file__).resolve().parents[1] / "db/actions")

    def __init__(self, interaction: discord.Interaction, member_snowflake: int, modal):
        super().__init__(timeout=120)
        self.infraction_information = {}
        self.author_snowflake = interaction.user.id
        self.categories = dir_to_classes(dir_paths=self.category_paths)
        self.interaction = interaction
        self.member_snowflake = member_snowflake
        self.modal = modal

    async def interinfraction_check(self, interaction):
        return interaction.user.id == self.author_snowflake

    async def setup(self):
        self.channel_select.options = await self._build_channel_options()
        self.category_select.options = await self._build_category_options()

    async def _build_channel_options(self):
        channels = []
        for category in self.categories:
            action = await category.select(
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

    async def _build_category_options(self):
        actions = []
        for category in self.categories:
            action = await category.select(
                guild_snowflake=self.interaction.guild.id,
                member_snowflake=self.member_snowflake,
                singular=True,
            )
            if action:
                actions.append(category)
        return [
            discord.SelectOption(label=category.ACT, value=category.__name__)
            for category in actions
            if category.ACT != "smute"
        ]

    @discord.ui.select(
        placeholder="Select channel",
        options=[],
    )
    async def channel_select(self, interaction, select):
        channel = interaction.guild.get_channel(int(select.values[0]))
        self.information["snowflake_kwargs"]["channel_snowflake"] = channel.id
        self.channel_select.placeholder = channel.name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select category",
        options=[],
    )
    async def category_select(self, interaction, select):
        category_name = select.values[0]
        category = next(c for c in self.categories if c.__name__ == category_name)
        self.category_select.placeholder = category.ACT
        self.infraction_information["alias_class"] = category
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        executor_role = await has_equal_or_lower_role_wrapper(
            source=interaction,
            member_snowflake=self.member_snowflake,
            sender_snowflake=interaction.user.id,
        )
        infraction_existing = await self.infraction_information.get(
            "alias_class", None
        ).select(
            channel_snowflake=self.infraction_information.get(
                "infraction_channel_snowflake", None
            ),
            member_snowflake=self.member_snowflake,
            singular=True,
        )
        self.information["executor_role"] = executor_role
        self.information["snowflake_kwargs"]["member_snowflake"] = self.member_snowflake
        self.information["existing"] = infraction_existing
        if hasattr(infraction_existing, "expires_in"):
            if (
                DurationObject.from_expires_in_to_str(infraction_existing.expires_in)
                == 0
            ):
                await interaction.response.send_message(
                    content="This moderation is permanent and can only be undone, not modified.",
                    ephemeral=True,
                )
                await interaction.message.delete()
                self.stop()
        modal = self.modal(infraction_information=self.infraction_information)
        await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        await interaction.message.delete()
        self.stop()
