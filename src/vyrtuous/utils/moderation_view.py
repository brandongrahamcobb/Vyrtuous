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
from vyrtuous.utils.alias_information import AliasInformation
from vyrtuous.utils.check import has_equal_or_lower_role_wrapper
from vyrtuous.utils.dir_to_classes import dir_to_classes


class ModerationView(discord.ui.View):

    category_paths = [Path(__file__).resolve().parents[1] / "db/infractions"]

    def __init__(self, interaction: discord.Interaction, modal):
        super().__init__(timeout=120)
        self.information = {}
        self.author_snowflake = interaction.user.id
        self.interaction = interaction
        self.modal = modal
        self.categories = dir_to_classes(dir_paths=self.category_paths)

    async def interaction_check(self, interaction: discord.Interaction) -> bool:
        return interaction.user.id == self.author_snowflake

    async def setup(self):
        # Populate all fields from AliasInformation.build
        self.information = await AliasInformation.build(self.interaction.message)
        if not self.information:
            await self.interaction.response.send_message(
                "Failed to load alias information.", ephemeral=True
            )
            self.stop()
            return

        self.information["category"] = self.information["alias"]  # convenience
        self.channel_select.options = await self._build_channel_options()
        self.category_select.options = await self._build_category_options()

    async def _build_channel_options(self):
        channels = []
        category = self.information.get("category")
        if category:
            existing_action = await category.select(
                guild_snowflake=self.information["snowflake_kwargs"]["guild_snowflake"],
                member_snowflake=self.information["snowflake_kwargs"].get(
                    "member_snowflake"
                ),
                singular=True,
            )
            if existing_action:
                ch = self.interaction.guild.get_channel(
                    existing_action.channel_snowflake
                )
                if ch:
                    channels.append(ch)
        return [
            discord.SelectOption(label=ch.name, value=str(ch.id)) for ch in channels
        ]

    async def _build_category_options(self):
        category = self.information.get("category")
        if category and category.category != "smute":
            return [
                discord.SelectOption(label=category.category, value=category.__name__)
            ]
        return []

    @discord.ui.select(placeholder="Select channel", options=[])
    async def channel_select(
        self, interaction: discord.Interaction, select: discord.ui.Select
    ):
        self.information["snowflake_kwargs"]["channel_snowflake"] = int(
            select.values[0]
        )
        channel = interaction.guild.get_channel(int(select.values[0]))
        if channel:
            self.channel_select.placeholder = channel.name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(placeholder="Select category", options=[])
    async def category_select(
        self, interaction: discord.Interaction, select: discord.ui.Select
    ):
        selected_name = select.values[0]
        category_class = next(c for c in self.categories if c.__name__ == selected_name)
        self.information["category"] = category_class
        self.category_select.placeholder = category_class.category
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction: discord.Interaction, button: discord.ui.Button):
        executor_role = await has_equal_or_lower_role_wrapper(
            source=interaction,
            member_snowflake=self.information["snowflake_kwargs"].get(
                "member_snowflake"
            ),
            sender_snowflake=interaction.user.id,
        )
        self.information["executor_role"] = executor_role

        existing = await self.information["category"].select(
            channel_snowflake=self.information["snowflake_kwargs"].get(
                "channel_snowflake"
            ),
            member_snowflake=self.information["snowflake_kwargs"].get(
                "member_snowflake"
            ),
            singular=True,
        )
        self.information["existing"] = existing

        if hasattr(existing, "expires_in") and existing.expires_in is None:
            await interaction.response.send_message(
                "This moderation is permanent and can only be undone, not modified.",
                ephemeral=True,
            )
            await interaction.message.delete()
            self.stop()
            return

        modal = self.modal(information=self.information)
        await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction: discord.Interaction, button: discord.ui.Button):
        await interaction.message.delete()
        self.stop()
