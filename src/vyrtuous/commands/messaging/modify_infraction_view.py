""" "!/bin/python3
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

from vyrtuous.commands.messaging.view_context import ViewContext
from vyrtuous.db.infractions.ban.ban import Ban
from vyrtuous.db.infractions.flag.flag import Flag
from vyrtuous.db.infractions.tmute.text_mute import TextMute
from vyrtuous.db.infractions.vmute.voice_mute import VoiceMute


class ModifyInfractionView(discord.ui.View):
    def __init__(
        self,
        ctx,
        modal,
        state,
    ):
        super().__init__(timeout=120)
        self.ctx = ctx
        self.infractions = [Ban, Flag, TextMute, VoiceMute]
        self.modal = modal
        self._record_map = {}
        self.state = state

    async def interaction_check(self, interaction):
        return interaction.user.id == self.ctx.source_member_snowflake

    async def setup(self):
        channel_options = await self._build_channel_options()
        self.channel_select.options = channel_options

    async def _build_channel_options(self):
        channel_options = [
            discord.SelectOption(label=c.name, value=str(c.id))
            for c in self.ctx.available_channels
            if c != "all"
        ]
        if "all" in self.ctx.available_channels:
            channel_options.append(discord.SelectOption(label="All", value="all"))
        return channel_options

    async def _build_infraction_options(self):
        for infraction in self.infractions:
            record = await infraction.select(
                channel_snowflake=self.ctx.target_channel_snowflake,
                guild_snowflake=self.ctx.source_guild_snowflake,
                member_snowflake=self.ctx.target_member_snowflake,
                singular=True,
            )
            if record:
                key = str(record.identifier)
                self._record_map[key] = record
        return [
            discord.SelectOption(label=r.identifier, value=k)
            for k, r in self._record_map.items()
        ]

    @discord.ui.select(
        placeholder="Select channel",
        options=[],
    )
    async def channel_select(self, interaction, select):
        channel = interaction.guild.get_channel(int(select.values[0]))
        self.ctx.target_channel_snowflake = channel.id
        self.channel_select.placeholder = channel.name
        infraction_options = await self._build_infraction_options()
        self.infraction_select.options = infraction_options
        self.infraction_select.disabled = False
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select an infraction",
        options=[discord.SelectOption(label="Select channel first", value="noop")],
        disabled=True,
    )
    async def infraction_select(self, interaction, select):
        key = select.values[0]
        self.ctx.record = self._record_map[key]
        self.infraction_select.placeholder = self.ctx.record.identifier
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        if not self.has_the_user_selected_all_fields():
            return await interaction.response.send_message(
                content="Please select all fields.", ephemeral=True
            )
        modal = self.modal(ctx=self.ctx, state=self.state)
        await modal.setup()
        await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        self.stop()
        return await self.state.end(success="Cancelled action.")

    def has_the_user_selected_all_fields(self):
        if not self.ctx.target_channel_snowflake or not self.ctx.record:
            return False
        return True
