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


class InfractionView(discord.ui.View):
    def __init__(
        self,
        *,
        cap_service=None,
        ctx=None,
        duration_service=None,
        modal=None,
        state=None,
    ):
        super().__init__(timeout=120)
        self.__cap_service = cap_service
        self.__ctx = ctx
        self.__duration_service = duration_service
        self.__modal = modal
        self.__state = state

    async def interaction_check(self, interaction):
        return interaction.user.id == self.__ctx.source_member_snowflake

    async def setup(self):
        channel_options = await self._build_channel_options()
        duration_options = self._build_duration_options()
        self.channel_select.options = channel_options
        self.__duration_select.options = duration_options

    async def _build_channel_options(self):
        channel_options = [
            discord.SelectOption(label=c.name, value=str(c.id))
            for c in self.__ctx.available_channels
            if c != "all"
        ]
        if "all" in self.__ctx.available_channels:
            channel_options.append(discord.SelectOption(label="All", value="all"))
        return channel_options

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
        self.__ctx.target_channel_snowflake = channel.id
        self.channel_select.placeholder = channel.name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select duration",
        options=[],
    )
    async def duration_select(self, interaction, select):
        duration_name = select.values[0]
        self.__duration_select.placeholder = duration_name
        self.__duration = self.__duration_service.parse(duration=duration_name)
        self.__ctx.expires_in = (
            self.__duration_service.to_expires_in(duration=self.__duration)
            if self.__duration.number != 0
            else None
        )
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        if not self.has_the_user_selected_all_fields():
            return await interaction.response.send_message(
                content="Please select all fields.", ephemeral=True
            )
        if await self.__cap_duration.assert_duration_exceeds_cap(
            category=self.__ctx.infraction.identifier,
            duration=self.__duration,
            source_kwargs=self.__ctx.source_kwargs,
        ):
            return await interaction.response.send_message(
                content="Duration exceeds the channel cap.", ephemeral=True
            )
        modal = self.__modal(
            ctx=self.__ctx,
            state=self.__state,
        )
        await modal.setup()
        await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        await interaction.message.delete()
        self.stop()
        return await self.__state.end(success="Cancelled action.")

    def has_the_user_selected_all_fields(self):
        if not self.__ctx.expires_in or not self.__ctx.target_channel_snowflake:
            return False
        return True
