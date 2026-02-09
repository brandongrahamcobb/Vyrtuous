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

from datetime import datetime, timezone

import discord

from vyrtuous.commands.fields.duration import DurationObject
from vyrtuous.commands.messaging.reason_modal import ReasonModal
from vyrtuous.db.mgmt.cap.cap_service import CapService


class InfractionView(discord.ui.View):
    def __init__(
        self,
        ctx,
        state,
    ):
        super().__init__(timeout=120)
        self.ctx = ctx
        self.duration: DurationObject | None = None
        self.selected_duration = None
        self.selected_channel_snowflake: int | None = None
        self.state = state

    async def interaction_check(self, interaction):
        return interaction.user.id == self.ctx.source_member_snowflake

    async def setup(self):
        channel_options = await self._build_channel_options()
        duration_options = self._build_duration_options()
        self.channel_select.options = channel_options
        self.duration_select.options = duration_options

    async def _build_channel_options(self):
        channel_options = [
            discord.SelectOption(label=c.name, value=str(c.id))
            for c in self.ctx.available_channels
            if c != "all"
        ]
        if "all" in self.ctx.available_channels:
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
        self.ctx.target_channel_snowflake = channel.id
        self.channel_select.placeholder = channel.name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select duration",
        options=[],
    )
    async def duration_select(self, interaction, select):
        duration_name = select.values[0]
        self.duration_select.placeholder = duration_name
        self.duration = DurationObject(duration_name)
        self.ctx.expires_in = (
            datetime.now(timezone.utc) + self.duration.to_timedelta()
            if self.duration.number != 0
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
        if await CapService.assert_duration_exceeds_cap(
            category=self.ctx.record.identifier,
            duration=self.duration,
            source_kwargs=self.ctx.source_kwargs,
        ):
            return await interaction.response.send_message(
                content="Duration exceeds the channel cap.", ephemeral=True
            )
        modal = ReasonModal(
            ctx=self.ctx,
            state=self.state,
        )
        await modal.setup()
        await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        return await self.state.end(success="Cancelled action.")

    def has_the_user_selected_all_fields(self):
        if not self.ctx.expires_in or not self.ctx.target_channel_snowflake:
            return False
        return True

        # executor_role = await PermissionService.resolve_highest_role(
        #     channel_snowflake=channel.id,
        #     guild_snowflake=interaction.guild.id,
        #     member_snowflake=interaction.user.id,
        # )
        # self.information["executor_role"] = executor_role
        # existing = await self.infraction.select(
        #     channel_snowflake=self.information.get("updated_kwargs", None).get(
        #         "channel_snowflake", None
        #     ),
        #     member_snowflake=self.member_snowflake,
        #     singular=True,
        # )
        # if existing:
        #     dir_paths = []
        #     dir_paths.append(Path("/app/vyrtuous/db/infractions"))
        #     infraction_services = dir_to_classes(
        #         dir_paths=dir_paths, parent=InfractionService
        #     )
        #     if service := next(
        #         (
        #             s
        #             for s in infraction_services
        #             if self.information["infraction"] is s.model
        #         ),
        #         None,
        #     ):
        #         await service.undo(
        #             information=self.information, source=interaction, state=self.state
        #         )
        #     self.stop()
        #     return

        # self.information.setdefault("updated_kwargs", {})
        # self.information["category"] = self.infraction.identifier
        # self.information["infraction"] = self.infraction
        # self.information["updated_kwargs"]["guild_snowflake"] = int(
        #     self.interaction.guild.id
        # )
        # self.information["updated_kwargs"]["member_snowflake"] = int(
        #     self.member_snowflake
        # )
