"""!/bin/python3
duration_modal.py The purpose of this program is to provide the duration utility modal.

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


class DurationModal(discord.ui.Modal):
    def __init__(
        self, ctx, state, *, cap_service=None, duration=None, duration_service=None
    ):
        super().__init__(title="Duration")
        self.ctx = ctx
        self.channel_snowflake = ctx.target_channel_snowflake
        self.guild_snowflake = ctx.source_channel_snowflake
        self.member_snowflake = ctx.target_member_snowflake
        self.record_snowflakes = {
            "channel_snowflake": self.channel_snowflake,
            "guild_snowflake": self.guild_snowflake,
            "member_snowflake": self.member_snowflake,
        }
        self.state = state
        self.__duration = duration
        self.__duration_service = duration_service
        self.__cap_service = cap_service

    async def setup(self):
        self.duration_selection = discord.ui.TextInput(
            label="Type the duration",
            style=discord.TextStyle.paragraph,
            required=True,
            default=self.__duration_service.from_expires_in_to_str(
                self.ctx.record.expires_in
            )
            or "",
        )
        self.add_item(self.duration_selection)

    async def on_submit(self, interaction):
        await interaction.response.defer()
        self.state.interaction = interaction
        self.duration = self.__duration_service.parse(self.duration_selection.value)
        self.ctx.expires_in = (
            datetime.now(timezone.utc)
            + self.__duration_service.to_timedelta(duration=self.duration)
            if self.duration.number != 0
            else None
        )
        if await self.__cap_service.assert_duration_exceeds_cap(
            category=self.ctx.record.identifier,
            duration=self.duration,
            source_kwargs=self.ctx.source_kwargs,
        ):
            return await interaction.response.send_message(
                content="Duration exceeds the channel cap.", ephemeral=True
            )
        if self.ctx.record:
            set_kwargs = {"expires_in": self.ctx.expires_in}
            await self.ctx.record.__class__.update(
                where_kwargs=self.record_snowflakes, set_kwargs=set_kwargs
            )
            await self.state.end(
                success=f"Duration has been updated to {self.duration}."
            )
        # else:
        #     for service in RecordService.__subclasses__():
        #         if isinstance(self.ctx.record, service.model):
        #             await service.enforce(
        #                 ctx=self.ctx, source=interaction, state=self.state
        #             )
