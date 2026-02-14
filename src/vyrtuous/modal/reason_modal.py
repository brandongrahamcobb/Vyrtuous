"""!/bin/python3
reason_modal.py The purpose of this program is to provide the reason utility modal.

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


class ReasonModal(discord.ui.Modal):
    def __init__(self, ctx, state):
        super().__init__(title="Reason")
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

    async def setup(self):
        self.reason_selection = discord.ui.TextInput(
            label="Type the reason",
            style=discord.TextStyle.paragraph,
            required=True,
            default=(self.ctx.record.reason if self.ctx.record else ""),
        )
        self.add_item(self.reason_selection)

    async def on_submit(self, interaction):
        await interaction.response.defer()
        self.state.interaction = interaction
        if self.ctx.record:
            set_kwargs = {"reason": self.reason_selection.value}
            await self.ctx.record.__class__.update(
                where_kwargs=self.record_snowflakes, set_kwargs=set_kwargs
            )
            await self.state.end(
                success=f"Existing infraction reason has been updated to {self.reason_selection.value}."
            )
        # else:
        #     for service in RecordService.__subclasses__():
        #         if isinstance(self.ctx.record, service.model):
        #             await service.enforce(
        #                 ctx=self.ctx, source=interaction, state=self.state
        #             )
