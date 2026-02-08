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
    def __init__(self, infraction, infraction_service, selections, state):
        super().__init__(title="Reason")
        self.channel_snowflake = selections["channel_snowflake"]
        self.guild_snowflake = selections["guild_snowflake"]
        self.member_snowflake = selections["member_snowflake"]
        self.record_snowflakes = {
            "channel_snowflake": self.channel_snowflake,
            "guild_snowflake": self.guild_snowflake,
            "member_snowflake": self.member_snowflake,
        }
        self.infraction = infraction
        self.infraction_service = infraction_service
        self.selections = selections
        self.state = state

    async def setup(self):
        self.record = await self.infraction.select(**self.record_snowflakes)
        self.reason = discord.ui.TextInput(
            label="Type the reason",
            style=discord.TextStyle.paragraph,
            required=True,
            default=(self.record.reason if self.record else ""),
        )
        self.add_item(self.reason)

    async def on_submit(self, interaction):
        self.state.interaction = interaction
        if self.record:
            set_kwargs = {"reason": self.reason.value}
            await self.infraction.update(
                where_kwargs=self.record_snowflakes, set_kwargs=set_kwargs
            )
            await self.state.end(success="Existing infraction reason has been updated")
        else:
            await self.infraction_service.enforce(state=self.state)
