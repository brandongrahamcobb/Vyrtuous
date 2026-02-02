"""reason_modal.py The purpose of this program is to provide the reason utility modal.

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

    def __init__(self, information):
        super().__init__(title=f'{information["alias"].category} Reason')
        self.information = information
        self.reason = discord.ui.TextInput(
            label="Type the reason",
            style=discord.TextStyle.paragraph,
            required=True,
            default=self.information.get("reason") or "",
        )
        self.add_item(self.reason)

    async def on_submit(self, interaction: discord.Interaction):
        await self.information["alias"].update(
            where_kwargs={
                "channel_snowflake": self.information["snowflake_kwargs"][
                    "channel_snowflake"
                ],
                "member_snowflake": self.information["snowflake_kwargs"][
                    "member_snowflake"
                ],
            },
            set_kwargs={"reason": self.reason.value},
        )
        await interaction.response.send_message(
            "Reason has been updated.", ephemeral=True
        )
