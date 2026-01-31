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

    def __init__(self, infraction_information):
        super().__init__(title=f'{infraction_information["alias_class"].SINGULAR} Reason')
        self.infraction_information = infraction_information
        self.reason = discord.ui.TextInput(
            label="Type the reason",
            style=discord.TextStyle.paragraph,
            required=True,
            default=self.infraction_information.get("infraction_existing", None).reason or "",
        )
        self.add_item(self.reason)

    async def on_submit(self, interaction):
        where_kwargs = {
            "channel_snowflake": self.infraction_information.get(
                "infraction_channel_snowflake", None
            ),
            "member_snowflake": self.infraction_information.get(
                "infraction_member_snowflake", None
            ),
        }
        set_kwargs = {"reason": self.reason.value}
        await self.infraction_information.get("alias_class", None).update(
            where_kwargs=where_kwargs, set_kwargs=set_kwargs
        )
        await interaction.response.send_message(
            content="Reason has been updated.", ephemeral=True
        )
