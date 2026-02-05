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

from vyrtuous.commands.messaging.state_service import StateService


class ReasonModal(discord.ui.Modal):

    def __init__(self, information):
        super().__init__(
            title=f'{information.get("infraction", None).identifier.capitalize()} Reason'
        )
        self.information = information
        self.reason = discord.ui.TextInput(
            label="Type the reason",
            style=discord.TextStyle.paragraph,
            required=True,
            default=(
                self.information.get("existing", None).reason
                if self.information.get("existing", None)
                else ""
            ),
        )
        self.add_item(self.reason)

    async def on_submit(self, interaction):
        state = StateService(interaction=interaction)
        if self.information.get("existing", None):
            where_kwargs = {
                "channel_snowflake": self.information.get("channel_snowflake", None),
                "member_snowflake": self.information.get("member_snowflake", None),
            }
            set_kwargs = {"reason": self.reason.value}
            await self.information.get("infraction", None).update(
                where_kwargs=where_kwargs, set_kwargs=set_kwargs
            )
            await interaction.response.send_message(
                content="Reason has been updated.", ephemeral=True
            )
        else:
            where_kwargs = {
                "channel_snowflake": self.information.get("channel_snowflake", None),
                "expires_in": self.information.get("expires_in", None),
                "guild_snowflake": self.information.get("guild_snowflake", None),
                "member_snowflake": self.information.get("member_snowflake", None),
                "reason": self.reason.value,
            }
            infraction = self.information.get("infraction", None)(**where_kwargs)
            await infraction.create()
            await interaction.response.send_message(
                content=f"{self.information.get('infraction', None).identifier.capitalize()} updated.",
                ephemeral=True,
            )
