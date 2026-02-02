"""duration_modal.py The purpose of this program is to provide the duration utility modal.

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

from vyrtuous.db.mgmt.cap import Cap
from vyrtuous.fields.duration import DurationObject


from datetime import datetime, timezone
import discord
from vyrtuous.fields.duration import DurationObject


class DurationModal(discord.ui.Modal):

    def __init__(self, information):
        super().__init__(title=f'{information["alias"].category} Duration')
        self.information = information
        existing_duration = getattr(information.get("existing"), "expires_in", None)
        self.duration = discord.ui.TextInput(
            label="Type the duration",
            style=discord.TextStyle.paragraph,
            required=True,
            default=DurationObject.from_expires_in_to_str(existing_duration) or "",
        )
        self.add_item(self.duration)

    async def on_submit(self, interaction: discord.Interaction):
        duration_obj = DurationObject(self.duration.value)
        cap_seconds = self.information.get(
            "cap_duration", DurationObject("8h").to_seconds()
        )

        if (
            duration_obj.to_timedelta().total_seconds() > cap_seconds
            and self.information.get("executor_role") == "Moderator"
        ):
            duration_str = DurationObject.from_seconds(cap_seconds)
            channel = interaction.guild.get_channel(
                self.information["snowflake_kwargs"]["channel_snowflake"]
            )
            await interaction.response.send_message(
                f"Cannot set {self.information['alias'].category} beyond {duration_str} as a Moderator in {channel.mention}.",
                ephemeral=True,
            )
            return

        expires_in = (
            None
            if duration_obj.number == 0
            else datetime.now(timezone.utc) + duration_obj.to_timedelta()
        )

        await self.information["category"].update(
            where_kwargs={
                "channel_snowflake": self.information["snowflake_kwargs"][
                    "channel_snowflake"
                ],
                "member_snowflake": self.information["snowflake_kwargs"][
                    "member_snowflake"
                ],
            },
            set_kwargs={"expires_in": expires_in},
        )
        await interaction.response.send_message(
            f"Duration has been updated to {duration_obj}.", ephemeral=True
        )
