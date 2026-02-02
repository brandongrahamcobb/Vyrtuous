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


class DurationModal(discord.ui.Modal):

    def __init__(self, information):
        super().__init__(title=f'{information["category"].category.capitalize()} Duration')
        self.information = information
        self.duration = discord.ui.TextInput(
            label="Type the duration",
            style=discord.TextStyle.paragraph,
            required=True,
            default=DurationObject.from_expires_in_to_str(
                self.information.get("existing", None).expires_in
            )
            or "",
        )
        self.add_item(self.duration)

    async def on_submit(self, interaction):
        channel = interaction.guild.get_channel(
            self.information.get("channel_snowflake", None)
        )
        duration_obj = DurationObject(self.duration.value)
        cap = await Cap.select(
            category=self.information.get("category", None).category,
            channel_snowflake=self.information.get(
                "channel_snowflake", None
            ),
            guild_snowflake=interaction.guild.id,
            singular=True,
        )
        if cap:
            channel_cap = cap.duration_seconds
        else:
            channel_cap = DurationObject("8h").to_seconds()
        expires_in_timedelta = DurationObject(self.duration.value).to_timedelta()
        if (
            self.information["existing"]
            and expires_in_timedelta.total_seconds() > channel_cap
            or duration_obj.number == 0
        ):
            if self.information.get("executor_role", None) == "Moderator":
                duration_str = DurationObject.from_seconds(channel_cap)
                await interaction.response.send_message(
                    content=f"Cannot set the "
                    f"{self.information['category'].category.capitalize()} beyond {duration_str} as a "
                    f"{self.information.get("executor_role", None)} in {channel.mention}."
                )
        where_kwargs = {
            "channel_snowflake": self.information.get(
                "channel_snowflake", None
            ),
            "member_snowflake": self.information.get(
                "member_snowflake", None
            ),
        }
        set_kwargs = {
            "expires_in": (
                datetime.now(timezone.utc)
                + DurationObject(self.duration.value).to_timedelta()
                if duration_obj.number != 0
                else None
            )
        }
        await self.information.get("category", None).update(
            where_kwargs=where_kwargs, set_kwargs=set_kwargs
        )
        await interaction.response.send_message(
            content=f"Duration has been updated to {DurationObject(self.duration.value)}.",
            ephemeral=True,
        )