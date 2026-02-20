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
        self,
        *,
        ctx=None,
        state=None,
        cap_service=None,
        default_ctx=None,
        duration_builder=None,
    ):
        super().__init__(title="Duration")
        self.__ctx = ctx
        self.__channel_snowflake = ctx.channel.id
        self.__guild_snowflake = ctx.source_channel_snowflake
        self.__member_snowflake = ctx.member.id
        self.__record_snowflakes = {
            "channel_snowflake": self.__channel_snowflake,
            "guild_snowflake": self.__guild_snowflake,
            "member_snowflake": self.__member_snowflake,
        }
        self.__state = state
        self.__duration_builder = duration_builder
        self.__cap_service = cap_service
        self.__d_ctx = default_ctx

    async def setup(self):
        self.duration_selection = discord.ui.TextInput(
            label="Type the duration",
            style=discord.TextStyle.paragraph,
            required=True,
            default=self.__duration_builder.from_timestamp(
                self.__ctx.record.expires_in
            ).build(as_str=True)
            or "",
        )
        self.add_item(self.duration_selection)

    async def on_submit(self, interaction):
        await interaction.response.defer()
        self.__state.interaction = interaction
        duration_str = self.__duration_builder.parse(
            self.duration_selection.value
        ).build(as_str=True)
        self.__ctx.duration = self.__duration_builder.parse(
            self.duration_selection.value
        ).build()
        self.__ctx.expires_in = self.__duration_builder.parse(
            self.duration_selection.value
        ).to_expires_in()
        if await self.__cap_service.assertion(ctx=self.__ctx, default_ctx=self.__d_ctx):
            return await interaction.response.send_message(
                content=f"Duration {duration_str} exceeds the channel cap.",
                ephemeral=True,
            )
        if self.__ctx.record:
            set_kwargs = {"expires_in": self.__ctx.expires_in}
            await self.__ctx.record.__class__.update(
                where_kwargs=self.__record_snowflakes, set_kwargs=set_kwargs
            )
            await self.__state.end(
                success=f"Duration has been updated to {duration_str}."
            )
        # else:
        #     for service in RecordService.__subclasses__():
        #         if isinstance(self.__ctx.record, service.model):
        #             await service.enforce(
        #                 ctx=self.__ctx, source=interaction, state=self.__state
        #             )
