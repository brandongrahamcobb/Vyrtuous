"""!/bin/python3
duration_modal.py The purpose of this program is to provide the duration modal which is used for changing infraction durations.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import CheckFailure
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import cap_service


class DurationModal(discord.ui.Modal):
    def __init__(
        self,
        author_snowflake: int,
        model: type,
        record: type,
        tick: Tick,
    ):
        super().__init__(title="Duration", timeout=120)
        self.__author_snowflake = author_snowflake
        self.__database_factory: DatabaseFactory = DatabaseFactory(model)
        self.__model = model
        self.__record = record
        self.__tick = tick

    async def setup(self):
        duration_builder: DurationBuilder = DurationBuilder()
        default_duration_value = duration_builder.from_timestamp(
            self.__record.expires_in
        ).as_str()
        self.duration_selection = discord.ui.TextInput(
            label="Type the duration",
            style=discord.TextStyle.paragraph,
            required=True,
            default=str(default_duration_value) or "",
        )
        self.add_item(self.duration_selection)

    async def on_submit(self, interaction) -> None:
        self.__tick.update_source(interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        duration_builder: DurationBuilder = DurationBuilder()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=self.__record.channel_snowflake,
            guild_snowflake=self.__record.guild_snowflake,
            author_snowflake=self.__author_snowflake,
            member_snowflake=self.__record.member_snowflake,
        )
        expires_in = duration_builder.parse(
            self.duration_selection.value
        ).to_expires_in()
        duration = duration_builder.parse(self.duration_selection.value).build()
        if await cap_service.exceeds_cap(
            category=self.__model.identifier,
            channel_snowflake=self.__record.channel_snowflake,
            duration=duration,
            guild_snowflake=self.__record.guild_snowflake,
        ):
            try:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    channel_snowflake=self.__record.channel_snowflake,
                    guild_snowflake=self.__record.guild_snowflake,
                    member_snowflake=self.__author_snowflake,
                    requested=["command.moderation.uncapped"],
                )
            except CheckFailure:
                return await interaction.response.send_message(
                    content=f"Duration {duration_builder.parse(self.duration_selection.value).as_str()} exceeds the channel cap.",
                    ephemeral=True,
                )
        set_kwargs = {"expires_in": expires_in}
        where_kwargs = {
            "channel_snowflake": self.__record.channel_snowflake,
            "guild_snowflake": self.__record.guild_snowflake,
            "member_snowflake": self.__record.member_snowflake,
        }
        await self.__database_factory.update(
            where_kwargs=where_kwargs, set_kwargs=set_kwargs
        )
        await self.__tick.end(
            success=f"Duration has been updated to {duration_builder.parse(self.duration_selection.value).to_unix_ts()}.",
            ephemeral=True,
        )
        self.stop()

    async def on_error(
        self,
        interaction: discord.Interaction,
        error: Exception,
    ):
        await self.__tick.end(error=str(error))
