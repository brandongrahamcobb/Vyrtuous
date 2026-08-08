"""!/bin/python3
reason_modal.py The purpose of this program is to provide the reason utility modal which is used to finalize infractions.

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
from vyrtuous.permissions import permission_service
from vyrtuous.utils.messaging.tick import Tick


class ReasonModal(discord.ui.Modal):
    def __init__(
        self,
        author_snowflake: int,
        model: type,
        record: type,
        tick: Tick,
    ):
        super().__init__(title="Reason", timeout=120)
        self.__author_snowflake = author_snowflake
        self.__database_factory: DatabaseFactory = DatabaseFactory(model)
        self.__record = record
        self.__tick = tick

    async def setup(self):
        self.reason_selection = discord.ui.TextInput(
            label="Type the reason",
            style=discord.TextStyle.paragraph,
            required=True,
            default=self.__record.reason,
        )
        self.add_item(self.reason_selection)

    async def on_submit(self, interaction) -> None:
        self.__tick.update_source(interaction=interaction)
        await self.__tick.defer()
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=self.__record.channel_snowflake,
            guild_snowflake=self.__record.guild_snowflake,
            author_snowflake=self.__author_snowflake,
            member_snowflake=self.__record.member_snowflake,
        )
        set_kwargs = {"reason": self.reason_selection.value}
        where_kwargs = {
            "channel_snowflake": self.__record.channel_snowflake,
            "guild_snowflake": self.__record.guild_snowflake,
            "member_snowflake": self.__record.member_snowflake,
        }
        await self.__database_factory.update(
            where_kwargs=where_kwargs, set_kwargs=set_kwargs
        )
        await self.__tick.end(
            success=f"Reason has been updated to {self.reason_selection.value}.",
            ephemeral=True,
        )
        self.stop()

    async def on_error(
        self,
        interaction: discord.Interaction,
        error: Exception,
    ):
        await self.__tick.end(error=str(error))
