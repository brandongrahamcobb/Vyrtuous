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

from vyrtuous.aliases import (
    ban_alias_service,
    flag_alias_service,
    text_mute_alias_service,
    voice_mute_alias_service,
)
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)

INFRACTION_MODELS = [
    ban_service.MODEL,
    flag_service.MODEL,
    text_mute_service.MODEL,
    voice_mute_service.MODEL,
]


class ReasonModal(discord.ui.Modal):
    def __init__(
        self,
        category: str,
        channel_snowflake: int,
        duration_value: str,
        guild_snowflake: int,
        member_snowflake: int,
        tick: Tick,
        *,
        is_channel_scope: bool = False,
    ):
        super().__init__(title="Reason", timeout=120)
        self.__category = category
        self.__channel_snowflake = channel_snowflake
        self.__database_factory: DatabaseFactory
        self.__duration_value = duration_value
        self.__guild_snowflake = guild_snowflake
        self.__is_channel_scope = is_channel_scope
        self.__member_snowflake = member_snowflake
        self.__tick = tick

    async def setup(self):
        model = next(
            (
                model
                for model in INFRACTION_MODELS
                if model.identifier == self.__category
            ),
            None,
        )
        self.__database_factory: DatabaseFactory = DatabaseFactory(model)
        record = await self.__database_factory.select(
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
            member_snowflake=self.__member_snowflake,
            singular=True,
        )
        self.__record = record
        self.reason_selection = discord.ui.TextInput(
            label="Type the reason",
            style=discord.TextStyle.paragraph,
            required=True,
            default=(self.__record.reason if self.__record else ""),
        )
        self.add_item(self.reason_selection)

    async def on_submit(self, interaction):
        await interaction.response.defer()
        if self.__record:
            set_kwargs = {"reason": self.reason_selection.value}
            where_kwargs = {
                "channel_snowflake": self.__channel_snowflake,
                "guild_snowflake": self.__guild_snowflake,
                "member_snowflake": self.__member_snowflake,
            }
            await self.__database_factory.update(
                where_kwargs=where_kwargs, set_kwargs=set_kwargs
            )
            await self.__tick.end(
                success=f"Existing infraction reason has been updated to {self.reason_selection.value}.",
            )
        else:
            match self.__category:
                case "ban":
                    await ban_alias_service.ban(
                        channel_snowflake=self.__channel_snowflake,
                        duration_value=self.__duration_value,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=self.reason_selection.value,
                    )
                    await ban_alias_service.log_ban(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        duration_value=self.__duration_value,
                        guild_snowflake=self.__guild_snowflake,
                        is_channel_scope=self.__is_channel_scope,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=await interaction.original_response().id,
                        message_channel_snowflake=interaction.channel.id,
                        reason=self.reason_selection.value,
                    )
                case "flag":
                    await flag_alias_service.flag(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=self.reason_selection.value,
                    )
                    await flag_alias_service.log_flag(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=await interaction.original_response().id,
                        message_channel_snowflake=interaction.channel.id,
                        reason=self.reason_selection.value,
                    )
                case "tmute":
                    await text_mute_alias_service.text_mute(
                        channel_snowflake=self.__channel_snowflake,
                        duration_value=self.__duration_value,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=self.reason_selection.value,
                    )
                    await text_mute_alias_service.log_text_mute(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        duration_value=self.__duration_value,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=await interaction.original_response().id,
                        message_channel_snowflake=interaction.channel.id,
                        reason=self.reason_selection.value,
                    )
                case "vmute":
                    target = "user"
                    await voice_mute_alias_service.voice_mute(
                        channel_snowflake=self.__channel_snowflake,
                        duration_value=self.__duration_value,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=self.reason_selection.value,
                        target=target,
                    )
                    await voice_mute_alias_service.log_voice_mute(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        duration_value=self.__duration_value,
                        guild_snowflake=self.__guild_snowflake,
                        is_channel_scope=self.__is_channel_scope,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=await interaction.original_response().id,
                        message_channel_snowflake=interaction.channel.id,
                        reason=self.reason_selection.value,
                        target=target,
                    )
