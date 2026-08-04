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
from discord.ext import commands

from vyrtuous.aliases import (
    ban_alias_service,
    flag_alias_service,
    text_mute_alias_service,
    unban_alias_service,
    unflag_alias_service,
    untext_mute_alias_service,
    unvoice_mute_alias_service,
    voice_mute_alias_service,
)
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import ChannelNotFound, GuildNotFound
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
        author_snowflake: int,
        category: str,
        channel_snowflake: int,
        duration: DurationObject | None,
        guild_snowflake: int,
        member_snowflake: int,
        tick: Tick,
        *,
        is_channel_scope: bool = False,
        is_modification: bool = False,
    ):
        super().__init__(title="Reason", timeout=120)
        self.__author_snowflake = author_snowflake
        self.__category = category
        self.__channel_snowflake = channel_snowflake
        self.__database_factory: DatabaseFactory
        self.__guild_snowflake = guild_snowflake
        self.__is_channel_scope = is_channel_scope
        self.__is_modification = is_modification
        self.__member_snowflake = member_snowflake
        self.__record = None
        self.__tick = tick
        self.__duration = duration

    async def setup(self, is_new: bool):
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
        if record:
            self.__record = record
        if not is_new and not record:
            bot: DiscordBot = DiscordBot.get_instance()
            guild = bot.get_guild(self.__guild_snowflake)
            if guild is None:
                raise GuildNotFound(str(self.__guild_snowflake))
            channel = guild.get_channel(self.__channel_snowflake)
            if channel is None:
                raise ChannelNotFound(str(self.__channel_snowflake))
            await self.__tick.end(
                warning="No infraction exists under this category ({self.__category}) for channel ({channel.mention}) in guild ({guild.name}).",
                ephemeral=True,
            )
        else:
            self.reason_selection = discord.ui.TextInput(
                label="Type the reason",
                style=discord.TextStyle.paragraph,
                required=True,
                default=(self.__record.reason if self.__record else ""),
            )
            self.add_item(self.reason_selection)

    async def on_submit(self, interaction) -> None:
        self.__tick.update_source(interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
            author_snowflake=self.__author_snowflake,
            member_snowflake=self.__member_snowflake,
        )
        if self.__is_modification:
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
                    success=f"Reason has been updated to {self.reason_selection.value}.",
                    ephemeral=True,
                )
                return None
            else:
                await self.__tick.end(
                    warning=f"No infraction exists for category ({self.__category}) on member ({self.__member_snowflake}).",
                    ephemeral=True,
                )
                return None
        if embed:
            await self.__tick.end(success=embed)
        self.stop()

    async def on_error(
        self,
        interaction: discord.Interaction,
        error: Exception,
    ):
        if isinstance(error, commands.BadArgument):
            await self.__tick.end(error=str(error))
        elif isinstance(error, commands.CheckFailure):
            await self.__tick.end(error=str(error))
        else:
            raise error
