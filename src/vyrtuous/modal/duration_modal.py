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
from discord.ext import commands
from discord.ext.commands.errors import ChannelNotFound

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import CheckFailure, GuildNotFound
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    cap_service,
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


class DurationModal(discord.ui.Modal):
    def __init__(
        self,
        author_snowflake: int,
        category: str,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
        tick: Tick,
    ):
        super().__init__(title="Duration", timeout=120)
        self.__author_snowflake = author_snowflake
        self.__category = category
        self.__channel_snowflake = channel_snowflake
        self.__guild_snowflake = guild_snowflake
        self.__member_snowflake = member_snowflake
        self.__tick = tick

    async def setup(self):
        duration_builder: DurationBuilder = DurationBuilder()
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
            default_duration_value = duration_builder.from_timestamp(
                record.expires_in
            ).as_str()
            self.__record = record
            self.duration_selection = discord.ui.TextInput(
                label="Type the duration",
                style=discord.TextStyle.paragraph,
                required=True,
                default=str(default_duration_value) or "",
            )
            self.add_item(self.duration_selection)
        else:
            bot: DiscordBot = DiscordBot.get_instance()
            guild = bot.get_guild(self.__guild_snowflake)
            if guild is None:
                raise GuildNotFound(str(self.__guild_snowflake))
            channel = guild.get_channel(self.__channel_snowflake)
            if channel is None:
                raise ChannelNotFound(str(self.__channel_snowflake))
            await self.__tick.end(
                warning=f"No infraction exists under this category ({self.__category}) for channel ({channel.mention}) in guild ({guild.name}).",
                ephemeral=True,
            )

    async def on_submit(self, interaction):
        await interaction.response.defer()
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        duration_builder: DurationBuilder = DurationBuilder()
        expires_in = duration_builder.parse(
            self.duration_selection.value
        ).to_expires_in()
        duration = duration_builder.parse(self.duration_selection.value).build()
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
            author_snowflake=self.__author_snowflake,
            member_snowflake=self.__member_snowflake,
        )
        if await cap_service.exceeds_cap(
            category=self.__category,
            channel_snowflake=self.__channel_snowflake,
            duration=duration,
            guild_snowflake=self.__guild_snowflake,
        ):
            try:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    requested=["command.moderation.uncapped"],
                )
            except CheckFailure:
                return await interaction.response.send_message(
                    content=f"Duration {duration_builder.parse(self.duration_selection.value).as_str()} exceeds the channel cap.",
                    ephemeral=True,
                )
        if self.__record:
            set_kwargs = {"expires_in": expires_in}
            where_kwargs = {
                "channel_snowflake": self.__channel_snowflake,
                "guild_snowflake": self.__guild_snowflake,
                "member_snowflake": self.__member_snowflake,
            }
            await self.__database_factory.update(
                where_kwargs=where_kwargs, set_kwargs=set_kwargs
            )
            await self.__tick.end(
                success=f"Duration has been updated to {self.duration_selection.value}.",
                ephemeral=True,
            )
            return None
        else:
            await self.__tick.end(
                warning=f"No infraction exists for category ({self.__category}) on member ({self.__member_snowflake}).",
                ephemeral=True,
            )
            return None

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
