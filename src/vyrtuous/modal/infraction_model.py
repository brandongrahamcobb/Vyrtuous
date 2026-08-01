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

from enum import Flag

import discord
from discord.ext import commands

from vyrtuous.aliases import (
    ban_alias_service,
    flag_alias_service,
    text_mute_alias_service,
    voice_mute_alias_service,
)
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.ban import Ban
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.text_mute import TextMute
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.tracking import stream_service


class InfractionModal(discord.ui.Modal):
    def __init__(
        self,
        author_snowflake: int,
        channel_snowflake: int,
        duration: DurationObject,
        guild_snowflake: int,
        member_snowflake: int,
        model: type,
        tick: Tick,
    ):
        super().__init__(title="Reason", timeout=120)
        self.__author_snowflake = author_snowflake
        self.__channel_snowflake = channel_snowflake
        self.__duration = duration
        self.__guild_snowflake = guild_snowflake
        self.__member_snowflake = member_snowflake
        self.__model = model
        self.__tick = tick

    async def setup(self):
        self.reason_selection = discord.ui.TextInput(
            label="Type the reason",
            style=discord.TextStyle.paragraph,
            required=True,
            default="",
        )
        self.add_item(self.reason_selection)

    async def on_submit(self, interaction) -> None:
        bot: DiscordBot = DiscordBot.get_instance()
        database_factory: DatabaseFactory = DatabaseFactory(self.__model)
        duration_builder: DurationBuilder = DurationBuilder()
        expires_in = (
            duration_builder.load(self.__duration).to_expires_in()
            if self.__duration
            else None
        )
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
            author_snowflake=self.__author_snowflake,
            member_snowflake=self.__member_snowflake,
        )
        display = True
        embed = None
        is_channel_scope = False
        message = interaction.message
        target = None
        match self.__model.identifier:
            case "ban":
                ban = Ban(
                    channel_snowflake=self.__channel_snowflake,
                    expires_in=expires_in,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                )
                await database_factory.create(ban)
                is_channel_scope = await ban_alias_service.enable(
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                )
                embed = ban_alias_service.build_ban_embed(
                    channel_snowflake=self.__channel_snowflake,
                    duration=self.__duration,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                )
            case "flag":
                flag = Flag(
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                )
                await database_factory.create(flag)
                await flag_alias_service.enable(
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                )
                embed = flag_alias_service.build_flag_embed(
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                )
            case "tmute":
                tmute = TextMute(
                    channel_snowflake=self.__channel_snowflake,
                    expires_in=expires_in,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                )
                await database_factory.create(tmute)
                await text_mute_alias_service.enable(
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                )
                embed = text_mute_alias_service.build_text_mute_embed(
                    channel_snowflake=self.__channel_snowflake,
                    duration=self.__duration,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                )
            case "vmute":
                automute_database_factory: DatabaseFactory = DatabaseFactory(AutoMute)
                automute = await automute_database_factory.select(
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__guild_snowflake,
                    singular=True,
                )
                if automute:
                    target = "auto"
                else:
                    target = "command"
                vmute = VoiceMute(
                    channel_snowflake=self.__channel_snowflake,
                    expires_in=expires_in,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                    target=target,
                )
                await database_factory.create(vmute)
                await voice_mute_alias_service.enable(
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                )
                embed = voice_mute_alias_service.build_voice_mute_embed(
                    channel_snowflake=self.__channel_snowflake,
                    duration=self.__duration,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__member_snowflake,
                    reason=self.reason_selection.value,
                )
        await stream_service.log(
            author_snowflake=self.__author_snowflake,
            channel_snowflake=self.__channel_snowflake,
            display=display,
            duration=DurationObject(number=0, prefix="", sign=1, unit=""),
            guild_snowflake=self.__guild_snowflake,
            identifier=self.__model.identifier,
            is_channel_scope=is_channel_scope,
            member_snowflake=self.__member_snowflake,
            message_snowflake=None,
            message_channel_snowflake=None,
            reason=self.reason_selection.value,
            role_snowflake=None,
            target=target,
        )
        await interaction.response.defer()
        await interaction.delete_original_response()
        if embed:
            await self.__tick.end(success=embed)
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
