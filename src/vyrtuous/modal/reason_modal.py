"""!/bin/python3
reason_modal.py The purpose of this program is to provide the reason utility modal which is used to finalize infractions.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationObject
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.users import moderator_service

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
                raise commands.GuildNotFound(str(self.__guild_snowflake))
            channel = guild.get_channel(self.__channel_snowflake)
            if channel is None:
                raise commands.ChannelNotFound(str(self.__channel_snowflake))
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
        await interaction.response.defer()
        await moderator_service.has_equal_or_lower_role(
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
            member_snowflake=self.__author_snowflake,
            target_member_snowflake=self.__member_snowflake,
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
        embed = None
        match self.__category:
            case "ban":
                if self.__record:
                    await unban_alias_service.unban(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                    )
                    await unban_alias_service.log_unban(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=None,
                        message_channel_snowflake=interaction.channel.id,
                        reason=self.reason_selection.value,
                    )
                    embed = unban_alias_service.build_unban_embed(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                    )
                else:
                    if self.__duration is None:
                        return
                    await ban_alias_service.ban(
                        channel_snowflake=self.__channel_snowflake,
                        duration=self.__duration,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=self.reason_selection.value,
                    )
                    await ban_alias_service.log_ban(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        duration=self.__duration,
                        guild_snowflake=self.__guild_snowflake,
                        is_channel_scope=self.__is_channel_scope,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=None,
                        message_channel_snowflake=interaction.channel.id,
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
                if self.__record:
                    await unflag_alias_service.unflag(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                    )
                    await unflag_alias_service.log_unflag(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=None,
                        message_channel_snowflake=interaction.channel.id,
                        reason=self.reason_selection.value,
                    )
                    embed = unflag_alias_service.build_unflag_embed(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                    )
                else:
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
                        message_snowflake=None,
                        message_channel_snowflake=interaction.channel.id,
                        reason=self.reason_selection.value,
                    )
                    embed = flag_alias_service.build_flag_embed(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=self.reason_selection.value,
                    )
            case "tmute":
                if self.__record:
                    await untext_mute_alias_service.untext_mute(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                    )
                    await untext_mute_alias_service.log_untext_mute(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=None,
                        message_channel_snowflake=interaction.channel.id,
                        reason=self.reason_selection.value,
                    )
                    embed = untext_mute_alias_service.build_untext_mute_embed(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                    )
                else:
                    if self.__duration is None:
                        return
                    await text_mute_alias_service.text_mute(
                        channel_snowflake=self.__channel_snowflake,
                        duration=self.__duration,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=self.reason_selection.value,
                    )
                    await text_mute_alias_service.log_text_mute(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        duration=self.__duration,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=None,
                        message_channel_snowflake=interaction.channel.id,
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
                target = "command"
                if self.__record:
                    await unvoice_mute_alias_service.unvoice_mute(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        target=target,
                    )
                    await unvoice_mute_alias_service.log_unvoice_mute(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        guild_snowflake=self.__guild_snowflake,
                        is_channel_scope=self.__is_channel_scope,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=None,
                        message_channel_snowflake=interaction.channel.id,
                        reason=self.reason_selection.value,
                        target=target,
                    )
                    embed = unvoice_mute_alias_service.build_unvoice_mute_embed(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                    )
                else:
                    if self.__duration is None:
                        return
                    await voice_mute_alias_service.voice_mute(
                        channel_snowflake=self.__channel_snowflake,
                        duration=self.__duration,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=self.reason_selection.value,
                        target=target,
                    )
                    await voice_mute_alias_service.log_voice_mute(
                        author_snowflake=interaction.user.id,
                        channel_snowflake=self.__channel_snowflake,
                        display=True,
                        duration=self.__duration,
                        guild_snowflake=self.__guild_snowflake,
                        is_channel_scope=self.__is_channel_scope,
                        member_snowflake=self.__member_snowflake,
                        message_snowflake=None,
                        message_channel_snowflake=interaction.channel.id,
                        reason=self.reason_selection.value,
                        target=target,
                    )
                    embed = voice_mute_alias_service.build_voice_mute_embed(
                        channel_snowflake=self.__channel_snowflake,
                        duration=self.__duration,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=self.reason_selection.value,
                    )
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
