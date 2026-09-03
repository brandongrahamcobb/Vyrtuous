"""!/bin/python3
cancel_confirm.py The purpose of this program is to provide an embed utility with cancellation and confirmation buttons.

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

from dataclasses import dataclass
from typing import Any, Union

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.db.autoassign import AutoAssignRole
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.modal.duration_modal import DurationModal
from vyrtuous.modal.reason_modal import ReasonModal
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import CheckFailure, GuildNotFound, MemberNotFound
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)

MODELS = [
    ban_service.MODEL,
    flag_service.MODEL,
    text_mute_service.MODEL,
    voice_mute_service.MODEL,
]
PERMISSIONS = {
    "command.moderation.ban": ban_service.MODEL,
    "command.moderation.flag": flag_service.MODEL,
    "command.moderation.text-mute": text_mute_service.MODEL,
    "command.moderation.voice-mute": voice_mute_service.MODEL,
}

SCOPES = {
    "command.moderation.voice-mute.auto": "Auto",
    "command.moderation.voice-mute.click": "Click",
    "command.moderation.voice-mute.command": "Command",
    "command.moderation.voice-mute.server": "Server",
}


@dataclass(frozen=True)
class ModifyRecord:
    model: type
    guild_snowflake: int
    channel_snowflake: int | None
    member_snowflake: int | None
    scope: str | None


class ModifyView(discord.ui.View):
    def __init__(
        self,
        author_snowflake: int,
        ctx: SnowflakeContext,
        interaction: discord.Interaction,
        modal: type,
        tick: Tick,
    ):
        super().__init__(timeout=120)
        self.__author_snowflake = author_snowflake
        self.__ctx = ctx
        self.__interaction = interaction
        self.__modal = modal
        self.__tick = tick
        self.available_channels = set()
        self.available_guilds = set()
        self.selected_model: type | None = None
        self.selected_guild: int | None = None
        self.selected_channel: int | None = None
        self.selected_scope: str | None = None
        self.records: list[ModifyRecord] = []

    async def setup(self):
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        guild = bot.get_guild(self.__ctx.guild_snowflake)
        if guild is None:
            raise GuildNotFound(str(self.__ctx.guild_snowflake))
        member = guild.get_member(self.__ctx.member_snowflake)
        if member is None:
            simplified_member = bot.registry.get(MemberState).active.get(
                self.__ctx.member_snowflake, None
            )
            if simplified_member is None:
                raise MemberNotFound(str(self.__ctx.member_snowflake))
        for guild in bot.guilds:
            for permission, model in PERMISSIONS.items():
                try:
                    database_factory: DatabaseFactory = DatabaseFactory(model)
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        channel_snowflake=self.__ctx.channel_snowflake,
                        guild_snowflake=self.__ctx.guild_snowflake,
                        member_snowflake=self.__author_snowflake,
                        requested=[permission],
                    )
                    items = await database_factory.select(
                        guild_snowflake=guild.id,
                        member_snowflake=self.__ctx.member_snowflake,
                        singular=False,
                    )
                    if items:
                        self.available_guilds.add(guild)
                        for item in items:
                            channel_snowflake = getattr(item, "channel_snowflake", None)
                            channel = (
                                guild.get_channel(channel_snowflake)
                                if channel_snowflake
                                else None
                            )
                            if channel:
                                self.available_channels.add(channel)
                            self.records.append(
                                ModifyRecord(
                                    model=model,
                                    guild_snowflake=guild.id,
                                    channel_snowflake=(channel.id if channel else None),
                                    member_snowflake=self.__ctx.member_snowflake,
                                    scope=None,
                                )
                            )
                except CheckFailure:
                    continue
            for permission, scope in SCOPES.items():
                try:
                    database_factory: DatabaseFactory = DatabaseFactory(VoiceMute)
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        channel_snowflake=self.__ctx.channel_snowflake,
                        guild_snowflake=self.__ctx.guild_snowflake,
                        member_snowflake=self.__author_snowflake,
                        requested=[permission],
                    )
                    items = await database_factory.select(
                        guild_snowflake=guild.id,
                        member_snowflake=self.__ctx.member_snowflake,
                        target=scope.lower(),
                        singular=False,
                    )
                    for item in items:
                        channel_snowflake = getattr(item, "channel_snowflake", None)
                        channel = (
                            guild.get_channel(channel_snowflake)
                            if channel_snowflake
                            else None
                        )
                        self.records.append(
                            ModifyRecord(
                                model=VoiceMute,
                                guild_snowflake=guild.id,
                                channel_snowflake=(channel.id if channel else None),
                                member_snowflake=self.__ctx.member_snowflake,
                                scope=scope,
                            )
                        )
                except CheckFailure:
                    continue
        if not self.records:
            raise CheckFailure("No records found.")
        guild_objs = self.limit_available_to_top_24_by_member_count(
            available=[
                g
                for g in (
                    bot.get_guild(gs)
                    for gs in {r.guild_snowflake for r in self.records}
                )
                if g is not None
            ]
        )
        channel_objs = self.limit_channels_to_top_24(
            available={
                c
                for c in (
                    bot.get_channel(cs)
                    for cs in {
                        r.channel_snowflake for r in self.records if r.channel_snowflake
                    }
                )
                if isinstance(
                    c, (discord.VoiceChannel, discord.TextChannel, discord.StageChannel)
                )
            }
        )
        allowed_guilds = {g.id for g in guild_objs}
        allowed_channels = {c.id for c in channel_objs}
        self.records = [
            r
            for r in self.records
            if r.guild_snowflake in allowed_guilds
            and (r.channel_snowflake is None or r.channel_snowflake in allowed_channels)
        ]
        self._update_disabled_state()
        self._refresh_options()
        if not self.guild_select.options or not self.channel_select.options:
            await self.__tick.end(
                warning="You have insufficient privileges to do that."
            )

    def _update_disabled_state(self) -> None:
        self.guild_select.disabled = self.selected_model is None
        guild_ready = self.selected_guild is not None
        self.channel_select.disabled = not guild_ready
        channel_ready = self.selected_channel is not None
        self.scope_select.disabled = (
            not channel_ready or self.selected_model is not VoiceMute
        )

    def limit_channels_to_top_24(self, available: set[discord.abc.GuildChannel]):
        relevant = [
            c
            for c in available
            if isinstance(
                c, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
            )
        ]
        relevant.sort(
            key=lambda c: (
                len(c.members)
                if isinstance(c, (discord.VoiceChannel, discord.StageChannel))
                else 0
            ),
            reverse=True,
        )
        return relevant[:24]

    async def interaction_check(self, interaction):
        return interaction.user.id == self.__author_snowflake

    def _label_for(self, dim: str) -> str:
        value = getattr(self, f"selected_{dim}")
        if value is None:
            return {
                "model": "Select a model",
                "guild": "Select a server",
                "channel": "Select a channel",
                "scope": "Select a scope",
            }[dim]
        if dim == "model":
            return value.name
        if dim == "guild":
            return self._guild_label(value)
        if dim == "channel":
            return self._channel_label(value)
        if dim == "scope":
            return value
        return "Unknown"

    def _visible_records(self, *, exclude: str | None) -> list[ModifyRecord]:
        records = self.records
        if exclude != "model" and self.selected_model is not None:
            records = [r for r in records if r.model == self.selected_model]
        if exclude != "guild" and self.selected_guild is not None:
            records = [r for r in records if r.guild_snowflake == self.selected_guild]
        if exclude != "channel" and self.selected_channel is not None:
            records = [
                r
                for r in records
                if r.channel_snowflake in (None, self.selected_channel)
            ]
        if exclude != "scope" and self.selected_scope is not None:
            records = [r for r in records if r.scope in (None, self.selected_scope)]
        records = [
            r for r in records if r.member_snowflake == self.__ctx.member_snowflake
        ]
        return records

    def limit_available_to_top_24_by_member_count(self, available):
        items = list(available)
        items.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
        top_24 = items[:24]
        return set(top_24)

    def _refresh_options(self) -> None:
        self._build_model_options(
            models={r.model for r in self._visible_records(exclude="model") if r.model}
        )
        self.model_select.placeholder = self._label_for("model")
        self._build_guild_options(
            guild_snowflakes={
                r.guild_snowflake for r in self._visible_records(exclude="guild")
            }
        )
        self.guild_select.placeholder = self._label_for("guild")
        self._build_channel_options(
            channel_snowflakes={
                r.channel_snowflake
                for r in self._visible_records(exclude="channel")
                if r.channel_snowflake
            }
        )
        self.channel_select.placeholder = self._label_for("channel")
        self._build_scope_options(
            scopes={r.scope for r in self._visible_records(exclude="scope") if r.scope}
        )
        self.scope_select.placeholder = self._label_for("scope")

    def _build_model_options(self, models):
        model_options = []
        model_options.extend(
            [discord.SelectOption(label=m.name, value=m.identifier) for m in models]
        )
        self.model_select.options = model_options

    def _guild_label(self, guild_snowflake: int):
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild:
            return guild.name
        else:
            return "Error"

    def _build_guild_options(self, guild_snowflakes: set[int], default: bool = False):
        bot = DiscordBot.get_instance()
        guilds = [
            g
            for g in (bot.get_guild(snowflake) for snowflake in guild_snowflakes)
            if g is not None
        ]
        guilds.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
        top_24 = guilds[:24]
        guild_options = []
        guild_options.extend(
            [
                discord.SelectOption(
                    label=self._guild_label(guild.id), value=str(guild.id)
                )
                for guild in top_24
            ]
        )
        self.guild_select.options = guild_options

    def _channel_label(self, channel_snowflake: int):
        bot: DiscordBot = DiscordBot.get_instance()
        channel = bot.get_channel(channel_snowflake)
        if isinstance(channel, discord.abc.GuildChannel):
            return channel.name
        else:
            return "Error"

    def _build_channel_options(
        self, channel_snowflakes: set[int], default: bool = False
    ):
        bot = DiscordBot.get_instance()
        channels = [
            c
            for c in (bot.get_channel(snowflake) for snowflake in channel_snowflakes)
            if isinstance(
                c,
                (discord.TextChannel, discord.VoiceChannel, discord.StageChannel),
            )
        ]
        channels.sort(
            key=lambda c: (
                len(c.members)
                if isinstance(c, (discord.VoiceChannel, discord.StageChannel))
                else 0
            ),
            reverse=True,
        )
        self.channel_select.options = [
            *[discord.SelectOption(label=c.name, value=str(c.id)) for c in channels],
        ]

    def _build_scope_options(self, scopes: set[str]):
        scope_options = []
        scope_options.extend(
            [
                discord.SelectOption(
                    label=s, value=s, default=(self.selected_scope == s)
                )
                for s in scopes
            ]
        )
        if scope_options:
            self.scope_select.options = scope_options

    def _get_model_by_identifier(self, model_identifier: str):
        for model in MODELS:
            if model.identifier == model_identifier:
                return model
        else:
            return None

    @discord.ui.select(placeholder="Select a model", options=[])
    async def model_select(self, interaction, select):
        raw = select.values[0]
        self.selected_model = self._get_model_by_identifier(raw)
        self.selected_guild = None
        self.selected_channel = None
        self.selected_scope = None
        self._update_disabled_state()
        self._refresh_options()
        await interaction.response.edit_message(view=self)

    @discord.ui.select(placeholder="Select a server", options=[], disabled=True)
    async def guild_select(self, interaction, select):
        raw = select.values[0]
        self.selected_guild = int(raw)
        self.selected_channel = None
        self.selected_scope = None
        self.scope_select.disabled = True
        self._update_disabled_state()
        self._refresh_options()
        await interaction.response.edit_message(view=self)

    @discord.ui.select(placeholder="Select a channel", options=[], disabled=True)
    async def channel_select(self, interaction, select):
        raw = select.values[0]
        self.selected_channel = int(raw)
        self.selected_scope = None
        self.scope_select.disabled = self.selected_model is not VoiceMute
        self._update_disabled_state()
        self._refresh_options()
        await interaction.response.edit_message(view=self)

    @discord.ui.select(
        placeholder="Select a scope",
        options=[
            discord.SelectOption(label="Only valid for voice mute.", value="default")
        ],
        disabled=True,
    )
    async def scope_select(self, interaction, select):
        raw = select.values[0]
        self.selected_scope = raw
        self._refresh_options()
        await interaction.response.edit_message(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        if self.selected_model is None:
            return await interaction.response.send_message(
                content="Please select a model.", ephemeral=True
            )
        if self.guild_select.disabled is False and self.selected_guild is None:
            return await interaction.response.send_message(
                content="Please select a server.", ephemeral=True
            )
        if self.channel_select.disabled is False and self.selected_channel is None:
            return await interaction.response.send_message(
                content="Please select a channel.", ephemeral=True
            )
        if self.scope_select.disabled is False and self.selected_scope is None:
            return await interaction.response.send_message(
                content="Please select a scope.", ephemeral=True
            )
        self.stop()
        database_factory: DatabaseFactory = DatabaseFactory(self.selected_model)
        select_kwargs: dict[str, Any] = {
            "guild_snowflake": self.selected_guild,
            "member_snowflake": self.__ctx.member_snowflake
        }
        if self.selected_channel is not None:
            select_kwargs["channel_snowflake"] = self.selected_channel
        if self.selected_scope is not None:
            select_kwargs["target"] = self.selected_scope 
        bot: DiscordBot = DiscordBot.get_instance()
        record = await database_factory.select(singular=True, **select_kwargs)
        if self.__modal == ReasonModal:
            modal = ReasonModal(
                author_snowflake=self.__author_snowflake,
                model=self.selected_model,
                record=record,
                tick=self.__tick,
            )
            await modal.setup()
            await interaction.response.send_modal(modal)
        elif self.__modal == DurationModal:
            modal = DurationModal(
                author_snowflake=self.__author_snowflake,
                model=self.selected_model,
                record=record,
                tick=self.__tick,
            )
            await modal.setup()
            await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        await interaction.response.defer()
        await interaction.delete_original_response()
        self.stop()

    async def on_error(
        self,
        interaction: discord.Interaction,
        error: Exception,
        item: discord.ui.Item,
    ) -> None:
        await self.__tick.end(error=str(error), ephemeral=True)
