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

from vyrtuous.aliases import alias_service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.db.autoassign import AutoAssignRole
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.permissions import permission_service
from vyrtuous.utils.channels import automute_channel_service, video_channel_service
from vyrtuous.utils.errors.error import CheckFailure, MemberNotFound
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    cap_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.tracking import stream_service
from vyrtuous.utils.users import autoassign_role_service, vegan_service
from vyrtuous.view.view_context import ViewContext

MODELS = {
    alias_service.MODEL,
    autoassign_role_service.MODEL,
    automute_channel_service.MODEL,
    ban_service.MODEL,
    cap_service.MODEL,
    flag_service.MODEL,
    permission_service.MODEL,
    text_mute_service.MODEL,
    vegan_service.MODEL,
    voice_mute_service.MODEL,
    stream_service.MODEL,
    video_channel_service.MODEL,
}

GUILD_PERMISSIONS = {
    "command.clear.category.alias": alias_service.MODEL,
    "command.clear.category.autoassign": autoassign_role_service.MODEL,
    "command.clear.category.automute": automute_channel_service.MODEL,
    "command.clear.category.ban": ban_service.MODEL,
    "command.clear.category.cap": cap_service.MODEL,
    "command.clear.category.flag": flag_service.MODEL,
    "command.clear.category.groups": permission_service.MODEL,
    "command.clear.category.text-mute": text_mute_service.MODEL,
    "command.clear.category.vegan": vegan_service.MODEL,
    "command.clear.category.voice-mute": voice_mute_service.MODEL,
    "command.clear.category.stream": stream_service.MODEL,
    "command.clear.category.video-channel": video_channel_service.MODEL,
}

CHANNEL_PERMISSIONS = {
    "command.clear.category.alias": alias_service.MODEL,
    "command.clear.category.autoassign": autoassign_role_service.MODEL,
    "command.clear.category.automute": automute_channel_service.MODEL,
    "command.clear.category.ban": ban_service.MODEL,
    "command.clear.category.cap": cap_service.MODEL,
    "command.clear.category.flag": flag_service.MODEL,
    "command.clear.category.groups": permission_service.MODEL,
    "command.clear.category.text-mute": text_mute_service.MODEL,
    "command.clear.category.voice-mute": voice_mute_service.MODEL,
    "command.clear.category.stream": stream_service.MODEL,
    "command.clear.category.video-channel": video_channel_service.MODEL,
}
MEMBER_PERMISSIONS = {
    "command.clear.category.ban": ban_service.MODEL,
    "command.clear.category.flag": flag_service.MODEL,
    "command.clear.category.groups": permission_service.MODEL,
    "command.clear.category.text-mute": text_mute_service.MODEL,
    "command.clear.category.vegan": vegan_service.MODEL,
    "command.clear.category.voice-mute": voice_mute_service.MODEL,
}

SCOPES = {
    "command.clear.category.voice-mute.auto": "Auto",
    "command.clear.category.voice-mute.click": "Click",
    "command.clear.category.voice-mute.command": "Command",
    "command.clear.category.voice-mute.server": "Server",
}


@dataclass(frozen=True)
class ClearRecord:
    model: type | None
    guild_snowflake: int
    channel_snowflake: int | None
    scope: str | None


class _All:
    def __repr__(self) -> str:
        return "ALL"


ALL = _All()


class ClearView(discord.ui.View):
    def __init__(
        self,
        author_snowflake: int,
        ctx: ViewContext,
        obj: Union[str, int, discord.Guild, discord.Member, discord.abc.GuildChannel],
        tick: Tick,
    ):
        super().__init__(timeout=120)
        self.__author_snowflake = author_snowflake
        self.__ctx = ctx
        self.__obj = obj
        self.__tick = tick
        self.available_channels = set()
        self.available_guilds = set()
        self.selected_model: type | _All | None = None
        self.selected_guild: int | _All | None = None
        self.selected_channel: int | _All | None = None
        self.selected_scope: str | _All | None = None
        self.records: list[ClearRecord] = []

    async def setup(self):
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        match self.__obj:
            case str() as all:
                if all.lower() != "all":
                    raise CheckFailure(
                        f"Invalid option `{all}`. Must be a `channel`, `member`, `server` or `all`."
                    )
                permissions: dict[str, Any] = {}
                permissions.update(CHANNEL_PERMISSIONS)
                permissions.update(GUILD_PERMISSIONS)
                permissions.update(MEMBER_PERMISSIONS)
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    channel_snowflake=self.__ctx.channel_snowflake,
                    guild_snowflake=self.__ctx.guild_snowflake,
                    member_snowflake=self.__author_snowflake,
                    requested=[
                        "command.clear.scope.all",
                        "other_guilds",
                    ],
                )
                for guild in bot.guilds:
                    for permission, model in permissions.items():
                        try:
                            database_factory: DatabaseFactory = DatabaseFactory(model)
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=self.__ctx.channel_snowflake,
                                guild_snowflake=guild.id,
                                member_snowflake=self.__author_snowflake,
                                requested=[permission],
                            )
                            items = await database_factory.select(
                                guild_snowflake=guild.id, singular=False
                            )
                            if items:
                                self.available_guilds.add(guild)
                                for item in items:
                                    channel_snowflake = getattr(
                                        item, "channel_snowflake", None
                                    )
                                    channel = (
                                        guild.get_channel(channel_snowflake)
                                        if channel_snowflake
                                        else None
                                    )
                                    if channel:
                                        self.available_channels.add(channel)
                                    self.records.append(
                                        ClearRecord(
                                            model=model,
                                            guild_snowflake=guild.id,
                                            channel_snowflake=(
                                                channel.id if channel else None
                                            ),
                                            scope=None,
                                        )
                                    )
                        except CheckFailure:
                            continue
                    for permission, scope in SCOPES.items():
                        try:
                            database_factory: DatabaseFactory = DatabaseFactory(
                                VoiceMute
                            )
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=self.__ctx.channel_snowflake,
                                guild_snowflake=self.__ctx.guild_snowflake,
                                member_snowflake=self.__author_snowflake,
                                requested=[permission],
                            )
                            items = await database_factory.select(
                                guild_snowflake=guild.id,
                                target=scope.lower(),
                                singular=False,
                            )
                            for item in items:
                                channel_snowflake = getattr(
                                    item, "channel_snowflake", None
                                )
                                channel = (
                                    guild.get_channel(channel_snowflake)
                                    if channel_snowflake
                                    else None
                                )
                                self.records.append(
                                    ClearRecord(
                                        model=VoiceMute,
                                        guild_snowflake=guild.id,
                                        channel_snowflake=(
                                            channel.id if channel else None
                                        ),
                                        scope=scope,
                                    )
                                )
                        except CheckFailure:
                            continue
            case discord.Guild() as guild:
                self.remove_item(self.channel_select)
                self.remove_item(self.guild_select)
                self.available_guilds.add(guild)
                self.available_channels = (
                    self.limit_available_to_top_24_by_member_count(
                        available=[channel for channel in guild.channels]
                    )
                )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    channel_snowflake=self.__ctx.channel_snowflake,
                    guild_snowflake=guild.id,
                    member_snowflake=self.__author_snowflake,
                    requested=[
                        "command.clear.scope.guild",
                    ],
                )
                for permission, model in GUILD_PERMISSIONS.items():
                    try:
                        database_factory: DatabaseFactory = DatabaseFactory(model)
                        await permission_service.has_permissions(
                            permission_state=permission_state,
                            channel_snowflake=self.__ctx.channel_snowflake,
                            guild_snowflake=guild.id,
                            member_snowflake=self.__author_snowflake,
                            requested=[permission],
                        )
                        items = await database_factory.select(
                            guild_snowflake=guild.id, singular=False
                        )
                        if items:
                            for item in items:
                                channel_snowflake = getattr(
                                    item, "channel_snowflake", None
                                )
                                channel = (
                                    guild.get_channel(channel_snowflake)
                                    if channel_snowflake
                                    else None
                                )
                                if channel:
                                    self.available_channels.add(channel)
                                self.records.append(
                                    ClearRecord(
                                        model=model,
                                        guild_snowflake=guild.id,
                                        channel_snowflake=(
                                            channel.id if channel else None
                                        ),
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
                            guild_snowflake=guild.id,
                            member_snowflake=self.__author_snowflake,
                            requested=[permission],
                        )
                        items = await database_factory.select(
                            guild_snowflake=guild.id,
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
                                ClearRecord(
                                    model=VoiceMute,
                                    guild_snowflake=guild.id,
                                    channel_snowflake=channel.id if channel else None,
                                    scope=scope,
                                )
                            )
                    except CheckFailure:
                        continue
                self.available_channels = self.available_channels | (
                    self.limit_available_to_top_24_by_member_count(
                        available=[channel for channel in guild.channels]
                    )
                )
            case discord.Member() as member:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    channel_snowflake=self.__ctx.channel_snowflake,
                    guild_snowflake=self.__ctx.guild_snowflake,
                    member_snowflake=self.__author_snowflake,
                    requested=[
                        "command.clear.scope.member",
                    ],
                )
                for guild in bot.guilds:
                    for permission, model in MEMBER_PERMISSIONS.items():
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
                                member_snowflake=member.id,
                                singular=False,
                            )
                            if items:
                                self.available_guilds.add(guild)
                                for item in items:
                                    channel_snowflake = getattr(
                                        item, "channel_snowflake", None
                                    )
                                    channel = (
                                        guild.get_channel(channel_snowflake)
                                        if channel_snowflake
                                        else None
                                    )
                                    if channel:
                                        self.available_channels.add(channel)
                                    self.records.append(
                                        ClearRecord(
                                            model=model,
                                            guild_snowflake=guild.id,
                                            channel_snowflake=(
                                                channel.id if channel else None
                                            ),
                                            scope=None,
                                        )
                                    )
                        except CheckFailure:
                            continue
                    for permission, scope in SCOPES.items():
                        try:
                            database_factory: DatabaseFactory = DatabaseFactory(
                                VoiceMute
                            )
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=self.__ctx.channel_snowflake,
                                guild_snowflake=self.__ctx.guild_snowflake,
                                member_snowflake=self.__author_snowflake,
                                requested=[permission],
                            )
                            items = await database_factory.select(
                                member_snowflake=member.id,
                                target=scope.lower(),
                                singular=False,
                            )
                            for item in items:
                                channel_snowflake = getattr(
                                    item, "channel_snowflake", None
                                )
                                channel = (
                                    guild.get_channel(channel_snowflake)
                                    if channel_snowflake
                                    else None
                                )
                                self.records.append(
                                    ClearRecord(
                                        model=VoiceMute,
                                        guild_snowflake=guild.id,
                                        channel_snowflake=(
                                            channel.id if channel else None
                                        ),
                                        scope=scope,
                                    )
                                )
                        except CheckFailure:
                            continue
            case int() as member_snowflake:
                simplified_member = bot.registry.get(MemberState).active.get(
                    member_snowflake, None
                )
                if simplified_member is None:
                    raise MemberNotFound(str(member_snowflake))
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    channel_snowflake=self.__ctx.channel_snowflake,
                    guild_snowflake=self.__ctx.guild_snowflake,
                    member_snowflake=self.__author_snowflake,
                    requested=[
                        "command.clear.scope.member",
                    ],
                )
                for guild in bot.guilds:
                    for permission, model in MEMBER_PERMISSIONS.items():
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
                                member_snowflake=member_snowflake,
                                singular=False,
                            )
                            if items:
                                self.available_guilds.add(guild)
                                for item in items:
                                    channel_snowflake = getattr(
                                        item, "channel_snowflake", None
                                    )
                                    channel = (
                                        guild.get_channel(channel_snowflake)
                                        if channel_snowflake
                                        else None
                                    )
                                    if channel:
                                        self.available_channels.add(channel)
                                    self.records.append(
                                        ClearRecord(
                                            model=model,
                                            guild_snowflake=guild.id,
                                            channel_snowflake=(
                                                channel.id if channel else None
                                            ),
                                            scope=None,
                                        )
                                    )
                        except CheckFailure:
                            continue
                    for permission, scope in SCOPES.items():
                        try:
                            database_factory: DatabaseFactory = DatabaseFactory(
                                VoiceMute
                            )
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=self.__ctx.channel_snowflake,
                                guild_snowflake=self.__ctx.guild_snowflake,
                                member_snowflake=self.__author_snowflake,
                                requested=[permission],
                            )
                            items = await database_factory.select(
                                member_snowflake=member_snowflake,
                                target=scope.lower(),
                                singular=False,
                            )
                            for item in items:
                                channel_snowflake = getattr(
                                    item, "channel_snowflake", None
                                )
                                channel = (
                                    guild.get_channel(channel_snowflake)
                                    if channel_snowflake
                                    else None
                                )
                                self.records.append(
                                    ClearRecord(
                                        model=VoiceMute,
                                        guild_snowflake=guild.id,
                                        channel_snowflake=(
                                            channel.id if channel else None
                                        ),
                                        scope=scope,
                                    )
                                )
                        except CheckFailure:
                            continue
            case discord.abc.GuildChannel() as channel:
                self.remove_item(self.guild_select)
                self.remove_item(self.channel_select)
                self.available_channels.add(channel)
                self.available_guilds.add(channel.guild)
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    channel_snowflake=channel.id,
                    guild_snowflake=channel.guild.id,
                    member_snowflake=self.__author_snowflake,
                    requested=[
                        "command.clear.scope.channel",
                    ],
                )
                for permission, model in CHANNEL_PERMISSIONS.items():
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
                            channel_snowflake=channel.id,
                            guild_snowflake=channel.guild.id,
                            singular=False,
                        )
                        if items:
                            self.records.append(
                                ClearRecord(
                                    model=model,
                                    guild_snowflake=channel.guild.id,
                                    channel_snowflake=channel.id,
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
                            channel_snowflake=channel.id,
                            guild_snowflake=channel.guild.id,
                            target=scope.lower(),
                            singular=False,
                        )
                        if items:
                            self.records.append(
                                ClearRecord(
                                    model=VoiceMute,
                                    guild_snowflake=channel.guild.id,
                                    channel_snowflake=channel.id,
                                    scope=scope,
                                )
                            )
                    except CheckFailure:
                        continue
        if not self.records:
            raise CheckFailure("No records to clear found.")
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
        channel_objs = self.limit_available_to_top_24_by_member_count(
            available=[
                c
                for c in (
                    bot.get_channel(cs)
                    for cs in {
                        r.channel_snowflake for r in self.records if r.channel_snowflake
                    }
                )
                if c is not None
            ]
        )
        allowed_guilds = {g.id for g in guild_objs}
        allowed_channels = {c.id for c in channel_objs}
        self.records = [
            r
            for r in self.records
            if r.guild_snowflake in allowed_guilds
            and (r.channel_snowflake is None or r.channel_snowflake in allowed_channels)
        ]
        self._refresh_options()

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
        if value is ALL:
            return "All"
        if dim == "model":
            return value.name
        if dim == "guild":
            return self._guild_label(value)
        if dim == "channel":
            return self._channel_label(value)
        if dim == "scope":
            return value
        return "All"

    def _visible_records(self, *, exclude: str) -> list[ClearRecord]:
        records = self.records
        if exclude != "model" and self.selected_model not in (None, ALL):
            records = [r for r in records if r.model == self.selected_model]
        if exclude != "guild" and self.selected_guild not in (None, ALL):
            records = [r for r in records if r.guild_snowflake == self.selected_guild]
        if exclude != "channel" and self.selected_channel not in (None, ALL):
            records = [
                r
                for r in records
                if r.channel_snowflake in (None, self.selected_channel)
            ]
        if exclude != "scope" and self.selected_scope not in (None, ALL):
            records = [r for r in records if r.scope in (None, self.selected_scope)]
        return records

    def limit_available_to_top_24_by_member_count(self, available):
        items = set()
        items = available
        items.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
        top_24 = items[:24]
        return set(top_24)

    def _refresh_options(self) -> None:
        self._build_model_options(
            models={r.model for r in self._visible_records(exclude="model") if r.model}
        )
        self._build_guild_options(
            guild_snowflakes={
                r.guild_snowflake for r in self._visible_records(exclude="guild")
            }
        )
        self._build_channel_options(
            channel_snowflakes={
                r.channel_snowflake
                for r in self._visible_records(exclude="channel")
                if r.channel_snowflake
            }
        )
        self._build_scope_options(
            scopes={r.scope for r in self._visible_records(exclude="scope") if r.scope}
        )
        self.model_select.placeholder = self._label_for("model")
        self.guild_select.placeholder = self._label_for("guild")
        self.channel_select.placeholder = self._label_for("channel")
        self.scope_select.placeholder = self._label_for("scope")

    def _build_model_options(self, models, default: bool = False):
        model_options = []
        model_options.append(
            discord.SelectOption(
                label="All",
                value="all",
                default=default,
            )
        )
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
        guild_options = []
        guild_options.append(
            discord.SelectOption(
                label="All",
                value="all",
                default=default,
            )
        )
        guild_options.extend(
            [
                discord.SelectOption(
                    label=self._guild_label(guild_snowflake), value=str(guild_snowflake)
                )
                for guild_snowflake in guild_snowflakes
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
        channel_options = []
        channel_options.append(
            discord.SelectOption(
                label="All",
                value="all",
                default=default,
            )
        )
        channel_options.extend(
            [
                discord.SelectOption(
                    label=self._channel_label(channel_snowflake),
                    value=str(channel_snowflake),
                )
                for channel_snowflake in channel_snowflakes
            ]
        )
        self.channel_select.options = channel_options

    def _build_scope_options(self, scopes: set[str], default: bool = False):
        scope_options = []
        scope_options.append(
            discord.SelectOption(
                label="All",
                value="all",
                default=default,
            )
        )
        scope_options.extend([discord.SelectOption(label=s, value=s) for s in scopes])
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
        self.selected_model = (
            ALL if raw == "all" else self._get_model_by_identifier(raw)
        )
        self.selected_guild = None
        self.selected_channel = None
        self.selected_scope = None
        self.guild_select.disabled = False
        self.channel_select.disabled = True
        self.scope_select.disabled = True
        self._refresh_options()
        await interaction.response.edit_message(view=self)

    @discord.ui.select(placeholder="Select a server", options=[], disabled=True)
    async def guild_select(self, interaction, select):
        raw = select.values[0]
        self.selected_guild = ALL if raw == "all" else int(raw)
        self.selected_channel = None
        self.selected_scope = None
        self.channel_select.disabled = self.selected_model is AutoAssignRole
        self.scope_select.disabled = True
        self._refresh_options()
        await interaction.response.edit_message(view=self)

    @discord.ui.select(placeholder="Select a channel", options=[], disabled=True)
    async def channel_select(self, interaction, select):
        raw = select.values[0]
        self.selected_channel = ALL if raw == "all" else int(raw)
        self.selected_scope = None
        self.scope_select.disabled = (
            self.selected_model is not VoiceMute or self.selected_model is not ALL
        )
        self._refresh_options()
        await interaction.response.edit_message(view=self)

    @discord.ui.select(placeholder="Select a scope", options=[], disabled=True)
    async def scope_select(self, interaction, select):
        raw = select.values[0]
        self.selected_scope = ALL if raw == "all" else raw
        self._refresh_options()
        await interaction.response.edit_message(view=self)

    # @discord.ui.button(label="Confirm", style=discord.ButtonStyle.green)
    # async def confirm(self, interaction, button):
    #     self._result = True
    #     await interaction.message.delete()
    #     self.stop()
    #
    # @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    # async def cancel(self, interaction, button):
    #     self._result = False
    #     await interaction.message.delete()
    #     self.stop()
    #
    # def build_embed(self):
    #     embed = discord.Embed(
    #         title="\U000026a0\U0000fe0f Clear Command Confirmation",
    #         description=f"**Action:** {self.__action}\n**Target:** {self.__mention}",
    #         color=discord.Color.orange(),
    #     )
    #     embed.set_footer(text="Please confirm or cancel this action.")
    #     return embed
    #
    # @property
    # def channel_snowflake(self):
    #     return self._channel_snowflake
    #
    # @property
    # def guild_snowflake(self):
    #     return self._guild_snowflake
    #
    # @property
    # def member_snowflake(self):
    #     return self._member_snowflake
    #
    # @property
    # def result(self):
    #     return self._result
