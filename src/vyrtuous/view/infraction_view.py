"""!/bin/python3
infraction_view.py The purpose of this program is to provide the view for creating an infraction.

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

import re

import discord

from vyrtuous.aliases import (
    unban_alias_service,
    unflag_alias_service,
    untext_mute_alias_service,
    unvoice_mute_alias_service,
)
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.permissions import PermissionScope
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.ban import Ban
from vyrtuous.db.cap import Cap
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.flag import Flag
from vyrtuous.db.text_mute import TextMute
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.modal.infraction_model import InfractionModal
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import CheckFailure, GuildNotFound
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import cap_service
from vyrtuous.utils.tracking import stream_service

DURATIONS = {
    "1hour",
    "8hour",
    "1day",
    "1week",
}


class InfractionView(discord.ui.View):
    def __init__(
        self,
        author_snowflake: int,
        ctx: SnowflakeContext,
        interaction: discord.Interaction,
        model: type,
        tick: Tick,
    ):
        super().__init__(timeout=120)
        self.available_channels = set()
        self.available_guilds = set()
        self.__all_available_channels = set()
        self.__author_snowflake = author_snowflake
        self.__channel_snowflake = ctx.channel_snowflake
        self.__ctx = ctx
        self.__duration = None
        self.__guild_snowflake = ctx.guild_snowflake
        self.__interaction = interaction
        self.__tick = tick
        self.__model = model
        if self.__model == Flag:
            self.remove_item(self.duration_select)

    async def interaction_check(self, interaction):
        return interaction.user.id == self.__author_snowflake

    async def setup(self):
        available_channels: set[discord.abc.GuildChannel] = set()
        available_guilds: set[discord.Guild] = set()
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        all_assigned_groups = await permission_service.resolve_all_assigned_groups(
            permission_state=permission_state, member_snowflake=self.__author_snowflake
        )
        for group, assigned_guild_snowflake, _ in all_assigned_groups:
            if group.scope == PermissionScope.GLOBAL:
                for guild in bot.guilds:
                    for channel in guild.channels:
                        available_channels.add(channel)
                    available_guilds.add(guild)
                try:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=self.__author_snowflake,
                        requested=["command.moderation.uncapped"],
                    )
                except:
                    continue
            elif group.scope == PermissionScope.GUILD:
                for guild in bot.guilds:
                    if guild.id != assigned_guild_snowflake:
                        try:
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                member_snowflake=self.__author_snowflake,
                                requested=["other_guilds"],
                            )
                        except:
                            continue
                    for channel in guild.channels:
                        available_channels.add(channel)
                    available_guilds.add(guild)
                    try:
                        await permission_service.has_permissions(
                            permission_state=permission_state,
                            guild_snowflake=guild.id,
                            member_snowflake=self.__author_snowflake,
                            requested=["command.moderation.uncapped"],
                        )
                    except:
                        continue
            elif group.scope == PermissionScope.CHANNEL:
                for guild in bot.guilds:
                    for channel in guild.channels:
                        effective_group = (
                            await permission_service.resolve_effective_group(
                                permission_state=permission_state,
                                member_snowflake=self.__author_snowflake,
                                channel_snowflake=channel.id,
                                guild_snowflake=guild.id,
                            )
                        )
                        if effective_group == group:
                            available_channels.add(channel)
                            available_guilds.add(guild)
                        try:
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=channel.id,
                                guild_snowflake=guild.id,
                                member_snowflake=self.__author_snowflake,
                                requested=["command.moderation.uncapped"],
                            )
                        except:
                            continue
            else:
                raise CheckFailure(
                    "You do not have sufficient privileges in this channel or server to use this command."
                )
        self.__all_available_channels = available_channels
        self.available_guilds = available_guilds
        self.available_channels = [
            channel
            for channel in available_channels
            if channel.guild.id == self.__ctx.guild_snowflake
        ]
        self._build_guild_options(available_guilds=available_guilds, default=True)
        limited_channels = self.limit_channels_to_top_24(
            available=set(self.available_channels)
        )
        database_factory: DatabaseFactory = DatabaseFactory(self.__model)
        if self.__model == VoiceMute:
            types = (discord.VoiceChannel, discord.StageChannel)
            records = []
            auto_records = await database_factory.select(
                channel_snowflake=self.__channel_snowflake,
                guild_snowflake=self.__guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                target="auto",
                singular=True,
            )
            command_records = await database_factory.select(
                channel_snowflake=self.__channel_snowflake,
                guild_snowflake=self.__guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                target="command",
                singular=True,
            )
            records.append(auto_records)
            records.append(command_records)
        else:
            types = (discord.VoiceChannel, discord.StageChannel, discord.TextChannel)
            records = await database_factory.select(
                channel_snowflake=self.__channel_snowflake,
                guild_snowflake=self.__guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                singular=True,
            )
        self._build_channel_options(
            limited_channels=limited_channels, default=True, types=types
        )
        await self.build_duration_options()
        if records:
            self.remove_item(self.duration_select)
        if (
            not self.guild_select.options
            or not self.channel_select.options
            or not self.duration_select.options
        ):
            await self.__tick.end(
                warning="You have insufficient privileges to do that."
            )

    async def _sync_duration_select(self, rebuild_options: bool = True):
        database_factory: DatabaseFactory = DatabaseFactory(self.__model)
        records = await database_factory.select(
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
            member_snowflake=self.__ctx.member_snowflake,
            singular=False,
        )
        should_show = self.__model != Flag and not records
        is_shown = self.duration_select in self.children
        if should_show:
            if not is_shown:
                self.add_item(self.duration_select)
            if rebuild_options:
                await self.build_duration_options()
            self.duration_select.disabled = False
        elif is_shown:
            self.remove_item(self.duration_select)

    async def build_duration_options(self, available_durations: set[str] | None = None):
        def _format_duration_label(duration: str):
            if duration == "0":
                return "Permanent"
            match = re.fullmatch(r"(\d+)([a-zA-Z]+)", duration)
            if match is None:
                return duration
            amount, unit = match.groups()
            return f"{amount} {unit}"

        if available_durations is None:
            bot: DiscordBot = DiscordBot.get_instance()
            cap = None
            duration_builder: DurationBuilder = DurationBuilder()
            database_factory: DatabaseFactory = DatabaseFactory(Cap)
            permission_state: PermissionState = bot.registry.get(PermissionState)
            try:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__guild_snowflake,
                    member_snowflake=self.__author_snowflake,
                    requested=["command.moderation.uncapped"],
                )
            except CheckFailure:
                cap = await database_factory.select(
                    category=self.__model.identifier,
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__guild_snowflake,
                    singular=True,
                )
                if cap:
                    available_durations = self.limit_available_to_cap(
                        duration_seconds=cap.duration_seconds
                    )
                else:
                    duration_seconds = duration_builder.load(
                        DurationObject(number=8, prefix="", sign=1, unit="")
                    ).to_seconds()
                    available_durations = self.limit_available_to_cap(
                        duration_seconds=duration_seconds
                    )
                duration_options = [
                    discord.SelectOption(
                        label=_format_duration_label(d),
                        value=d,
                    )
                    for d in available_durations
                ]
                self.duration_select.options = duration_options
            else:
                durations = set(DURATIONS)
                durations.add("0")
                duration_options = [
                    discord.SelectOption(
                        label=_format_duration_label(d),
                        value=d,
                    )
                    for d in durations
                ]
                self.duration_select.options = duration_options
        else:
            duration_options = [
                discord.SelectOption(
                    label=_format_duration_label(d),
                    value=d,
                )
                for d in available_durations
            ]
            self.duration_select.options = duration_options

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

    def limit_available_to_cap(self, duration_seconds: int):
        duration_builder: DurationBuilder = DurationBuilder()
        duration = duration_builder.from_seconds(duration_seconds).to_timedelta()
        options = set()
        durations = ["1hour", "8hours", "1day", "1week"]
        for cmp in durations:
            compare_duration = duration_builder.parse(cmp).to_timedelta()
            if duration >= compare_duration:
                options.add(cmp)
        return options

    def limit_available_to_top_24_by_member_count(self, available):
        items = list(available)
        items.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
        top_24 = items[:24]
        return set(top_24)

    def _build_channel_options(
        self,
        limited_channels: list[
            discord.TextChannel | discord.VoiceChannel | discord.StageChannel
        ],
        types: tuple[
            type[discord.VoiceChannel | discord.StageChannel | discord.TextChannel], ...
        ],
        default: bool = False,
    ):
        channel_options = []
        bot: DiscordBot = DiscordBot.get_instance()
        channel = bot.get_channel(self.__ctx.channel_snowflake)
        if (
            channel
            and isinstance(channel, types)
            and channel.guild.id == self.__guild_snowflake
            and channel in self.available_channels
        ):
            channel_options.append(
                discord.SelectOption(
                    label=channel.name,
                    value=str(self.__ctx.channel_snowflake),
                    default=default,
                )
            )
        channel_options.extend(
            [
                discord.SelectOption(label=c.name, value=str(c.id))
                for c in limited_channels
                if isinstance(c, types) and c.id != self.__ctx.channel_snowflake
            ]
        )
        self.channel_select.options = list(channel_options)

    def _build_guild_options(
        self, available_guilds: set[discord.Guild], default: bool = False
    ):
        guild_options = []
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(self.__ctx.guild_snowflake)
        if guild:
            guild_options.append(
                discord.SelectOption(
                    label=guild.name,
                    value=str(self.__ctx.guild_snowflake),
                    default=default,
                )
            )
        limited_guilds = self.limit_available_to_top_24_by_member_count(
            available=available_guilds
        )
        guild_options.extend(
            [
                discord.SelectOption(label=g.name, value=str(g.id))
                for g in limited_guilds
                if g.id != self.__ctx.guild_snowflake
            ]
        )
        self.guild_select.options = guild_options

    @discord.ui.select(
        placeholder="Select a server",
        options=[],
    )
    async def guild_select(self, interaction, select):
        await interaction.response.defer()
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(int(select.values[0]))
        if guild is None:
            raise GuildNotFound(str(select.values[0]))
        for option in self.guild_select.options:
            option.default = False
        self.__guild_snowflake = guild.id
        self.guild_select.placeholder = guild.name
        available_channels = {
            channel
            for channel in self.__all_available_channels
            if channel in guild.channels
        }
        limited_channels = self.limit_channels_to_top_24(
            available=available_channels,
        )
        if self.__model == VoiceMute:
            types = (discord.VoiceChannel, discord.StageChannel)
        else:
            types = (discord.VoiceChannel, discord.StageChannel, discord.TextChannel)
        self._build_channel_options(
            limited_channels=limited_channels, default=False, types=types
        )
        self.channel_select.placeholder = "Select a channel"
        self.__channel_snowflake = None
        self.channel_select.disabled = False
        await self._sync_duration_select(rebuild_options=True)
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select a channel",
        options=[discord.SelectOption(label="Select a guild first", value=str(None))],
    )
    async def channel_select(self, interaction, select):
        await interaction.response.defer()
        bot: DiscordBot = DiscordBot.get_instance()
        for option in self.channel_select.options:
            option.default = False
        channel = bot.get_channel(int(select.values[0]))
        if isinstance(channel, discord.abc.GuildChannel):
            self.channel_select.placeholder = channel.name
            self.__channel_snowflake = channel.id
            self.__guild_snowflake = channel.guild.id
        await self._sync_duration_select(rebuild_options=True)
        await interaction.edit_original_response(view=self)

    @discord.ui.select(placeholder="Select a duration", options=[])
    async def duration_select(self, interaction, select):
        await interaction.response.defer()
        duration_name = next(
            option.label
            for option in select.options
            if option.value == select.values[0]
        )
        duration_value = select.values[0]
        self.duration_select.placeholder = duration_name
        duration_builder = DurationBuilder()
        self.__duration = duration_builder.parse(duration_value).build()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        if self.__channel_snowflake is None or self.__guild_snowflake is None:
            return await interaction.response.send_message(
                content=f"Please select all fields.",
                ephemeral=True,
            )
        self.stop()
        self.__tick.update_source(interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
            author_snowflake=self.__author_snowflake,
            member_snowflake=self.__ctx.member_snowflake,
        )
        database_factory: DatabaseFactory = DatabaseFactory(self.__model)
        records = await database_factory.select(
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
            member_snowflake=self.__ctx.member_snowflake,
            singular=False,
        )
        duration_builder = DurationBuilder()
        if records:
            await interaction.response.defer()
            embed = None
            is_channel_scope = False
            target = None
            for record in records:
                if hasattr(record, "expires_in"):
                    duration = duration_builder.from_timestamp(
                        record.expires_in
                    ).build()
                    if await cap_service.exceeds_cap(
                        channel_snowflake=self.__channel_snowflake,
                        category=self.__model.identifier,
                        duration=duration,
                        guild_snowflake=self.__guild_snowflake,
                    ):
                        try:
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=self.__channel_snowflake,
                                guild_snowflake=self.__guild_snowflake,
                                member_snowflake=self.__author_snowflake,
                                requested=["command.moderation.uncapped"],
                            )
                        except CheckFailure:
                            return await interaction.response.send_message(
                                content=f"Duration {duration_builder.load(duration).as_str()} exceeds the channel cap.",
                                ephemeral=True,
                            )
                match record:
                    case Ban() as ban:
                        if ban.blacklisted:
                            return (
                                await unban_alias_service.build_blacklisted_block_embed(
                                    channel_snowflake=self.__channel_snowflake,
                                    guild_snowflake=self.__guild_snowflake,
                                    member_snowflake=self.__ctx.member_snowflake,
                                )
                            )
                        await database_factory.delete_by_cls(
                            ban,
                            channel_snowflake=self.__channel_snowflake,
                            guild_snowflake=self.__guild_snowflake,
                            member_snowflake=self.__ctx.member_snowflake,
                        )
                        is_channel_scope = await unban_alias_service.disable(
                            channel_snowflake=self.__channel_snowflake,
                            guild_snowflake=self.__guild_snowflake,
                            member_snowflake=self.__ctx.member_snowflake,
                            reason="Unban",
                        )
                        embed = unban_alias_service.build_unban_embed(
                            channel_snowflake=self.__channel_snowflake,
                            guild_snowflake=self.__guild_snowflake,
                            member_snowflake=self.__ctx.member_snowflake,
                        )
                    case Flag() as flag:
                        await database_factory.delete_by_cls(
                            flag,
                            channel_snowflake=self.__channel_snowflake,
                            guild_snowflake=self.__guild_snowflake,
                            member_snowflake=self.__ctx.member_snowflake,
                        )
                        is_channel_scope = await unflag_alias_service.disable(
                            channel_snowflake=self.__channel_snowflake,
                            guild_snowflake=self.__guild_snowflake,
                            member_snowflake=self.__ctx.member_snowflake,
                        )
                        embed = unflag_alias_service.build_unflag_embed(
                            channel_snowflake=self.__channel_snowflake,
                            guild_snowflake=self.__guild_snowflake,
                            member_snowflake=self.__ctx.member_snowflake,
                        )
                    case TextMute() as tmute:
                        await database_factory.delete_by_cls(
                            tmute,
                            channel_snowflake=self.__channel_snowflake,
                            guild_snowflake=self.__guild_snowflake,
                            member_snowflake=self.__ctx.member_snowflake,
                        )
                        is_channel_scope = await untext_mute_alias_service.disable(
                            channel_snowflake=self.__channel_snowflake,
                            guild_snowflake=self.__guild_snowflake,
                            member_snowflake=self.__ctx.member_snowflake,
                            reason="Untmute",
                        )
                        embed = untext_mute_alias_service.build_untext_mute_embed(
                            channel_snowflake=self.__channel_snowflake,
                            guild_snowflake=self.__guild_snowflake,
                            member_snowflake=self.__ctx.member_snowflake,
                        )
                    case VoiceMute() as vmute:
                        if vmute.target in ["auto", "command"]:
                            await database_factory.delete_by_cls(
                                vmute,
                                channel_snowflake=self.__channel_snowflake,
                                guild_snowflake=self.__guild_snowflake,
                                member_snowflake=self.__ctx.member_snowflake,
                                target=record.target,
                            )
                            target = record.target
                            is_channel_scope = await unvoice_mute_alias_service.disable(
                                channel_snowflake=self.__channel_snowflake,
                                guild_snowflake=self.__guild_snowflake,
                                member_snowflake=self.__ctx.member_snowflake,
                                reason="Unvmute",
                            )
                            embed = unvoice_mute_alias_service.build_unvoice_mute_embed(
                                channel_snowflake=self.__channel_snowflake,
                                guild_snowflake=self.__guild_snowflake,
                                member_snowflake=self.__ctx.member_snowflake,
                            )
                await stream_service.log(
                    author_snowflake=self.__author_snowflake,
                    channel_snowflake=self.__channel_snowflake,
                    display=True,
                    duration=DurationObject(number=0, prefix="", sign=1, unit=""),
                    guild_snowflake=self.__guild_snowflake,
                    identifier=record.identifier,
                    is_channel_scope=is_channel_scope,
                    member_snowflake=self.__ctx.member_snowflake,
                    message_snowflake=None,
                    message_channel_snowflake=None,
                    reason="No reason provided.",
                    role_snowflake=None,
                    target=target,
                )
            if embed:
                await self.__tick.end(success=embed)
        else:
            if self.__duration:
                if await cap_service.exceeds_cap(
                    channel_snowflake=self.__channel_snowflake,
                    category=self.__model.identifier,
                    duration=self.__duration,
                    guild_snowflake=self.__guild_snowflake,
                ):
                    try:
                        bot: DiscordBot = DiscordBot.get_instance()
                        permission_state: PermissionState = bot.registry.get(
                            PermissionState
                        )
                        await permission_service.has_permissions(
                            permission_state=permission_state,
                            channel_snowflake=self.__channel_snowflake,
                            guild_snowflake=self.__guild_snowflake,
                            member_snowflake=self.__author_snowflake,
                            requested=["command.moderation.uncapped"],
                        )
                    except CheckFailure:
                        return await interaction.response.send_message(
                            content=f"Duration {duration_builder.load(self.__duration).as_str()} exceeds the channel cap.",
                            ephemeral=True,
                        )
            modal = InfractionModal(
                author_snowflake=self.__author_snowflake,
                channel_snowflake=self.__channel_snowflake,
                duration=self.__duration,
                guild_snowflake=self.__guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                model=self.__model,
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
