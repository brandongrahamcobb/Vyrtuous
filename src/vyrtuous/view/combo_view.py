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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.permissions import PermissionScope
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.cap import Cap
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.flag import Flag
from vyrtuous.db.text_mute import TextMute
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.modal.combo_modal import ComboModal
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import CheckFailure, GuildNotFound
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import cap_service

DURATIONS = {"1hour", "8hour", "1day", "1week", "30day"}

MODELS = [Flag, TextMute, VoiceMute]

class ToggleButton(discord.ui.Button):
    def __init__(self, label: str, attribute: str, view: 'ComboView'):
        super().__init__(label=label, style=discord.ButtonStyle.success)
        self.__attribute = attribute
        self.__view = view

    async def callback(self, interaction: discord.Interaction):
        enabled = await getattr(self.__view, self.__attribute)()
        self.style = discord.ButtonStyle.success if enabled else discord.ButtonStyle.danger
        await interaction.response.edit_message(view=self.view)

class ComboView(discord.ui.View):
    def __init__(
        self,
        author_snowflake: int,
        ctx: SnowflakeContext,
        interaction: discord.Interaction,
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
        self.__records = []
        self.__tick = tick
        self.__flag_enabled = True
        self.__text_mute_enabled = True
        self.__voice_mute_enabled = True
        self.add_item(ToggleButton('Flag', 'toggle_flag', self))
        self.add_item(ToggleButton('Text-Mute', 'toggle_text_mute', self))
        self.add_item(ToggleButton('Voice-Mute', 'toggle_voice_mute', self))

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
        types = (
            discord.VoiceChannel,
            discord.StageChannel,
            discord.TextChannel,
        )
        self._build_channel_options(
            limited_channels=limited_channels, default=True, types=types
        )
        await self.build_duration_options()
        if (
            not self.guild_select.options
            or not self.channel_select.options
            or not self.duration_select.options
        ):
            await self.__tick.end(
                warning="You have insufficient privileges to do that."
            )

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
                caps = []
                for model in MODELS:
                    cap = await database_factory.select(
                        category=model.identifier,
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        singular=True,
                    )
                    if cap:
                        caps.append(cap.duration_seconds)
                if caps:
                    cap_seconds = min(caps)
                    available_durations = self.limit_available_to_cap(
                        duration_seconds=cap_seconds
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
        for cmp in DURATIONS:
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

    async def toggle_flag(self):
        self.__flag_enabled = not self.__flag_enabled
        return self.__flag_enabled

    async def toggle_text_mute(self):
        self.__text_mute_enabled = not self.__text_mute_enabled
        return self.__text_mute_enabled

    async def toggle_voice_mute(self):
        self.__voice_mute_enabled = not self.__voice_mute_enabled
        return self.__voice_mute_enabled

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
        types = (discord.VoiceChannel, discord.StageChannel, discord.TextChannel)
        self._build_channel_options(
            limited_channels=limited_channels, default=False, types=types
        )
        self.channel_select.placeholder = "Select a channel"
        self.__channel_snowflake = None
        self.channel_select.disabled = False
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
        if (
            self.__channel_snowflake is None
            or self.__duration is None
            or self.__guild_snowflake is None
        ):
            return await interaction.response.send_message(
                content=f"Please select all fields.",
                ephemeral=True,
            )
        self.__tick.update_source(interaction=interaction)
        if (
            not self.__flag_enabled and not self.__text_mute_enabled and not self.__voice_mute_enabled
        ):
            return await interaction.response.send_message(
                content=f"None of the available actions were selected.",
                ephemeral=True,
            )
        self.stop()
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
            author_snowflake=self.__author_snowflake,
            member_snowflake=self.__ctx.member_snowflake,
        )
        duration_builder = DurationBuilder()
        for model in MODELS:
            database_factory: DatabaseFactory = DatabaseFactory(model)
            if model == VoiceMute:
                if self.__voice_mute_enabled:
                    command_record = await database_factory.select(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__ctx.member_snowflake,
                        target="command",
                        singular=True,
                    )
                    if command_record:
                        self.__records.append(command_record)
            if model == Flag:
                if self.__flag_enabled:
                    to_add = await database_factory.select(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__ctx.member_snowflake,
                        singular=True,
                    )
                    if to_add:
                        self.__records.append(to_add)
            if model == TextMute:
                if self.__text_mute_enabled:
                    to_add = await database_factory.select(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__ctx.member_snowflake,
                        singular=True,
                    )
                    if to_add:
                        self.__records.append(to_add)
        if self.__records:
            for record in self.__records:
                if hasattr(record, "expires_in"):
                    duration = duration_builder.from_timestamp(
                        record.expires_in
                    ).build()
                    if await cap_service.exceeds_cap(
                        channel_snowflake=self.__channel_snowflake,
                        category=record.identifier,
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
        modal = ComboModal(
            author_snowflake=self.__author_snowflake,
            channel_snowflake=self.__channel_snowflake,
            duration=self.__duration,
            guild_snowflake=self.__guild_snowflake,
            flag_enabled=self.__flag_enabled,
            member_snowflake=self.__ctx.member_snowflake,
            records=self.__records,
            text_mute_enabled=self.__text_mute_enabled,
            tick=self.__tick,
            voice_mute_enabled=self.__voice_mute_enabled
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
