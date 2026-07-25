"""!/bin/python3
infraction_view.py The purpose of this program is to provide the view for creating an infraction.

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

import re

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.administrator import Administrator
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.cap import Cap
from vyrtuous.db.coordinator import Coordinator
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.flag import Flag
from vyrtuous.db.moderator import Moderator
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.modal.reason_modal import ReasonModal
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    cap_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.users import moderator_service
from vyrtuous.view.view_context import ViewContext

INFRACTION_MODELS = [
    ban_service.MODEL,
    flag_service.MODEL,
    text_mute_service.MODEL,
    voice_mute_service.MODEL,
]

UNDO_CATEGORIES = {
    "ban": "ban",
    "flag": "flag",
    "tmute": "text-mute",
    "vmute": "voice-mute",
}


class InfractionView(discord.ui.View):
    def __init__(
        self,
        author_snowflake: int,
        ctx: ViewContext,
        tick: Tick,
    ):
        super().__init__(timeout=120)
        self.__author_snowflake = author_snowflake
        self.__channel_snowflake = ctx.channel_snowflake
        self.__ctx = ctx
        self.__duration = None
        self.__guild_snowflake = ctx.guild_snowflake
        self.__is_channel_scope: bool = False
        self.__tick = tick
        self.__record = None
        self.__model = None

    async def interaction_check(self, interaction):
        return interaction.user.id == self.__author_snowflake

    async def setup(self):
        available_channels: list[discord.abc.GuildChannel] = []
        available_guilds: list[discord.Guild] = []
        if self.__ctx.category == "flag":
            self.duration_select.disabled = True
        bot: DiscordBot = DiscordBot.get_instance()
        highest_role = await moderator_service.check_minimum_role(
            member_snowflake=self.__author_snowflake, lowest_role="Developer"
        )
        if highest_role in ["Sysadmin", "Developer"]:
            for guild in bot.guilds:
                available_channels.extend(guild.channels)
                available_guilds.append(guild)
            self._build_guild_options(
                available_channels=available_channels,
                available_guilds=available_guilds,
            )
            limited_channels = self.limit_available_to_top_24_by_member_count(
                available=available_channels
            )
            self._build_channel_options(limited_channels=limited_channels)
            available_durations = [
                "1hour",
                "8hours",
                "1day",
                "1week",
                "0",
            ]
            await self.build_duration_options(available_durations=available_durations)
        else:
            is_guild_owner = False
            available_channels: list[discord.abc.GuildChannel] = []
            for guild in bot.guilds:
                if guild.owner_id == self.__author_snowflake:
                    is_guild_owner = True
                    available_channels.extend(guild.channels)
                    available_guilds.append(guild)
            if is_guild_owner:
                self._build_guild_options(
                    available_channels=available_channels,
                    available_guilds=available_guilds,
                )
                limited_channels = self.limit_available_to_top_24_by_member_count(
                    available=available_channels
                )
                self._build_channel_options(limited_channels=limited_channels)
                available_durations = [
                    "1hour",
                    "8hour",
                    "1day",
                    "1week",
                    "0",
                ]
                await self.build_duration_options(
                    available_durations=available_durations
                )
            else:
                database_factory: DatabaseFactory = DatabaseFactory(Coordinator)
                administrators: list[Administrator] = await database_factory.select(
                    member_snowflake=self.__author_snowflake,
                    singular=False,
                )
                if administrators:
                    for administrator in administrators:
                        guild = bot.get_guild(administrator.guild_snowflake)
                        if guild is None:
                            continue
                        available_guilds.append(guild)
                        available_channels.extend(guild.channels)
                    self._build_guild_options(
                        available_channels=available_channels,
                        available_guilds=available_guilds,
                    )
                    limited_channels = self.limit_available_to_top_24_by_member_count(
                        available=available_channels
                    )
                    self._build_channel_options(limited_channels=limited_channels)
                    available_durations = ["1hour", "8hours", "1day", "1week", "0"]
                    await self.build_duration_options(
                        available_durations=available_durations
                    )
                else:
                    database_factory: DatabaseFactory = DatabaseFactory(Coordinator)
                    coordinators: list[Coordinator] = await database_factory.select(
                        member_snowflake=self.__author_snowflake,
                        singular=False,
                    )
                    if coordinators:
                        for coordinator in coordinators:
                            guild = bot.get_guild(coordinator.guild_snowflake)
                            if guild is None:
                                continue
                            if guild not in available_guilds:
                                available_guilds.append(guild)
                            channel = guild.get_channel(coordinator.channel_snowflake)
                            if channel is None or not isinstance(
                                channel, discord.abc.GuildChannel
                            ):
                                continue
                            available_channels.append(channel)
                        self._build_guild_options(
                            available_guilds=available_guilds,
                            available_channels=available_channels,
                        )
                        available_durations = [
                            "1hour",
                            "8hours",
                            "1day",
                            "1week",
                            "0",
                        ]
                        await self.build_duration_options(
                            available_durations=available_durations
                        )
                    else:
                        database_factory: DatabaseFactory = DatabaseFactory(Moderator)
                        moderators: list[Moderator] = await database_factory.select(
                            member_snowflake=self.__author_snowflake,
                            singular=False,
                        )
                        if moderators:
                            for moderator in moderators:
                                guild = bot.get_guild(self.__ctx.guild_snowflake)
                                if guild is None:
                                    continue
                                if guild not in available_guilds:
                                    available_guilds.append(guild)
                                channel = guild.get_channel(moderator.channel_snowflake)
                                if channel is None or not isinstance(
                                    channel, discord.abc.GuildChannel
                                ):
                                    continue
                                available_channels.append(channel)
                            self._build_guild_options(
                                available_guilds=available_guilds,
                                available_channels=available_channels,
                            )
                            limited_channels = (
                                self.limit_available_to_top_24_by_member_count(
                                    available=available_channels
                                )
                            )
                            self._build_channel_options(
                                limited_channels=limited_channels
                            )
                            await self.build_duration_options()
                        else:
                            raise app_commands.CheckFailure(
                                "You do not have sufficient privileges in this channel or server to use this command."
                            )

    async def build_duration_options(
        self, available_durations: list[str] | None = None
    ):
        def _format_duration_label(duration: str):
            if duration == "0":
                return "Permanent"
            match = re.fullmatch(r"(\d+)([a-zA-Z]+)", duration)
            if match is None:
                return duration
            amount, unit = match.groups()
            return f"{amount} {unit}"

        if available_durations is None:
            duration_builder: DurationBuilder = DurationBuilder()
            database_factory: DatabaseFactory = DatabaseFactory(Cap)
            cap: Cap = await database_factory.select(
                category=self.__ctx.category,
                channel_snowflake=self.__channel_snowflake,
                guild_snowflake=self.__ctx.guild_snowflake,
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
            duration_options = [
                discord.SelectOption(
                    label=_format_duration_label(d),
                    value=d,
                )
                for d in available_durations
            ]
            self.duration_select.options = duration_options

    def limit_available_to_cap(self, duration_seconds: int):
        duration_builder: DurationBuilder = DurationBuilder()
        duration = duration_builder.from_seconds(duration_seconds).to_timedelta()
        options = []
        durations = ["1hour", "8hours", "1day", "1week"]
        for cmp in durations:
            compare_duration = duration_builder.parse(cmp).to_timedelta()
            if duration >= compare_duration:
                options.append(cmp)
        return options

    def limit_available_to_top_24_by_member_count(self, available):
        items = []
        items.extend(available)
        items.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
        top_24 = items[:24]
        return top_24

    def _build_channel_options(
        self, limited_channels: list[discord.abc.GuildChannel], all: bool = False
    ):
        channel_options = []
        bot: DiscordBot = DiscordBot.get_instance()
        channel = bot.get_channel(self.__ctx.channel_snowflake)
        if channel:
            if isinstance(
                channel,
                (discord.TextChannel, discord.VoiceChannel, discord.StageChannel),
            ):
                channel_options.append(
                    discord.SelectOption(
                        label=channel.name,
                        value=str(self.__ctx.channel_snowflake),
                        default=True,
                    )
                )
        channel_options.extend(
            [
                discord.SelectOption(label=c.name, value=str(c.id))
                for c in limited_channels
                if isinstance(
                    c, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
                )
                and c.id != self.__ctx.channel_snowflake
            ]
        )
        if all:
            channel_options.append(discord.SelectOption(label="All", value="all"))
        self.channel_select.options = channel_options

    def _build_guild_options(
        self,
        available_channels: list[discord.abc.GuildChannel],
        available_guilds: list[discord.Guild],
    ):
        guild_options = []
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(self.__ctx.guild_snowflake)
        if guild:
            guild_options.append(
                discord.SelectOption(
                    label=guild.name,
                    value=str(self.__ctx.guild_snowflake),
                    default=True,
                )
            )
        if len(available_guilds) == 1:
            self.remove_item(self.guild_select)
            top_24_channels = self.limit_available_to_top_24_by_member_count(
                available=available_channels
            )
            self._build_channel_options(limited_channels=top_24_channels, all=False)
        else:
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
            guild_options.append(discord.SelectOption(label="All", value="all"))
            self.guild_select.options = guild_options

    @discord.ui.select(
        placeholder="Select a guild",
        options=[],
    )
    async def guild_select(self, interaction, select):
        await interaction.response.defer()
        bot: DiscordBot = DiscordBot.get_instance()
        if select.values[0] == "all":
            self.guild_select.placeholder = "All"
            all_channels = []
            for guild in bot.guilds:
                all_channels.extend(guild.channels)
            limited_channels = self.limit_available_to_top_24_by_member_count(
                available=all_channels,
            )
            self._build_channel_options(limited_channels=limited_channels, all=True)
        else:
            guild = bot.get_guild(int(select.values[0]))
            if guild is None:
                raise commands.GuildNotFound(str(select.values[0]))
            limited_channels = self.limit_available_to_top_24_by_member_count(
                available=guild.channels,
            )
            self._build_channel_options(limited_channels=limited_channels, all=True)
        self.channel_select.disabled = False
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select a channel",
        options=[discord.SelectOption(label="Select a guild first", value=str(None))],
    )
    async def channel_select(self, interaction, select):
        await interaction.response.defer()
        bot: DiscordBot = DiscordBot.get_instance()
        duration_builder: DurationBuilder = DurationBuilder()
        channel = bot.get_channel(int(select.values[0]))
        if isinstance(channel, discord.abc.GuildChannel):
            self.channel_select.placeholder = channel.name
            self.__channel_snowflake = channel.id
            self.__guild_snowflake = channel.guild.id
            await moderator_service.has_equal_or_lower_role(
                channel_snowflake=self.__channel_snowflake,
                guild_snowflake=self.__guild_snowflake,
                member_snowflake=self.__author_snowflake,
                target_member_snowflake=self.__ctx.member_snowflake,
            )
        if isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
            for member in channel.members:
                if self.__ctx.member_snowflake == member.id:
                    self.__is_channel_scope = True
        self.__model = next(
            (
                model
                for model in INFRACTION_MODELS
                if model.identifier == self.__ctx.category
            ),
            None,
        )
        database_factory: DatabaseFactory = DatabaseFactory(AutoMute)
        automute = await database_factory.select(
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__ctx.guild_snowflake,
            singular=True,
        )
        if self.__model == VoiceMute:
            if automute:
                target = "auto"
            else:
                target = "command"
            database_factory: DatabaseFactory = DatabaseFactory(self.__model)
            record = await database_factory.select(
                channel_snowflake=self.__channel_snowflake,
                guild_snowflake=self.__ctx.guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                target=target,
                singular=True,
            )
        else:
            database_factory: DatabaseFactory = DatabaseFactory(self.__model)
            record = await database_factory.select(
                channel_snowflake=self.__channel_snowflake,
                guild_snowflake=self.__ctx.guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                singular=True,
            )
        if not record and not automute and self.__model != Flag:
            highest_role = await moderator_service.resolve_highest_role(
                channel_snowflake=self.__channel_snowflake,
                guild_snowflake=self.__guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
            )
            if highest_role == "Moderator":
                await self.build_duration_options()
            else:
                available_durations = [
                    "1hour",
                    "8hours",
                    "1day",
                    "1week",
                    "0",
                ]
                await self.build_duration_options(
                    available_durations=available_durations
                )
        if record:
            self.__record = record
            self.remove_item(self.duration_select)
            self.__duration = duration_builder.parse(0).build()
            if category := UNDO_CATEGORIES.get(self.__ctx.category, None):
                await interaction.edit_original_response(
                    content=f"Are you sure you want to undo their {category}?",
                    view=self,
                )
                return None
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
        duration_builder = DurationBuilder()
        if self.__duration is None:
            if self.__model != Flag:
                return await interaction.response.send_message(
                    content="Please select all fields.", ephemeral=True
                )
        else:
            if self.__record is None:
                if await cap_service.exceeds_cap(
                    channel_snowflake=self.__channel_snowflake,
                    category=self.__ctx.category,
                    duration=self.__duration,
                    guild_snowflake=self.__guild_snowflake,
                ):
                    role = await moderator_service.resolve_highest_role(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__author_snowflake,
                    )
                    if role in ["Moderator", "Everyone"]:
                        return await interaction.response.send_message(
                            content=f"Duration {duration_builder.load(self.__duration).as_str()} exceeds the channel cap.",
                            ephemeral=True,
                        )
        modal = ReasonModal(
            author_snowflake=self.__author_snowflake,
            category=self.__ctx.category,
            channel_snowflake=self.__channel_snowflake,
            duration=self.__duration,
            guild_snowflake=self.__guild_snowflake,
            member_snowflake=self.__ctx.member_snowflake,
            is_channel_scope=self.__is_channel_scope,
            tick=self.__tick,
        )
        await modal.setup(is_new=True)
        await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        self.stop()
        await interaction.response.edit_message(content="Cancelled action.", view=None)

    async def on_error(
        self,
        interaction: discord.Interaction,
        error: Exception,
        item: discord.ui.Item,
    ) -> None:
        await self.__tick.end(error=str(error), ephemeral=True)
