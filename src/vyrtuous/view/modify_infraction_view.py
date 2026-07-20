"""!/bin/python3
new_infraction_view.py The purpose of this program is to provide the view for creating an infraction.

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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.administrator import Administrator
from vyrtuous.db.coordinator import Coordinator
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.developer import Developer
from vyrtuous.db.moderator import Moderator
from vyrtuous.modal.duration_modal import DurationModal
from vyrtuous.modal.reason_modal import ReasonModal
from vyrtuous.models.duration import DurationObject
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.users import guild_owner_service, sysadmin_service
from vyrtuous.view.view_context import ViewContext

INFRACTION_MODELS = [
    ban_service.MODEL,
    flag_service.MODEL,
    text_mute_service.MODEL,
    voice_mute_service.MODEL,
]


class ModifyInfractionView(discord.ui.View):
    def __init__(
        self,
        author_snowflake: int,
        ctx: ViewContext,
        modal: str,
        tick: Tick,
    ):
        super().__init__(timeout=120)
        self.__author_snowflake = author_snowflake
        self.__all_available_channels: list[discord.abc.GuildChannel | str] = []
        self.__available_channels: list[discord.abc.GuildChannel | str] = []
        self.__available_guilds: list[discord.Guild | str] = []
        self.__category: str
        self.__channel_snowflake: int
        self.__ctx = ctx
        self.__duration: DurationObject
        self.__guild_snowflake: int
        self.__is_channel_scope: bool = False
        self.__modal = modal
        self.__tick = tick

    async def interaction_check(self, interaction):
        return interaction.user.id == self.__author_snowflake

    async def setup(self):
        guild_select_disabled = True
        bot: DiscordBot = DiscordBot.get_instance()
        database_factory: DatabaseFactory = DatabaseFactory(Developer)
        developer: Developer = await database_factory.select(
            member_snowflake=self.__author_snowflake, singular=True
        )
        while True:
            if developer is not None or sysadmin_service.is_sysadmin(
                self.__author_snowflake
            ):
                guild_select_disabled = False
                self.__available_guilds.extend(bot.guilds)
                self.__available_guilds.append("all")
                for guild in bot.guilds:
                    self.__all_available_channels.extend(guild.channels)
                    self.__all_available_channels.append("all")
                break
            if await guild_owner_service.is_guild_owner(
                guild_snowflake=self.__ctx.guild_snowflake,
                member_snowflake=self.__author_snowflake,
            ):
                guild = bot.get_guild(self.__ctx.guild_snowflake)
                if guild is None:
                    raise commands.GuildNotFound(str(self.__ctx.guild_snowflake))
                self.__available_guilds.append(guild)
                self.__all_available_channels.extend(guild.channels)
                break
            database_factory: DatabaseFactory = DatabaseFactory(Administrator)
            administrators: list[Administrator] = await database_factory.select(
                member_snowflake=self.__author_snowflake,
                guild_snowflake=self.__ctx.guild_snowflake,
                singular=False,
            )
            if administrators:
                for administrator in administrators:
                    guild = bot.get_guild(administrator.guild_snowflake)
                    if guild is None:
                        continue
                    self.__available_guilds.append(guild)
                    self.__all_available_channels.extend(guild.channels)
                break
            database_factory: DatabaseFactory = DatabaseFactory(Coordinator)
            coordinators: list[Coordinator] = await database_factory.select(
                member_snowflake=self.__author_snowflake, singular=False
            )
            if coordinators:
                for coordinator in coordinators:
                    guild = bot.get_guild(coordinator.guild_snowflake)
                    if guild is None:
                        continue
                    channel = bot.get_channel(coordinator.channel_snowflake)
                    if channel is None or not isinstance(
                        channel, discord.abc.GuildChannel
                    ):
                        continue
                    self.__available_guilds.append(guild)
                    self.__all_available_channels.append(channel)
                break
            database_factory: DatabaseFactory = DatabaseFactory(Moderator)
            moderators: list[Moderator] = await database_factory.select(
                member_snowflake=self.__author_snowflake, singular=False
            )
            if moderators:
                for moderator in moderators:
                    guild = bot.get_guild(moderator.guild_snowflake)
                    if guild is None:
                        continue
                    channel = bot.get_channel(moderator.channel_snowflake)
                    if channel is None or not isinstance(
                        channel, discord.abc.GuildChannel
                    ):
                        continue
                    self.__available_guilds.append(guild)
                    self.__all_available_channels.append(channel)
                break
            break
        if guild_select_disabled:
            self.remove_item(self.guild_select)
            self.__available_channels = self.limit_available_to_top_25_by_member_count(
                available=self.__all_available_channels
            )
            channel_options = self._build_channel_options()
            self.channel_select.options = channel_options
            self.channel_select.disabled = False
        else:
            self.__available_guilds = self.limit_available_to_top_25_by_member_count(
                available=self.__available_guilds
            )
            self.channel_select.disabled = True
            guild_options = self._build_guild_options()
            self.guild_select.options = guild_options

    def limit_available_to_top_25_by_member_count(self, available):
        items = []
        items.extend(available)
        items.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
        top_25 = items[:25]
        return top_25

    def _build_guild_options(self):
        guild_options = [
            discord.SelectOption(label=g.name, value=str(g.id))
            for g in self.__available_guilds
            if not isinstance(g, str) and g is not None
        ]
        return guild_options

    def _build_channel_options(self):
        channel_options = [
            discord.SelectOption(label=c.name, value=str(c.id))
            for c in self.__available_channels
            if isinstance(c, (discord.VoiceChannel, discord.StageChannel))
        ]
        if "all" in self.__available_channels:
            channel_options.append(discord.SelectOption(label="All", value="all"))
        return channel_options

    async def _build_category_options(
        self, channel_snowflake: int, guild_snowflake: int
    ):
        categories = []
        for model in INFRACTION_MODELS:
            database_factory: DatabaseFactory = DatabaseFactory(model)
            infraction = await database_factory.select(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                singular=True,
            )
            if infraction:
                categories.append(model.identifier)
        category_options = [
            discord.SelectOption(label=c, value=str(c)) for c in categories
        ]
        return category_options

    @discord.ui.select(
        placeholder="Select a guild",
        options=[],
    )
    async def guild_select(self, interaction, select):
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(int(select.values[0]))
        if guild is None:
            raise commands.GuildNotFound(str(select.values[0]))
        self.guild_select.placeholder = guild.name
        self.__guild_snowflake = guild.id
        available_channels = [
            c
            for c in self.__available_channels
            if isinstance(c, (discord.VoiceChannel, discord.StageChannel))
            and c.guild.id == self.__guild_snowflake
        ]
        self.__available_channels = self.limit_available_to_top_25_by_member_count(
            available=available_channels
        )
        channel_options = self._build_channel_options()
        self.channel_select.options = channel_options
        self.channel_select.disabled = False
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select a channel",
        options=[discord.SelectOption(label="Select a guild", value="None")],
        disabled=True,
    )
    async def channel_select(self, interaction, select):
        channel = interaction.guild.get_channel(int(select.values[0]))
        self.channel_select.placeholder = channel.name
        self.__channel_snowflake = channel.id
        bot: DiscordBot = DiscordBot.get_instance()
        channel = bot.get_channel(self.__channel_snowflake)
        if channel is None:
            raise commands.ChannelNotFound(str(self.__channel_snowflake))
        if isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
            for member in channel.members:
                if self.__ctx.member_snowflake == member.id:
                    self.__is_channel_scope = True
        category_options = await self._build_category_options(
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
        )
        self.category_select.options = category_options
        self.category_select.disabled = False
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select a category",
        options=[discord.SelectOption(label="Select a channel", value="None")],
        disabled=True,
    )
    async def category_select(self, interaction, select):
        self.category_select.placeholder = select.values[0]
        self.__category = select.values[0]
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        if not self.has_the_user_selected_all_fields():
            return await interaction.response.send_message(
                content="Please select all fields.", ephemeral=True
            )
        match self.__modal:
            case "duration":
                modal = DurationModal(
                    author_snowflake=self.__author_snowflake,
                    category=self.__category,
                    channel_snowflake=self.__channel_snowflake,
                    guild_snowflake=self.__ctx.guild_snowflake,
                    member_snowflake=self.__ctx.member_snowflake,
                    tick=self.__tick,
                )
                await modal.setup()
                await interaction.response.send_modal(modal)
            case "reason":
                modal = ReasonModal(
                    author_snowflake=self.__author_snowflake,
                    category=self.__category,
                    channel_snowflake=self.__channel_snowflake,
                    duration=None,
                    guild_snowflake=self.__ctx.guild_snowflake,
                    member_snowflake=self.__ctx.member_snowflake,
                    is_channel_scope=self.__is_channel_scope,
                    is_modification=True,
                    tick=self.__tick,
                )
                await modal.setup(is_new=False)
                await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        self.stop()
        await interaction.response.edit_message(content="Cancelled action.", view=None)

    def has_the_user_selected_all_fields(self):
        if not self.__channel_snowflake:
            return False
        return True
