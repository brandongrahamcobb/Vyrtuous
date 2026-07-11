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

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.administrator import Administrator
from vyrtuous.db.coordinator import Coordinator
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.developer import Developer
from vyrtuous.db.moderator import Moderator
from vyrtuous.modal.reason_modal import ReasonModal
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import cap_service
from vyrtuous.view.view_context import ViewContext


class InfractionView(discord.ui.View):
    def __init__(
        self,
        author_snowflake: int,
        ctx: ViewContext,
        tick: Tick,
    ):
        super().__init__(timeout=120)
        self.__author_snowflake = author_snowflake
        self.__available_channels: list[discord.abc.GuildChannel] = []
        self.__available_guilds: list[discord.Guild] = []
        self.__channel_snowflake: int
        self.__ctx = ctx
        self.__duration_value: str
        self.__is_channel_scope: bool = False
        self.__tick = tick

    async def interaction_check(self, interaction):
        return interaction.user.id == self.__author_snowflake

    async def setup(self):
        bot: DiscordBot = DiscordBot.get_instance()
        database_factory: DatabaseFactory = DatabaseFactory(Developer)
        developer: Developer = await database_factory.select(
            member_snowflake=self.__ctx.member_snowflake, singular=True
        )
        if developer is not None:
            self.__available_guilds.extend(bot.guilds)
            for guild in bot.guilds:
                for channel in guild.channels:
                    self.__available_channels.append(channel)
        database_factory: DatabaseFactory = DatabaseFactory(Administrator)
        administrators: list[Administrator] = await database_factory.select(
            member_snowflake=self.__ctx.member_snowflake, singular=False
        )
        if administrators is not None:
            for administrator in administrators:
                guild = bot.get_guild(administrator.guild_snowflake)
                if guild is None:
                    continue
                self.__available_guilds.append(guild)
                self.__available_channels.extend(guild.channels)
        database_factory: DatabaseFactory = DatabaseFactory(Coordinator)
        coordinators: list[Coordinator] = await database_factory.select(
            member_snowflake=self.__ctx.member_snowflake, singular=False
        )
        if coordinators is not None:
            for coordinator in coordinators:
                guild = bot.get_guild(coordinator.guild_snowflake)
                if guild is None:
                    continue
                channel = bot.get_channel(coordinator.channel_snowflake)
                if channel is None or not isinstance(channel, discord.abc.GuildChannel):
                    continue
                self.__available_guilds.append(guild)
                self.__available_channels.append(channel)

        database_factory: DatabaseFactory = DatabaseFactory(Moderator)
        moderators: list[Moderator] = await database_factory.select(
            member_snowflake=self.__ctx.member_snowflake, singular=False
        )
        if moderators is not None:
            for moderator in moderators:
                guild = bot.get_guild(moderator.guild_snowflake)
                if guild is None:
                    continue
                channel = bot.get_channel(moderator.channel_snowflake)
                if channel is None or not isinstance(channel, discord.abc.GuildChannel):
                    continue
                self.__available_guilds.append(guild)
                self.__available_channels.append(channel)
        self.__available_channels = self.limit_available_to_top_25_by_member_count(
            available=self.__available_channels
        )
        self.__available_guilds = self.limit_available_to_top_25_by_member_count(
            available=self.__available_guilds
        )
        channel_options = self._build_channel_options()
        duration_options = self._build_duration_options()
        self.channel_select.options = channel_options
        self.duration_select.options = duration_options

    def limit_available_to_top_25_by_member_count(self, available):
        all_key = "all"
        items = []
        if all_key in available:
            items.append(all_key)
        else:
            items.extend(available)
        items.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
        top_25 = items[:25]
        return top_25

    def _build_channel_options(self):
        channel_options = [
            discord.SelectOption(label=c.name, value=str(c.id))
            for c in self.__available_channels
            if c != "all" and c is not None
        ]
        if "all" in self.__available_channels:
            channel_options.append(discord.SelectOption(label="All", value="all"))
        bot: DiscordBot = DiscordBot.get_instance()
        bot.logger.info(len(channel_options))
        return channel_options

    def _build_duration_options(self):
        durations = ["0", "1h", "8h", "1d", "1w"]
        return [
            discord.SelectOption(
                label=duration if duration != "0" else "Permanent", value=duration
            )
            for duration in durations
        ]

    @discord.ui.select(
        placeholder="Select channel",
        options=[],
    )
    async def channel_select(self, interaction, select):
        channel = interaction.guild.get_channel(int(select.values[0]))
        self.channel_select.placeholder = channel.name
        self.__channel_snowflake = channel.id
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(self.__ctx.guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(self.__ctx.guild_snowflake))
        channel = guild.get_channel(self.__channel_snowflake)
        if channel is None:
            raise commands.ChannelNotFound(str(self.__channel_snowflake))
        if isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
            for member in channel.members:
                if self.__ctx.member_snowflake == member.id:
                    self.__is_channel_scope = True
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select duration",
        options=[],
    )
    async def duration_select(self, interaction, select):
        duration_name = select.values[0]
        self.duration_select.placeholder = duration_name
        self.__duration_value = duration_name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        if not self.has_the_user_selected_all_fields():
            return await interaction.response.send_message(
                content="Please select all fields.", ephemeral=True
            )
        if await cap_service.exceeds_cap(
            category=self.__ctx.category,
            channel_snowflake=self.__channel_snowflake,
            duration_value=self.__duration_value,
            guild_snowflake=self.__ctx.guild_snowflake,
        ):
            return await interaction.response.send_message(
                content="Duration exceeds the channel cap.", ephemeral=True
            )
        modal = ReasonModal(
            category=self.__ctx.category,
            channel_snowflake=self.__channel_snowflake,
            duration_value=self.__duration_value,
            guild_snowflake=self.__ctx.guild_snowflake,
            member_snowflake=self.__ctx.member_snowflake,
            is_channel_scope=self.__is_channel_scope,
            tick=self.__tick,
        )
        await modal.setup()
        await interaction.response.send_modal(modal)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        await interaction.message.delete()
        self.stop()
        return await self.__tick.end(warning="Cancelled action.")

    def has_the_user_selected_all_fields(self):
        if not self.__duration_value or not self.__channel_snowflake:
            return False
        return True
