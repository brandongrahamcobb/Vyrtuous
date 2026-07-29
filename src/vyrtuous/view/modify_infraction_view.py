"""!/bin/python3
modify_infraction_view.py The purpose of this program is to provide the view for modifying an infraction.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.permissions import PermissionScope
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.modal.duration_modal import DurationModal
from vyrtuous.modal.reason_modal import ReasonModal
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)
from vyrtuous.permissions import permission_service
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
        self.__channel_snowflake = None
        self.__ctx = ctx
        self.__guild_snowflake: int
        self.__is_channel_scope: bool = False
        self.__modal = modal
        self.__tick = tick

    # async def interaction_check(self, interaction):
    #     return interaction.user.id == self.__author_snowflake
    #
    # async def setup(self):
    #     available_channels: list[discord.abc.GuildChannel] = []
    #     available_guilds: list[discord.Guild] = []
    #     bot: DiscordBot = DiscordBot.get_instance()
    #     permission_state: PermissionState = bot.registry.get(PermissionState)
    #     all_assigned_groups = await permission_service.resolve_all_assigned_groups(
    #         permission_state=permission_state, member_snowflake=self.__author_snowflake
    #     )
    #     for group, _, _ in all_assigned_groups:
    #         if group.scope == PermissionScope.GLOBAL:
    #             for guild in bot.guilds:
    #                 available_channels.extend(guild.channels)
    #                 available_guilds.append(guild)
    #             break
    #         elif group.scope == PermissionScope.GUILD:
    #             for guild in bot.guilds:
    #                 if permission_service.has_permissions(
    #                     permission_state=permission_state,
    #                     member_snowflake=self.__author_snowflake,
    #                     requested=["other_guilds"],
    #                 ):
    #                     available_channels.extend(guild.channels)
    #                     available_guilds.append(guild)
    #                 elif guild == self.__ctx.guild_snowflake:
    #                     available_channels.extend(guild.channels)
    #                     available_guilds.append(guild)
    #             break
    #         elif group.scope == PermissionScope.CHANNEL:
    #             for guild in bot.guilds:
    #                 for channel in guild.channels:
    #                     effective_group = permission_service.resolve_effective_group(
    #                         permission_state=permission_state,
    #                         member_snowflake=self.__author_snowflake,
    #                         channel_snowflake=channel.id,
    #                         guild_snowflake=guild.id,
    #                     )
    #                     if effective_group == group:
    #                         available_guilds.append(guild)
    #                         available_channels.append(channel)
    #             break
    #         else:
    #             raise app_commands.CheckFailure(
    #                 "You do not have sufficient privileges in this channel or server to use this command."
    #             )
    #     self._build_guild_options(
    #         available_channels=available_channels,
    #         available_guilds=available_guilds,
    #     )
    #     limited_channels = self.limit_available_to_top_24_by_member_count(
    #         available=available_channels
    #     )
    #     self._build_channel_options(limited_channels=limited_channels)
    #     await self.build_duration_options(available_durations=available_durations)
    #
    # def limit_available_to_top_25_by_member_count(self, available):
    #     items = []
    #     items.extend(available)
    #     items.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
    #     top_25 = items[:25]
    #     return top_25
    #
    # def _build_guild_options(self):
    #     guild_options = [
    #         discord.SelectOption(label=g.name, value=str(g.id))
    #         for g in self.__available_guilds
    #         if not isinstance(g, str) and g is not None
    #     ]
    #     return guild_options
    #
    # def _build_channel_options(self):
    #     channel_options = [
    #         discord.SelectOption(label=c.name, value=str(c.id))
    #         for c in self.__available_channels
    #         if isinstance(c, (discord.VoiceChannel, discord.StageChannel))
    #     ]
    #     if "all" in self.__available_channels:
    #         channel_options.append(discord.SelectOption(label="All", value="all"))
    #     return channel_options
    #
    # async def _build_category_options(
    #     self, channel_snowflake: int, guild_snowflake: int
    # ):
    #     categories = []
    #     for model in INFRACTION_MODELS:
    #         if model == VoiceMute:
    #             database_factory: DatabaseFactory = DatabaseFactory(AutoMute)
    #             automute = await database_factory.select(
    #                 channel_snowflake=channel_snowflake,
    #                 guild_snowflake=guild_snowflake,
    #                 singular=True,
    #             )
    #             if automute:
    #                 target = "auto"
    #             else:
    #                 target = "command"
    #             database_factory: DatabaseFactory = DatabaseFactory(model)
    #             infraction = await database_factory.select(
    #                 channel_snowflake=channel_snowflake,
    #                 guild_snowflake=guild_snowflake,
    #                 member_snowflake=self.__ctx.member_snowflake,
    #                 target=target,
    #                 singular=True,
    #             )
    #         else:
    #             database_factory: DatabaseFactory = DatabaseFactory(model)
    #             infraction = await database_factory.select(
    #                 channel_snowflake=channel_snowflake,
    #                 guild_snowflake=guild_snowflake,
    #                 member_snowflake=self.__ctx.member_snowflake,
    #                 singular=True,
    #             )
    #         if infraction:
    #             categories.append(model.identifier)
    #     category_options = [
    #         discord.SelectOption(label=c, value=str(c)) for c in categories
    #     ]
    #     return category_options
    #
    # @discord.ui.select(
    #     placeholder="Select a guild",
    #     options=[],
    # )
    # async def guild_select(self, interaction, select):
    #     bot: DiscordBot = DiscordBot.get_instance()
    #     guild = bot.get_guild(int(select.values[0]))
    #     if guild is None:
    #         raise commands.GuildNotFound(str(select.values[0]))
    #     self.guild_select.placeholder = guild.name
    #     self.__guild_snowflake = guild.id
    #     available_channels = [
    #         c
    #         for c in self.__available_channels
    #         if isinstance(c, (discord.VoiceChannel, discord.StageChannel))
    #         and c.guild.id == self.__guild_snowflake
    #     ]
    #     self.__available_channels = self.limit_available_to_top_25_by_member_count(
    #         available=available_channels
    #     )
    #     channel_options = self._build_channel_options()
    #     self.channel_select.options = channel_options
    #     self.channel_select.disabled = False
    #     await interaction.response.defer()
    #     await interaction.edit_original_response(view=self)
    #
    # @discord.ui.select(
    #     placeholder="Select a channel",
    #     options=[discord.SelectOption(label="Select a guild", value="None")],
    #     disabled=True,
    # )
    # async def channel_select(self, interaction, select):
    #     bot: DiscordBot = DiscordBot.get_instance()
    #     permission_state: PermissionState = bot.registry.get(PermissionState)
    #     channel = interaction.guild.get_channel(int(select.values[0]))
    #     self.channel_select.placeholder = channel.name
    #     self.__channel_snowflake = channel.id
    #     bot: DiscordBot = DiscordBot.get_instance()
    #     channel = bot.get_channel(self.__channel_snowflake)
    #     if channel is None:
    #         raise commands.ChannelNotFound(str(self.__channel_snowflake))
    #     if isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
    #         for member in channel.members:
    #             if self.__ctx.member_snowflake == member.id:
    #                 self.__is_channel_scope = True
    #     if isinstance(channel, discord.abc.GuildChannel):
    #         await permission_service.has_equal_or_lower_role(
    #             permission_state=permission_state,
    #             channel_snowflake=self.__channel_snowflake,
    #             guild_snowflake=self.__guild_snowflake,
    #             author_snowflake=self.__author_snowflake,
    #             member_snowflake=self.__ctx.member_snowflake,
    #         )
    #     category_options = await self._build_category_options(
    #         channel_snowflake=self.__channel_snowflake,
    #         guild_snowflake=self.__guild_snowflake,
    #     )
    #     self.category_select.options = category_options
    #     self.category_select.disabled = False
    #     await interaction.response.defer()
    #     await interaction.edit_original_response(view=self)
    #
    # @discord.ui.select(
    #     placeholder="Select a category",
    #     options=[discord.SelectOption(label="Select a channel", value="None")],
    #     disabled=True,
    # )
    # async def category_select(self, interaction, select):
    #     self.category_select.placeholder = select.values[0]
    #     self.__category = select.values[0]
    #     await interaction.response.defer()
    #     await interaction.edit_original_response(view=self)
    #
    # @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    # async def submit(self, interaction, button):
    #     if self.__channel_snowflake is None:
    #         return await interaction.response.send_message(
    #             content="Please select all fields.", ephemeral=True
    #         )
    #     match self.__modal:
    #         case "duration":
    #             modal = DurationModal(
    #                 author_snowflake=self.__author_snowflake,
    #                 category=self.__category,
    #                 channel_snowflake=self.__channel_snowflake,
    #                 guild_snowflake=self.__ctx.guild_snowflake,
    #                 member_snowflake=self.__ctx.member_snowflake,
    #                 tick=self.__tick,
    #             )
    #             await modal.setup()
    #             await interaction.response.send_modal(modal)
    #         case "reason":
    #             modal = ReasonModal(
    #                 author_snowflake=self.__author_snowflake,
    #                 category=self.__category,
    #                 channel_snowflake=self.__channel_snowflake,
    #                 duration=None,
    #                 guild_snowflake=self.__ctx.guild_snowflake,
    #                 member_snowflake=self.__ctx.member_snowflake,
    #                 is_channel_scope=self.__is_channel_scope,
    #                 is_modification=True,
    #                 tick=self.__tick,
    #             )
    #             await modal.setup(is_new=False)
    #             await interaction.response.send_modal(modal)
    #
    # @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    # async def cancel(self, interaction, button):
    #     self.stop()
    #     await interaction.response.edit_message(content="Cancelled action.", view=None)
    #
    # async def on_error(
    #     self,
    #     interaction: discord.Interaction,
    #     error: Exception,
    #     item: discord.ui.Item,
    # ) -> None:
    #     if isinstance(error, commands.BadArgument):
    #         await self.__tick.end(error=str(error), ephemeral=True)
    #     elif isinstance(error, commands.CheckFailure):
    #         await self.__tick.end(error=str(error), ephemeral=True)
    #     else:
    #         raise error
