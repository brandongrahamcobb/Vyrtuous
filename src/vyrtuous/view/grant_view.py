"""!/bin/python3
new_infraction_view.py The purpose of this program is to provide the view for creating an infraction.

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

from dataclasses import dataclass, field

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.permissions import PermissionGroup, PermissionScope
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.permission_entry import PermissionEntry
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import CheckFailure, MemberNotFound
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick

MODEL = PermissionEntry


@dataclass
class GroupScope:
    group: PermissionGroup
    guilds: dict[int, discord.Guild] = field(default_factory=dict)
    channels: dict[int, discord.abc.GuildChannel] = field(default_factory=dict)


SCOPE_REQUIREMENTS: dict[PermissionScope, frozenset[str]] = {
    PermissionScope.CHANNEL: frozenset({"guild", "channel"}),
    PermissionScope.GUILD: frozenset({"guild"}),
    PermissionScope.GLOBAL: frozenset(),
}


class GrantView(discord.ui.View):
    def __init__(
        self,
        author_snowflake: int,
        ctx: SnowflakeContext,
        interaction: discord.Interaction,
        tick: Tick,
    ):
        super().__init__(timeout=120)
        if ctx.member_snowflake is None:
            raise CheckFailure("Member not found.")
        self.__author_snowflake = author_snowflake
        self.__ctx = ctx
        self.__interaction = interaction
        self.__member_snowflake = ctx.member_snowflake
        self.__tick = tick
        self.__groups: dict[str, PermissionGroup] = {}
        self.__scopes: dict[str, GroupScope] = {}
        self.__selected_group: PermissionGroup | None = None
        self.__selected_guild: discord.Guild | None = None
        self.__selected_channel: discord.abc.GuildChannel | None = None
        self.remove_item(self.channel_select)
        self.remove_item(self.guild_select)

    async def interaction_check(self, interaction):
        return interaction.user.id == self.__author_snowflake

    def limit_available_to_top_24_by_member_count(self, available):
        items = []
        items.extend(available)
        items.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
        top_24 = items[:24]
        return top_24

    async def setup(self):
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        self.__groups.clear()
        self.__scopes.clear()
        assigned = await permission_service.resolve_all_assigned_groups(
            permission_state=permission_state,
            member_snowflake=self.__author_snowflake,
        )
        for group, guild_snowflake, channel_snowflake in assigned:
            scope = GroupScope(group=group)
            if guild_snowflake is not None:
                guild = bot.get_guild(guild_snowflake)
                if guild is not None:
                    scope.guilds[guild.id] = guild
            elif guild_snowflake is None and group.scope == PermissionScope.GLOBAL:
                for guild in bot.guilds:
                    for channel in guild.channels:
                        scope.channels[channel.id] = channel
                    scope.guilds[guild.id] = guild
            if channel_snowflake is not None:
                channel = bot.get_channel(channel_snowflake)
                if channel is not None and isinstance(
                    channel,
                    (discord.VoiceChannel, discord.TextChannel, discord.StageChannel),
                ):
                    scope.channels[channel.id] = channel
            self.add_group_scope(scope)
            self.add_selectable_group(group, scope)
        if not self.__groups:
            raise CheckFailure(
                "You do not have sufficient privileges in this channel or server to use this command."
            )
        self._build_group_options(available_groups=list(self.__groups.keys()))

    def add_group_scope(self, scope: GroupScope):
        existing = self.__scopes.setdefault(
            scope.group.alias,
            GroupScope(group=scope.group),
        )
        existing.guilds.update(scope.guilds)
        existing.channels.update(scope.channels)

    def add_selectable_group(self, group: PermissionGroup, scope: GroupScope):
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        ancestors = permission_service.resolve_ancestors(
            group_alias=group.alias,
            groups=permission_state.groups,
        )
        for ancestor_alias in ancestors:
            ancestor = permission_state.groups[ancestor_alias]
            if ancestor.is_sysadmin or ancestor.is_guild_owner or ancestor.default:
                continue
            self.__groups.setdefault(ancestor.alias, ancestor)
            ancestor_scope = GroupScope(group=ancestor)
            ancestor_scope.guilds.update(scope.guilds)
            ancestor_scope.channels.update(scope.channels)
            self.add_group_scope(ancestor_scope)

    def _build_group_options(self, available_groups: list[str]):
        group_options = []
        groups = [self.__groups[alias] for alias in available_groups]
        groups.sort(
            key=lambda group: len(group.ancestors),
            reverse=True,
        )
        group_options.extend(
            [discord.SelectOption(label=g.name, value=g.alias) for g in groups]
        )
        self.group_select.options = group_options

    def _build_channel_options(
        self, available_channels: list[discord.abc.GuildChannel], default: bool
    ):
        limited_channels = self.limit_available_to_top_24_by_member_count(
            available=available_channels
        )
        channel_options = []
        channel_options.extend(
            [
                discord.SelectOption(label=c.name, value=str(c.id))
                for c in limited_channels
                if isinstance(
                    c, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
                )
            ]
        )
        for option in channel_options:
            if option.value.isdigit():
                if int(option.value) == self.__ctx.channel_snowflake and default:
                    option.default = True
        self.channel_select.options = channel_options

    def _build_guild_options(
        self, available_guilds: list[discord.Guild], default: bool
    ):
        limited_guilds = self.limit_available_to_top_24_by_member_count(
            available=available_guilds
        )
        guild_options = []
        guild_options.extend(
            [
                discord.SelectOption(label=g.name, value=str(g.id))
                for g in limited_guilds
            ]
        )
        for option in guild_options:
            if option.value.isdigit():
                if int(option.value) == self.__ctx.guild_snowflake and default:
                    option.default = True
        self.guild_select.options = guild_options

    @discord.ui.select(placeholder="Select a group", options=[])
    async def group_select(self, interaction, select):
        await interaction.response.defer()
        group_alias = select.values[0]
        group = self.__groups[group_alias]
        self.group_select.placeholder = group.name
        self.__selected_group = group
        bot: DiscordBot = DiscordBot.get_instance()
        scope = self.__scopes.get(group_alias, GroupScope(group=group))
        guilds = list(scope.guilds.values())
        channels = list(scope.channels.values())
        self.__selected_guild = None
        self.__selected_channel = None
        if self.guild_select in self.children:
            self.remove_item(self.guild_select)
        if self.channel_select in self.children:
            self.remove_item(self.channel_select)
        required = SCOPE_REQUIREMENTS[group.scope]
        if "guild" in required:
            if not guilds:
                guilds = list(bot.guilds)
            if interaction.guild and interaction.guild.id in {g.id for g in guilds}:
                self.__selected_guild = interaction.guild
            elif len(guilds) == 1:
                self.__selected_guild = guilds[0]
            self._build_guild_options(guilds, default=True)
            self.add_item(self.guild_select)
        if "channel" in required:
            if self.__selected_guild:
                if scope.channels:
                    channels = [
                        c
                        for c in scope.channels.values()
                        if c.guild.id == self.__selected_guild.id
                    ]
                else:
                    channels = list(self.__selected_guild.channels)
            elif not channels:
                channels = [c for g in guilds for c in g.channels]
            if interaction.channel and interaction.channel.id in {
                c.id for c in channels
            }:
                self.__selected_channel = interaction.channel
            elif len(channels) == 1:
                self.__selected_channel = channels[0]
            self._build_channel_options(channels, default=True)
            self.add_item(self.channel_select)
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select a guild",
        options=[discord.SelectOption(label="Select a group first", value=str(None))],
    )
    async def guild_select(self, interaction, select):
        await interaction.response.defer()
        guild_snowflake = int(select.values[0])
        if self.__selected_group:
            scope = self.__scopes.get(
                self.__selected_group.alias, GroupScope(group=self.__selected_group)
            )
            self.__selected_guild = scope.guilds[guild_snowflake]
            self.guild_select.placeholder = self.__selected_guild.name
            channels = [
                channel
                for channel in scope.channels.values()
                if channel.guild.id == guild_snowflake
            ]
            if len(channels) == 1:
                self.__selected_channel = channels[0]
                self.channel_select.disabled = True
            else:
                self.__selected_channel = None
                self.channel_select.disabled = False
            self._build_channel_options(channels, default=False)
            for option in self.guild_select.options:
                option.default = False
            await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select a channel",
        options=[discord.SelectOption(label="Select a guild first", value=str(None))],
    )
    async def channel_select(self, interaction, select):
        await interaction.response.defer()
        channel_snowflake = int(select.values[0])
        if self.__selected_group:
            scope = self.__scopes.get(
                self.__selected_group.alias, GroupScope(group=self.__selected_group)
            )
            self.__selected_channel = scope.channels[channel_snowflake]
            self.channel_select.placeholder = self.__selected_channel.name
            if self.__selected_guild is None:
                self.__selected_guild = self.__selected_channel.guild
            for option in self.channel_select.options:
                option.default = False
            await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        if self.__selected_group is None:
            return await interaction.response.send_message(
                content="Please select all fields.", ephemeral=True
            )
        self.__tick.update_source(interaction=interaction)
        record = None
        bot: DiscordBot = DiscordBot.get_instance()
        database_factory: DatabaseFactory = DatabaseFactory(MODEL)
        permission_state: PermissionState = bot.registry.get(PermissionState)
        member = bot.get_user(self.__member_snowflake)
        if member:
            member_str = member.mention
        else:
            simplified_member = bot.registry.get(MemberState).active.get(
                self.__member_snowflake
            )
            if not simplified_member:
                raise MemberNotFound(str(self.__member_snowflake))
            display_name = simplified_member[0]
            member_str = display_name
        if self.__selected_group and self.__selected_channel and self.__selected_guild:
            group = await permission_service.resolve_effective_group(
                permission_state=permission_state,
                member_snowflake=self.__member_snowflake,
                guild_snowflake=self.__selected_guild.id,
                channel_snowflake=self.__selected_channel.id,
            )
            if group:
                if self.__selected_group.alias in group.ancestors:
                    return await interaction.response.send_message(
                        content=f"You cannot grant {member_str} a group they inherit from.",
                        ephemeral=True,
                    )
            record = await database_factory.select(
                channel_snowflake=self.__selected_channel.id,
                guild_snowflake=self.__selected_guild.id,
                member_snowflake=self.__member_snowflake,
                group_alias=self.__selected_group.alias,
                singular=True,
            )
            entry = PermissionEntry(
                channel_snowflake=self.__selected_channel.id,
                group_alias=self.__selected_group.alias,
                guild_snowflake=self.__selected_guild.id,
                member_snowflake=self.__member_snowflake,
            )
            embed = self.build_grant_embed(
                group=self.__selected_group,
                member_snowflake=self.__member_snowflake,
                channel_snowflake=self.__selected_channel.id,
                guild_snowflake=self.__selected_guild.id,
            )
        elif self.__selected_group and self.__selected_guild:
            group = await permission_service.resolve_effective_group(
                permission_state=permission_state,
                member_snowflake=self.__member_snowflake,
                guild_snowflake=self.__selected_guild.id,
            )
            if group:
                if self.__selected_group.alias in group.ancestors:
                    return await interaction.response.send_message(
                        content=f"You cannot grant {member_str} a group they inherit from.",
                        ephemeral=True,
                    )
            record = await database_factory.select(
                guild_snowflake=self.__selected_guild.id,
                member_snowflake=self.__member_snowflake,
                group_alias=self.__selected_group.alias,
                singular=True,
            )
            entry = PermissionEntry(
                group_alias=self.__selected_group.alias,
                guild_snowflake=self.__selected_guild.id,
                member_snowflake=self.__member_snowflake,
            )
            embed = self.build_grant_embed(
                group=self.__selected_group,
                member_snowflake=self.__member_snowflake,
                guild_snowflake=self.__selected_guild.id,
            )
        else:
            group = await permission_service.resolve_effective_group(
                permission_state=permission_state,
                member_snowflake=self.__member_snowflake,
            )
            if group:
                if self.__selected_group.alias in group.ancestors:
                    return await interaction.response.send_message(
                        content=f"You cannot grant {member_str} a group they inherit from.",
                        ephemeral=True,
                    )
            record = await database_factory.select(
                member_snowflake=self.__member_snowflake,
                group_alias=self.__selected_group.alias,
                singular=True,
            )
            bot.logger.info(record)
            entry = PermissionEntry(
                group_alias=self.__selected_group.alias,
                member_snowflake=self.__member_snowflake,
            )
            embed = self.build_grant_embed(
                group=self.__selected_group,
                member_snowflake=self.__member_snowflake,
            )
        if record:
            return await self.__tick.end(
                warning=f"This member is already apart of this group (`{self.__selected_group.name}`).",
                ephemeral=True,
            )
        await database_factory.create(entry)
        await self.__tick.end(success=embed)
        self.stop()

    def build_grant_embed(
        self,
        group: PermissionGroup,
        member_snowflake: int,
        channel_snowflake: int | None = None,
        guild_snowflake: int | None = None,
    ) -> discord.Embed:
        bot: DiscordBot = DiscordBot.get_instance()
        description = ""
        author = bot.get_user(self.__author_snowflake)
        if author:
            description += f"**Author:** {author.mention}\n"
        member = bot.get_user(member_snowflake)
        if member:
            display_name = member.display_name
            member_str = member.mention
        else:
            simplified_member = bot.registry.get(MemberState).active.get(
                member_snowflake
            )
            if not simplified_member:
                raise commands.MemberNotFound(str(member_snowflake))
            display_name = simplified_member[0]
            member_str = display_name
        description += f"**User:** {member_str}\n"
        if guild_snowflake:
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                raise commands.GuildNotFound(str(guild_snowflake))
            description += f"**Server:** {guild.name}\n"
            if channel_snowflake:
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    raise commands.ChannelNotFound(str(channel_snowflake))
                description += f"**Channel:** {channel.mention}"
        embed = discord.Embed(
            title=f"{emojis.get_random_emoji()} {display_name} has been granted {group.name}",
            description=description,
            color=discord.Color.blue(),
        )
        if member:
            embed.set_thumbnail(url=member.display_avatar.url)
        return embed

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
