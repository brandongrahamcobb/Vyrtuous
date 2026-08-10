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

import discord

from vyrtuous.aliases import role_alias_service, unrole_alias_service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.permissions import PermissionScope
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.role import Role
from vyrtuous.models.duration import DurationObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import (
    CheckFailure,
    GuildNotFound,
    MemberNotFound,
    RoleNotFound,
)
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.tracking import stream_service
from vyrtuous.utils.users import vegan_service


class RoleView(discord.ui.View):
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
        self.__author_snowflake = author_snowflake
        self.__channel_snowflake = ctx.channel_snowflake
        self.__ctx = ctx
        self.__guild_snowflake = ctx.guild_snowflake
        self.__role_snowflake = None
        self.__interaction = interaction
        self.__tick = tick
        self.__model = Role

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
        for group, _, _ in all_assigned_groups:
            if group.scope == PermissionScope.GLOBAL:
                for guild in bot.guilds:
                    for channel in guild.channels:
                        available_channels.add(channel)
                    available_guilds.add(guild)
                break
            elif group.scope == PermissionScope.GUILD:
                for guild in bot.guilds:
                    if guild.id != self.__ctx.guild_snowflake:
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
                break
            elif group.scope == PermissionScope.CHANNEL:
                for guild in bot.guilds:
                    for channel in guild.channels:
                        effective_group = permission_service.resolve_effective_group(
                            permission_state=permission_state,
                            member_snowflake=self.__author_snowflake,
                            channel_snowflake=channel.id,
                            guild_snowflake=guild.id,
                        )
                        if effective_group == group:
                            for channel in guild.channels:
                                available_channels.add(channel)
                            available_guilds.add(guild)
                break
            else:
                raise CheckFailure(
                    "You do not have sufficient privileges in this channel or server to use this command."
                )
        self.available_guilds = available_guilds
        self.available_channels = available_channels
        self._build_guild_options(available_guilds=available_guilds, default=True)
        await self._build_role_options(
            available_guilds=available_guilds, available_channels=available_channels
        )

    def limit_available_to_top_24_by_member_count(self, available):
        items = list(available)
        items.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
        top_24 = items[:24]
        return set(top_24)

    async def _build_role_options(
        self,
        available_guilds: set[discord.Guild],
        available_channels: set[discord.abc.GuildChannel],
    ):
        role_options = []
        database_factory: DatabaseFactory = DatabaseFactory(self.__model)
        for guild in available_guilds:
            for channel in available_channels:
                if channel in guild.channels:
                    records = await database_factory.select(
                        channel_snowflake=channel.id,
                        guild_snowflake=guild.id,
                        singular=False,
                    )
                    for record in records:
                        role = guild.get_role(record.role_snowflake)
                        if role is None:
                            raise RoleNotFound(str(record.role_snowflake))
                        role_options.extend(
                            [discord.SelectOption(label=role.name, value=str(role.id))]
                        )
        if role_options:
            self.role_select.placeholder = "Select a role."
            self.role_select.options = role_options
            self.role_select.disabled = False
        else:
            self.role_select.placeholder = "No roles found."
            self.role_select.options = [
                discord.SelectOption(label="No roles found.", value=str(None))
            ]
            self.role_select.disabled = True

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
        self.guild_select.placeholder = guild.name
        self.__guild_snowflake = guild.id
        await self._build_role_options(
            available_guilds=set([guild]), available_channels=self.available_channels
        )
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select a role",
        options=[],
    )
    async def role_select(self, interaction, select):
        await interaction.response.defer()
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(int(self.__guild_snowflake))
        if guild is None:
            raise GuildNotFound(str(self.__guild_snowflake))
        role = guild.get_role(int(select.values[0]))
        if role is None:
            raise RoleNotFound(select.values[0])
        self.role_select.placeholder = role.name
        self.__role_snowflake = role.id
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        if (
            self.__channel_snowflake is None
            or self.__guild_snowflake is None
            or self.__role_snowflake is None
        ):
            return await interaction.response.send_message(
                content=f"Please select all fields.",
                ephemeral=True,
            )
        self.stop()
        self.__tick.update_source(interaction=interaction)
        await self.__tick.defer()
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=self.__channel_snowflake,
            guild_snowflake=self.__guild_snowflake,
            author_snowflake=self.__author_snowflake,
            member_snowflake=self.__ctx.member_snowflake,
        )
        guild = bot.get_guild(self.__guild_snowflake)
        if guild is None:
            raise GuildNotFound(str(self.__guild_snowflake))
        member = guild.get_member(self.__ctx.member_snowflake)
        if member is None:
            raise MemberNotFound(str(self.__ctx.member_snowflake))
        role = guild.get_role(self.__role_snowflake)
        if role is None:
            raise RoleNotFound(str(self.__role_snowflake))
        if role in member.roles:
            embed = None
            target = None
            is_channel_scope = await unrole_alias_service.disable(
                channel_snowflake=self.__channel_snowflake,
                guild_snowflake=self.__guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                role_snowflake=self.__role_snowflake,
                reason="Revoking role.",
            )
            embed = unrole_alias_service.build_unrole_embed(
                guild_snowflake=self.__guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                role_snowflake=self.__role_snowflake,
            )
            await stream_service.log(
                author_snowflake=self.__author_snowflake,
                channel_snowflake=self.__channel_snowflake,
                display=True,
                duration=DurationObject(number=0, prefix="", sign=1, unit=""),
                guild_snowflake=self.__guild_snowflake,
                identifier="unrole",
                is_channel_scope=is_channel_scope,
                member_snowflake=self.__ctx.member_snowflake,
                message_snowflake=None,
                message_channel_snowflake=None,
                reason="No reason provided.",
                role_snowflake=self.__role_snowflake,
                target=target,
            )
            if embed:
                await self.__tick.end(success=embed)
        else:
            embed = None
            target = None
            is_channel_scope = await role_alias_service.enable(
                channel_snowflake=self.__channel_snowflake,
                guild_snowflake=self.__guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                role_snowflake=self.__role_snowflake,
                reason="Granting role.",
            )
            embed = role_alias_service.build_enrole_embed(
                guild_snowflake=self.__guild_snowflake,
                member_snowflake=self.__ctx.member_snowflake,
                role_snowflake=self.__role_snowflake,
            )
            await stream_service.log(
                author_snowflake=self.__author_snowflake,
                channel_snowflake=self.__channel_snowflake,
                display=True,
                duration=DurationObject(number=0, prefix="", sign=1, unit=""),
                guild_snowflake=self.__guild_snowflake,
                identifier="role",
                is_channel_scope=is_channel_scope,
                member_snowflake=self.__ctx.member_snowflake,
                message_snowflake=None,
                message_channel_snowflake=None,
                reason="No reason provided.",
                role_snowflake=self.__role_snowflake,
                target=target,
            )
            if embed:
                await self.__tick.end(success=embed)

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
