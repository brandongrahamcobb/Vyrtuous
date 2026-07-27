"""!/bin/python3
list_app_commands.py A discord.py cog containing listing commands for the Vyrtuous bot.

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
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.listing import list_bans, list_flags, list_text_mutes, list_voice_mutes
from vyrtuous.models.scope import AppScope, ScopeObject
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.permissions import permission_service


class ListAppCommands(commands.Cog):

    def __init__(
        self,
        bot: DiscordBot,
    ):
        self.__bot = bot

    @app_commands.command(name="bans", description="List bans.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_bans_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.scope.guild"],
                )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning=f"This command must be used in a server channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.guild"],
            )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.channel"],
            )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.member"],
            )
        pages = await list_bans.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @app_commands.command(name="flags", description="List flags.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_flags_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.scope.guild"],
                )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning=f"This command must be used in a server channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.guild"],
            )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.channel"],
            )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.member"],
            )
        pages = await list_flags.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @app_commands.command(name="mutes", description="List mutes.")
    @app_commands.describe(
        target="Specify one of: a channel ID/mention, member ID/mention or server ID.",
        scope="Specify one of: `all`, `auto`, `click`, or `command`.",
        guild="Specify a server ID.",
    )
    async def list_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        scope: app_commands.Transform[ScopeObject | None, AppScope] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.scope.guild"],
                )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning=f"This command must used in a server channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.guild"],
            )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.channel"],
            )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.member"],
            )
        if scope is None:
            mute_type = "all"
        else:
            mute_type = scope.scope
        match mute_type:
            case "all":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=[
                        "command.listing.voice-mutes.auto"
                        "command.listing.voice-mutes.click"
                        "command.listing.voice-mutes.command"
                    ],
                )
            case "auto":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.auto"],
                )
            case "click":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.click"],
                )
            case "command":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.command"],
                )
        pages = await list_voice_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj, mute_type=mute_type
        )
        return await tick.end(success=pages)

    @app_commands.command(name="summary", description="List user moderation.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
    )
    async def list_moderation_summary_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        scope: app_commands.Transform[ScopeObject | None, AppScope] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.listing.scope.member"],
        )
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.scope.guild"],
                )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        if scope is None:
            mute_type = "all"
        else:
            mute_type = scope.scope
        match mute_type:
            case "all":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=[
                        "command.listing.voice-mutes.auto"
                        "command.listing.voice-mutes.click"
                        "command.listing.voice-mutes.command"
                    ],
                )
            case "auto":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.auto"],
                )
            case "click":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.click"],
                )
            case "command":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.command"],
                )
        pages: list[discord.Embed] = []
        services = []
        services.append(list_bans)
        services.append(list_flags)
        services.append(list_text_mutes)
        for service in services:
            summary_pages = await service.build_pages(
                guild_snowflake=guild_snowflake,
                obj=member_snowflake,
            )
            if isinstance(summary_pages, list):
                for page in summary_pages:
                    if isinstance(page, discord.Embed):
                        pages.append(page)
        summary_pages = await list_voice_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=member_snowflake, mute_type=mute_type
        )
        if isinstance(summary_pages, list):
            for page in summary_pages:
                if isinstance(page, discord.Embed):
                    pages.append(page)
        if not pages:
            return await tick.end(success="No infractions found")
        return await tick.end(success=pages)

    @app_commands.command(name="tmutes", description="List text-mutes.")
    @app_commands.describe(
        target="Specify a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_text_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.scope.guild"],
                )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning=f"This command must be used in a server channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.guild"],
            )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.channel"],
            )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.member"],
            )
        pages = await list_text_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    async def cog_app_command_error(self, interaction, error):
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))


async def setup(bot: DiscordBot):
    await bot.add_cog(ListAppCommands(bot=bot))
