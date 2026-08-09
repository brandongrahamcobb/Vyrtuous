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

from datetime import timedelta

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.inc.helpers import DISCORD_COGS, DISCORD_COGS_CLASSES, PATH_LOG
from vyrtuous.listing import (list_autoassign_roles, list_automute_channels,
                              list_bans, list_caps, list_flags, list_heroes,
                              list_intents, list_overwrites, list_streams,
                              list_text_mutes, list_vegans,
                              list_video_channels, list_voice_mutes)
from vyrtuous.models.metadata import metadata
from vyrtuous.models.scope import AppScope, ScopeObject
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.statistics import system_monitoring_service
from vyrtuous.utils.users import vegan_service
from vyrtuous.view.groups_view import GroupsView


class InfoAppCommands(commands.Cog):

    def __init__(
        self,
        bot: DiscordBot,
    ):
        self.__bot = bot

    async def cog_app_command_error(self, interaction, error):
        self.__bot.logger.info(str(error))
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error), ephemeral=True)

    @metadata(permission="command.info.autoassigns")
    @app_commands.command(name="autoassigns", description="List group autoassignment roles.")
    @app_commands.describe(guild="Specify a server ID.")
    async def list_autoassignment_roles_app_command(
        self,
        interaction: discord.Interaction,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(warning="This command must target a valid server.", ephemeral=True)
        if interaction.channel is None:
            return await tick.end(warning="This command must target a valid server channel.", ephemeral=True)
        else:
            channel_snowflake = interaction.channel.id
        if guild is None:
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        await tick.defer()
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.scope.guild", "command.info.autoassigns"],
        )
        if interaction.guild.id != guild_snowflake:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["other_guilds"],
            )
        pages = await list_autoassign_roles.build_pages(
            guild_snowflake=guild_snowflake
        )
        return await tick.end(success=pages)

    @metadata(permission="command.info.automutes")
    @app_commands.command(name="automutes", description="List automute channels.")
    @app_commands.describe(
        target="Specify a channel ID/mention or server ID.",
    )
    async def list_automute_channels_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(warning="This command must target a valid server.", ephemeral=True)
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(warning="This command must target a valid server channel.", ephemeral=True)
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
            obj = interaction.guild
        else:
            obj = target.target
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.automutes"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.automutes"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=obj.id,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_automute_channels.build_pages(obj=obj)
        return await tick.end(success=pages)

    @metadata(permission="command.info.bans")
    @app_commands.command(name="bans", description="List bans.")
    @app_commands.describe(
        target="Specify one of: a channel ID/mention, member ID/mention or server ID.",
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
                return await tick.end(warning="This command must target a valid server.", ephemeral=True)
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must target a valid server.", ephemeral=True
                    )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        if interaction.channel is None:
            return await tick.end(
                warning=f"This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
            obj = interaction.channel
        else:
            obj = target.target
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.bans"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.bans"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=obj.id,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.bans"],
            )
            if (
                guild
                and isinstance(guild.target, discord.Guild)
                and guild.target.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.target.id,
                    requested=["other_guilds"],
                )
        pages = await list_bans.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @metadata(permission="command.info.blacklists")
    @app_commands.command(name="blacklists", description="List blacklists.")
    @app_commands.describe(
        target="Specify a channel ID/mention, member ID/mention or server ID."
    )
    async def list_blacklists_app_command(
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
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
           obj = interaction.channel
        else:
            obj = target.target
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.blacklists"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.blacklists"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=obj.id,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.blacklists"],
            )
            if (
                guild
                and isinstance(guild.target, discord.Guild)
                and guild.target.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.target.id,
                    requested=["other_guilds"],
                )
        pages = await list_bans.build_blacklist_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @metadata(permission="command.info.caps")
    @app_commands.command(name="caps", description="List caps.")
    @app_commands.describe(
        target="Specify a channel ID/mention or server ID.",
    )
    async def list_caps_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(
                warning="This command must target a valid server.", ephemeral=True
            )
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
           obj = interaction.channel
        else:
            obj = target.target if target.target != "all" else interaction.guild
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.caps"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.caps"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_caps.build_pages(obj=obj)
        return await tick.end(success=pages)

    @metadata(permission="command.info.cogs")
    @app_commands.command(name="cogs", description="Lists cogs.")
    async def list_cogs_app_command(
        self, interaction: discord.Interaction
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.", ephemeral=True)
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        await tick.defer()
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.cogs"],
        )
        loaded, not_loaded = [], []
        embed = discord.Embed(
            title=f"{emojis.get_random_emoji()} Cogs for {interaction.guild.me.name}",
            color=discord.Color.blurple(),
        )
        for representation, cog in zip(
            sorted(DISCORD_COGS), sorted(DISCORD_COGS_CLASSES)
        ):
            if cog in self.__bot.cogs:
                loaded.append(representation)
            else:
                not_loaded.append(representation)
        if loaded:
            embed.add_field(name="Loaded", value="\n".join(loaded), inline=False)
        if not_loaded:
            embed.add_field(
                name="Not Loaded", value="\n".join(not_loaded), inline=False
            )
        if not loaded and not not_loaded:
            embed.add_field(name="No cogs available.", value=None, inline=False)
        return await tick.end(success=embed)

    @metadata(permission="command.info.debug")
    @app_commands.command(
        name="debug", description="Shows the last `n` number of logging."
    )
    @app_commands.describe(lines="Specify the number of lines.")
    async def debug_app_command(
        self, interaction: discord.Interaction, lines: int | None = None
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.", ephemeral=True)
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if lines is None:
            lines = 25
        if 25 < lines <= 0:
            return await tick.end(warning="Lines must be between than 1 and 25", ephemeral=True)
        await tick.defer()
        await permission_service.has_permissions(
            permission_state=permission_state,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=interaction.user.id,
            requested=["command.info.debug"],
        )
        try:
            with open(PATH_LOG, "r") as f:
                content = f.readlines()[-lines:]
                content = [line.split(" - ", 3)[-1] for line in content]
        except FileNotFoundError:
            return await tick.end(warning="Log file not found")
        output = "".join(content)
        if len(output) > 1900:
            output = output[-1900:]
        return await tick.end(success=f"```log\n{output}\n```")

    @metadata(permission="command.info.flags")
    @app_commands.command(name="flags", description="List flags.")
    @app_commands.describe(
        target="Specify a channel ID/mention, member ID/mention or server ID.",
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
                return await tick.end(warning="This command must be used in a server.", ephemeral=True)
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must target a valid server.", ephemeral=True
                    )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        if interaction.channel is None:
            return await tick.end(
                warning=f"This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
           obj = interaction.channel
        else:
            obj = target.target
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.flags"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.flags"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.flags"],
            )
            if (
                guild
                and isinstance(guild.target, discord.Guild)
                and guild.target.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.target.id,
                    requested=["other_guilds"],
                )
        pages = await list_flags.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @metadata(permission="command.info.groups")
    @app_commands.command(name="groups", description="List groups.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def grant_group_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.", ephemeral=True)
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if member is None:
            member_snowflake = None
        elif isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.", ephemeral=True)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.groups"],
        )
        view = GroupsView(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            interaction=interaction,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the group", view=view, ephemeral=True
        )

    @metadata(permission="command.info.heroes")
    @app_commands.command(name="heroes", description="List heroes.")
    @app_commands.describe(
        target="Specify a member ID/mention or server ID.", guild="Specify a server ID."
    )
    async def list_heroes_app_command(
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
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
           obj = interaction.guild
        else:
            obj = target.target
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.heroes"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.heroes"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.heroes"],
            )
            if (
                guild
                and isinstance(guild.target, discord.Guild)
                and guild.target.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.target.id,
                    requested=["other_guilds"],
                )
        pages = await list_heroes.build_pages(
            guild_snowflake=guild_snowflake,
            obj=obj,
        )
        return await tick.end(success=pages)

    @metadata(permission="command.info.intents")
    @app_commands.command(name="intents", description="View intents.")
    @app_commands.describe(
        target="Specify a channel ID/mention or server ID.",
    )
    async def list_intents_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(
                warning="This command must be used in a server.", ephemeral=True
            )
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
            obj = interaction.channel
        else:
            obj = target.target
        await tick.defer()
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.intents"],
        )
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["command.info.intents"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.intents"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_intents.build_pages(
            obj=obj,
        )
        return await tick.end(success=pages)

    @metadata(permission="command.info.members")
    @app_commands.command(name="members", description="List role members.")
    @app_commands.describe(
        role="Specify a role ID/mention.", guild="Specify a server ID."
    )
    async def list_role_members_app_command(
        self,
        interaction: discord.Interaction,
        role: app_commands.Transform[TargetObject, AppTarget],
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(
                warning="This command must target a valid server.", ephemeral=True
            )
        if guild is None:
            guild_snowflake = interaction.guild.id
            members = interaction.guild.members
        else:
            if isinstance(guild.target, discord.Guild):
                members = guild.target.members
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if isinstance(role.target, discord.Role):
            role_name = role.target.name
            color = (
                role.target.color
                if role.target.color.value
                else discord.Color.blurple()
            )
        else:
            return await tick.end(warning="This command must target a valid role.", ephemeral=True)
        await tick.defer()
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.members", "command.info.scope.role"],
        )
        if interaction.guild.id != guild_snowflake:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["other_guilds"],
            )
        embeds = []
        members = [member for member in members if role in member.roles]
        chunk_size = 12
        for index in range(0, len(members), chunk_size):
            chunk = members[index : index + chunk_size]
            description = "\n".join(
                f"{position + 1}. {member.mention} ({member.id})"
                for position, member in enumerate(chunk, start=index)
            )
            embed = discord.Embed(
                title=f"{role_name} Members",
                description=description or "No members found.",
                color=color,
            )
            embeds.append(embed)
        if not embeds:
            embed = discord.Embed(
                title=f"{role_name} Members",
                description="No members found.",
                color=color,
            )
            embeds.append(embed)
        return await tick.end(success=embeds)



    @metadata(permission="command.info.voice-mutes")
    @app_commands.command(name="mutes", description="List mutes.")
    @app_commands.describe(
        target="Specify one of: a channel ID/mention, member ID/mention or server ID.",
        scope="Specify one of: `all`, `auto`, `click`, `command` or `server`.",
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
                return await tick.end(warning="This command must target a valid server.", ephemeral=True)
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must target a valid server.", ephemeral=True
                    )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        if interaction.channel is None:
            return await tick.end(
                warning=f"This command must used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
           obj = interaction.channel
        else:
            obj = target.target
        if scope is None:
            mute_type = "command"
        else:
            mute_type = scope.scope
            if mute_type == "server":
                obj = interaction.guild
        async def has_mute_permissions(permission_state: PermissionState, channel_snowflake: int, guild_snowflake: int, member_snowflake: int, mute_type: str):
            match mute_type:
                case "all":
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=member_snowflake,
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        requested=[
                            "command.info.voice-mutes.auto",
                            "command.info.voice-mutes.click",
                            "command.info.voice-mutes.command",
                            "command.info.voice-mutes.server",
                        ],
                    )
                case "auto":
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=member_snowflake,
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        requested=["command.info.voice-mutes.auto"],
                    )
                case "click":
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=member_snowflake,
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        requested=["command.info.voice-mutes.click"],
                    )
                case "command":
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=member_snowflake,
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        requested=["command.info.voice-mutes.command"],
                    )
                case "server":
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=member_snowflake,
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        requested=["command.info.voice-mutes.server"],
                    )
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
            await has_mute_permissions(permission_state=permission_state, channel_snowflake=channel_snowflake, guild_snowflake=obj.id, member_snowflake=interaction.user.id, mute_type=mute_type)
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
            await has_mute_permissions(permission_state=permission_state, channel_snowflake=obj.id, guild_snowflake=obj.guild.id, member_snowflake=interaction.user.id, mute_type=mute_type)
        elif isinstance(obj, discord.Member):
            if (
                guild
                and isinstance(guild.target, discord.Guild)
                and guild.target.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.target.id,
                    requested=["command.info.scope.member"],
                )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.target.id,
                    requested=["other_guilds"],
                )
                await has_mute_permissions(permission_state=permission_state, channel_snowflake=channel_snowflake, guild_snowflake=guild.target.id, member_snowflake=interaction.user.id, mute_type=mute_type)
            else:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                    requested=["command.info.scope.member"],
                )
                await has_mute_permissions(permission_state=permission_state, channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=interaction.user.id, mute_type=mute_type)
        pages = await list_voice_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj, mute_type=mute_type
        )
        return await tick.end(success=pages)

    @metadata(permission="command.info.overwrites")
    @app_commands.command(name="overwrites", description="Overwrite stats.")
    @app_commands.describe(
        target="Specify a channel ID/mention or server ID.",
    )
    async def list_overwrites_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(
                warning="This command must be used in a server.", ephemeral=True
            )
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
           obj = interaction.channel
        else:
            obj = target.target
        await tick.defer()
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.overwrites"],
        )
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.overwrites"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.overwrites"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        embed = list_overwrites.build_embed(obj=obj)
        return await tick.end(success=embed)

    @metadata(permission="command.info.roleid")
    @app_commands.command(name="roleid", description="Get role by name.")
    @app_commands.describe(
        role_name="Specify a role name.", guild="Specify a server ID."
    )
    async def get_role_snowflake_app_command(
        self,
        interaction: discord.Interaction,
        role_name: str,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(
                warning="This command must be used in a server.", ephemeral=True
            )
        if guild is None:
           guild_obj = interaction.guild
        else:
            if isinstance(guild.target, discord.Guild):
                guild_obj = guild.target
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        await tick.defer()
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_obj.id,
            requested=["command.info.roleid", "command.info.scope.role"],
        )
        if interaction.guild.id != guild_obj.id:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_obj.id,
                requested=["other_guilds"],
            )
        role = discord.utils.get(guild_obj.roles, name=role_name)
        if role:
            return await tick.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await tick.end(
                warning=f"No role named `{role_name}` found in this server."
            )

    @metadata(permission="command.info.stats")
    @app_commands.command(name="stats", description="Lists stats.")
    async def list_stats_app_command(
        self, interaction: discord.Interaction
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.", ephemeral=True)
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        await tick.defer()
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.stats"],
        )
        embed = discord.Embed(title="Statistics")
        cpu_usage = await system_monitoring_service.calculate_cpu_usage()
        embed.add_field(name="CPU %", value=cpu_usage, inline=True)
        with open("/sys/fs/cgroup/memory.current", "r") as file:
            bits = int(file.read())
            memory_usage = round((bits / 1024) / 1024, 0)
            embed.add_field(name="RAM", value=f"{memory_usage} MB", inline=True)
        with open("/proc/uptime", "r") as file:
            content = file.readline()
            fields = content.split()
            uptime_seconds = float(fields[0])
            time = timedelta(seconds=uptime_seconds)
        embed.add_field(name="Uptime", value=f"{str(time)}", inline=True)
        rx_usage = await system_monitoring_service.calculate_rx_usage()
        embed.add_field(name="RX MB", value=f"{rx_usage} MB", inline=True)
        tx_usage = await system_monitoring_service.calculate_tx_usage()
        embed.add_field(name="TX MB", value=f"{tx_usage} MB", inline=True)
        number_of_servers = len(self.__bot.guilds)
        embed.add_field(name="Servers", value=number_of_servers, inline=True)
        return await tick.end(success=embed)

    @metadata(permission="command.info.streams")
    @app_commands.command(name="streams", description="List streaming routes.")
    @app_commands.describe(
        target="Specify one of: a channel ID/mention or server ID",
    )
    async def list_streaming_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(
                warning="This command must be used in a server.", ephemeral=True
            )
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
            obj = interaction.guild
        else:
            obj = target.target
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.streams"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.streams"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_streams.build_pages(obj=obj)
        return await tick.end(success=pages)

    @metadata(permission="command.info.summary")
    @app_commands.command(name="summary", description="List user moderation.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
        scope="Specify `all`, `auto`, `click` or `command`",
        guild="Specify a server ID."
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
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.", ephemeral=True)
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.", ephemeral=True)
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server.", ephemeral=True
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
                    requested=["command.info.scope.guild"],
                )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.", ephemeral=True)
        if scope is None:
            mute_type = "all"
        else:
            mute_type = scope.scope
        async def has_mute_permissions(permission_state: PermissionState, channel_snowflake: int, guild_snowflake: int, member_snowflake: int, mute_type: str):
            match mute_type:
                case "all":
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=member_snowflake,
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        requested=[
                            "command.info.voice-mutes.auto",
                            "command.info.voice-mutes.click",
                            "command.info.voice-mutes.command",
                        ],
                    )
                case "auto":
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=member_snowflake,
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        requested=["command.info.voice-mutes.auto"],
                    )
                case "click":
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=member_snowflake,
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        requested=["command.info.voice-mutes.click"],
                    )
                case "command":
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=member_snowflake,
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        requested=["command.info.voice-mutes.command"],
                    )
        await tick.defer()
        if (
            guild
            and isinstance(guild.target, discord.Guild)
            and guild.target.id != guild_snowflake
        ):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild.target.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild.target.id,
                requested=["command.info.bans", "command.info.flags", "command.info.text-mutes"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild.target.id,
                requested=["other_guilds"],
            )
            await has_mute_permissions(permission_state=permission_state, channel_snowflake=channel_snowflake, guild_snowflake=guild.target.id, member_snowflake=interaction.user.id, mute_type=mute_type)
        else:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["command.info.bans", "command.info.flags", "command.info.text-mutes"],
            )
            await has_mute_permissions(permission_state=permission_state, channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=interaction.user.id, mute_type=mute_type)
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

    @metadata(permission="command.info.survey")
    @app_commands.command(name="survey", description="Show member permissions.")
    @app_commands.describe(
        channel="Specify a channel ID/mention.",
    )
    async def survey_app_command(
        self,
        interaction: discord.Interaction,
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.", ephemeral=True)
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        if channel is None:
            if not isinstance(interaction.channel, (discord.VoiceChannel, discord.StageChannel)):
                return await tick.end(
                    warning="This command must be used in a voice channel or stage channel.", ephemeral=True
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.", ephemeral=True)
        await tick.defer()
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.survey", "command.info.scope.channel"],
        )
        if interaction.guild.id != guild_snowflake:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["other_guilds"],
            )
        pages = await permission_service.survey(
            permission_state=permission_state,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
        )
        return await tick.end(success=pages)

    @metadata(permission="command.info.text-mutes")
    @app_commands.command(name="tmutes", description="List text-mutes.")
    @app_commands.describe(
        target="Specify a channel ID/mention, member ID/mention or server ID.",
        guild="Speficy a server ID."
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
                return await tick.end(warning="This command must be used in a server.", ephemeral=True)
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server.", ephemeral=True
                    )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        if interaction.channel is None:
            return await tick.end(
                warning=f"This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
           obj = interaction.channel
        else:
            obj = target.target
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.text-mutes"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.text-mutes"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.text-mutes"],
            )
            if (
                guild
                and isinstance(guild.target, discord.Guild)
                and guild.target.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.target.id,
                    requested=["other_guilds"],
                )
        pages = await list_text_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @metadata(permission="command.info.vegans")
    @app_commands.command(name="vegans", description="List new vegans.")
    @app_commands.describe(
        target="Specify a member ID/mention or server ID.",
        guild="Specify a server ID."
    )
    async def list_new_vegans_app_command(
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
                return await tick.end(
                    warning="This command must be used in a server.", ephemeral=True
                )
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
           obj = interaction.guild
        else:
            obj = target.target if target.target != "all" else interaction.guild
        if not vegan_service.is_vegan(
            guild_snowflake=guild_snowflake, member_snowflake=interaction.user.id
        ):
            return await tick.end(
                warning="You have insufficient privileges to do that (`command.info.vegans`).", ephemeral=True
            )
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.vegans"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.info.vegans"],
            )
            if (
                guild
                and isinstance(guild.target, discord.Guild)
                and guild.target.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.target.id,
                    requested=["other_guilds"],
                )
        pages = await list_vegans.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @metadata(permission="command.info.video-channels")
    @app_commands.command(
        name="video-only-channels",
        description="List video channels.",
    )
    @app_commands.describe(target="Specify one of: a channel ID/mention or server ID.")
    async def list_video_channels_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(
                warning="This command must be used in a server.", ephemeral=True
            )
        else:
            guild_snowflake = interaction.guild.id
        if interaction.channel is None:
            return await tick.end(
                warning="This command must used in a server channel.", ephemeral=True
            )
        else:
            channel_snowflake = interaction.channel.id
        if target is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
            obj = interaction.guild
        else:
            obj = target.target
        await tick.defer()
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.video-channels"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.video-channels"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_video_channels.build_pages(obj=obj)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(InfoAppCommands(bot=bot))
