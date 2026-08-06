"""!/bin/python3
sysadmin_text_commands.py A discord.py cog containing sysadmin commands for the Vyrtuous bot.

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
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.inc.helpers import DISCORD_COGS, DISCORD_COGS_CLASSES, PATH_LOG
from vyrtuous.listing import (list_aliases, list_autoassign_roles,
                              list_automute_channels, list_bans, list_caps,
                              list_flags, list_heroes, list_intents,
                              list_overwrites, list_streams, list_text_mutes,
                              list_vegans, list_video_channels,
                              list_voice_mutes)
from vyrtuous.models.metadata import metadata
from vyrtuous.models.scope import Scope, ScopeObject
from vyrtuous.models.target import Target, TargetObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.statistics import system_monitoring_service


class InfoTextCommands(commands.Cog):

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot


    @commands.command(name="autoassigns", help="List group autoassignment roles.")
    @metadata(permission="command.info.autoassigns")
    async def list_autoassignment_roles_text_command(
        self,
        ctx: commands.Context,
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if ctx.channel is None:
            return await tick.end(warning="This command must target a valid server channel.")
        else:
            channel_snowflake = ctx.channel.id
        if guild is None:
            guild_snowflake = ctx.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.scope.guild", "command.info.autoassigns"],
        )
        if ctx.guild.id != guild_snowflake:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["other_guilds"],
            )
        pages = await list_autoassign_roles.build_pages(
            guild_snowflake=guild_snowflake
        )
        return await tick.end(success=pages)

    @commands.command(name="automutes", help="List automute channels.")
    @metadata(permission="command.info.automutes")
    async def list_automute_channels_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must target a valid server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
            obj = ctx.guild
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.automutes"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.automutes"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=obj.id,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_automute_channels.build_pages(obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="bans", help="List bans.")
    @metadata(permission="command.info.bans")
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify one of: channel ID/mention, member ID/mention or server ID.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = ctx.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if ctx.guild is None:
                    return await tick.end(
                        warning="This command must target a valid server."
                    )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if ctx.channel is None:
            return await tick.end(
                warning=f"This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
            obj = ctx.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.bans"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.bans"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=obj.id,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.bans"],
            )
            if (
                guild
                and isinstance(guild, discord.Guild)
                and guild.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_bans.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="blacklists", help="List blacklists.")
    @metadata(permission="command.info.blacklists")
    async def list_blacklists_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention, member ID/mention or server ID.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = ctx.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
            obj = ctx.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.blacklists"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.blacklists"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=obj.id,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.blacklists"],
            )
            if (
                guild
                and isinstance(guild, discord.Guild)
                and guild.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_bans.build_blacklist_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @commands.command(name="caps", help="List caps.")
    @metadata(permission="command.info.caps")
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(
                warning="This command must target a valid server."
            )
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
            obj = ctx.guild
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.caps"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.caps"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_caps.build_pages(obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="aliases", help="List aliases.")
    @metadata(permission="command.info.aliases")
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify channel ID/mention, or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.aliases"],
        )
        if target is None:
            obj = ctx.channel
        else:
            obj = target.target if target.target != "all" else ctx.guild
        pages = await list_aliases.build_pages(obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="cogs", help="Lists cogs.")
    @metadata(permission="command.info.cogs")
    async def list_cogs_text_command(self, ctx: commands.Context) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.cogs"],
        )
        loaded, not_loaded = [], []
        embed = discord.Embed(
            title=f"{emojis.get_random_emoji()} Cogs for {ctx.guild.me.name}",
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

    @commands.command(name="debug", help="Shows the last `n` number of logging.")
    @metadata(permission="command.info.debug")
    async def debug_text_command(
        self,
        ctx,
        lines: int = commands.parameter(
            default=25, description="Specify the number of lines."
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        await permission_service.has_permissions(
            permission_state=permission_state,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=ctx.author.id,
            requested=["command.info.debug"],
        )
        if 25 < lines <= 0:
            return await tick.end(warning="Lines must be between 0 and 25")
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

    @commands.command(name="flags", help="List flags.")
    @metadata(permission="command.info.flags")
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify one of: channel ID/mention, member ID/mention, or server ID.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = ctx.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if ctx.guild is None:
                    return await tick.end(
                        warning="This command must target a valid server."
                    )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if ctx.channel is None:
            return await tick.end(
                warning=f"This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
           obj = ctx.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.flags"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.flags"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.flags"],
            )
            if (
                guild
                and isinstance(guild, discord.Guild)
                and guild.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_flags.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="heroes", help="List heroes.")
    @metadata(permission="command.info.heroes")
    async def list_heroes_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a member mention/ID or server ID.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = ctx.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
           obj = ctx.guild
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.heroes"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.heroes"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.heroes"],
            )
            if (
                guild
                and isinstance(guild, discord.Guild)
                and guild.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_heroes.build_pages(
            guild_snowflake=guild_snowflake,
            obj=obj,
        )
        return await tick.end(success=pages)

    @commands.command(name="intents", help="View intents.")
    @metadata(permission="command.info.intents")
    async def list_intents_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None= commands.parameter(
            converter=Target,
            default=None,
            description="Specify one of: channel ID/mention or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(
                warning="This command must be used in a server."
            )
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            obj = ctx.channel
        else:
            obj = target.target
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.intents"],
        )
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["command.info.intents"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.intents"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_intents.build_pages(
            obj=obj,
        )
        return await tick.end(success=pages)

    @commands.command(name="mutes", help="List mutes.")
    @metadata(permission="command.info.voice-mutes")
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention, member ID/mention, or server ID.",
        ),
        scope: ScopeObject | None = commands.parameter(
            converter=Scope,
            default=None,
            description="Specify `all`, `click` or `command`.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(warning="This command must target a valid server.")
            guild_snowflake = ctx.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if ctx.guild is None:
                    return await tick.end(
                        warning="This command must target a valid server."
                    )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if ctx.channel is None:
            return await tick.end(
                warning=f"This command must used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
           obj = ctx.channel
        else:
            obj = target.target
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
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
            await has_mute_permissions(permission_state=permission_state, channel_snowflake=channel_snowflake, guild_snowflake=obj.id, member_snowflake=ctx.author.id, mute_type=mute_type)
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
            await has_mute_permissions(permission_state=permission_state, channel_snowflake=obj.id, guild_snowflake=obj.guild.id, member_snowflake=ctx.author.id, mute_type=mute_type)
        elif isinstance(obj, discord.Member):
            if (
                guild
                and isinstance(guild, discord.Guild)
                and guild.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.id,
                    requested=["command.info.scope.member"],
                )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.id,
                    requested=["other_guilds"],
                )
                await has_mute_permissions(permission_state=permission_state, channel_snowflake=channel_snowflake, guild_snowflake=guild.id, member_snowflake=ctx.author.id, mute_type=mute_type)
            else:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                    requested=["command.info.scope.member"],
                )
                await has_mute_permissions(permission_state=permission_state, channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=ctx.author.id, mute_type=mute_type)
        pages = await list_voice_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj, mute_type=mute_type
        )
        return await tick.end(success=pages)

    @commands.command(name="overwrites", help="Overwrite stats.")
    @metadata(permission="command.info.overwrites")
    async def list_overwrites_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(
                warning="This command must be used in a server."
            )
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
           obj = ctx.channel
        else:
            obj = target.target
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.overwrites"],
        )
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.overwrites"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.overwrites"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        embed = list_overwrites.build_embed(obj=obj)
        return await tick.end(success=embed)

    @commands.command(name="roleid", help="Get role ID by name.")
    @metadata(permission="command.info.roleid")
    async def get_role_snowflake_text_command(
        self,
        ctx: commands.Context,
        role_name: str,
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(
                warning="This command must target a valid server."
            )
        if guild is None:
           guild_obj = ctx.guild
        else:
            if isinstance(guild.target, discord.Guild):
                guild_obj = guild.target
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_obj.id,
            requested=["command.info.roleid", "command.info.scope.role"],
        )
        if ctx.guild.id != guild_obj.id:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
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

    @commands.command(name="members", help="List role members.")
    @metadata(permission="command.info.members")
    async def list_roles_text_command(
        self,
        ctx: commands.Context,
        role: TargetObject = commands.parameter(
            converter=Target,
            description="Specify a role ID/mention.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(
                warning="This command must target a valid server."
            )
        if guild is None:
            guild_snowflake = ctx.guild.id
            members = ctx.guild.members
        else:
            if isinstance(guild.target, discord.Guild):
                members = guild.target.members
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if isinstance(role.target, discord.Role):
            role_name = role.target.name
            color = (
                role.target.color
                if role.target.color.value
                else discord.Color.blurple()
            )
        else:
            return await tick.end(warning="This command must target a valid role.")
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.members", "command.info.scope.role"],
        )
        if ctx.guild.id != guild_snowflake:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
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

    @commands.command(name="stats", help="Lists stats.")
    @metadata(permission="command.info.stats")
    async def list_stats_text_command(self, ctx: commands.Context) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
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

    @commands.command(name="streams", help="List streaming routes.")
    @metadata(permission="command.info.streams")
    async def list_streaming_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention or a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(
                warning="This command must be used in a server."
            )
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            obj = ctx.guild
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.streams"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.streams"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_streams.build_pages(obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="summary", help="List user moderation.")
    @metadata(permission="command.info.summary")
    async def list_moderation_summary_text_command(
        self,
        ctx: commands.Context,
        member: TargetObject = commands.parameter(
            converter=Target,
            description="Specify a member ID/mention.",
        ),
        scope: ScopeObject | None = commands.parameter(
            converter=Scope,
            default=None,
            description="Specify one of `all`, `click` or `command`.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            requested=["command.info.scope.member"],
        )
        if guild is None:
            if ctx.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = ctx.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if ctx.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != ctx.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=ctx.author.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    requested=["command.info.scope.guild"],
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
        if (
            guild
            and isinstance(guild, discord.Guild)
            and guild.id != guild_snowflake
        ):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild.id,
                requested=["command.info.bans", "command.info.flags", "command.info.text-mutes"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild.id,
                requested=["other_guilds"],
            )
            await has_mute_permissions(permission_state=permission_state, channel_snowflake=channel_snowflake, guild_snowflake=guild.id, member_snowflake=ctx.author.id, mute_type=mute_type)
        else:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["command.info.bans", "command.info.flags", "command.info.text-mutes"],
            )
            await has_mute_permissions(permission_state=permission_state, channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=ctx.author.id, mute_type=mute_type)
        pages: list[discord.Embed] = []
        services = []
        services.append(list_bans)
        services.append(list_flags)
        services.append(list_text_mutes)
        for service in services:
            summary_pages = await service.build_pages(
                guild_snowflake=guild_snowflake, obj=member_snowflake
            )
            if isinstance(summary_pages, list):
                for page in summary_pages:
                    if isinstance(page, discord.Embed):
                        pages.append(page)
        if scope is None:
            mute_type = "all"
        else:
            mute_type = scope.scope
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

    @commands.command(name="survey", help="Show member permissions.")
    @metadata(permission="command.info.survey")
    async def survey_text_command(
        self,
        ctx: commands.Context,
        channel: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(
                warning="This command must be used in a server."
            )
        if ctx.channel is None:
            return await tick.end(
                warning="This command must used in a server channel."
            )
        if channel is None:
            if not isinstance(ctx.channel, (discord.VoiceChannel, discord.StageChannel)):
                return await tick.end(
                    warning="This command must be used in a voice channel or stage channel."
                )
            channel_snowflake = ctx.channel.id
            guild_snowflake = ctx.guild.id
        elif isinstance(channel.target, discord.abc.GuildChannel):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(
                warning="This command must target a valid channel."
            )
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.info.survey", "command.info.scope.channel"],
        )
        if ctx.guild.id != guild_snowflake:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["other_guilds"],
            )
        pages = await permission_service.survey(
            permission_state=permission_state,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake
        )
        return await tick.end(success=pages)

    @commands.command(name="tmutes", help="List text-mutes.")
    @metadata(permission="command.info.text-mutes")
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention, member ID/mention or server ID.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = ctx.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if ctx.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if ctx.channel is None:
            return await tick.end(
                warning=f"This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
           obj = ctx.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.text-mutes"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.text-mutes"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.text-mutes"],
            )
            if (
                guild
                and isinstance(guild, discord.Guild)
                and guild.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_text_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @commands.command(name="vegans", help="List new vegans.")
    @metadata(permission="command.info.vegans")
    async def list_new_vegans_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention, member ID/mention, or server ID.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify one a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = ctx.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
           obj = ctx.guild
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.vegans"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.scope.member"],
            )
            await permission_service.has_permissions_at_all(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                requested=["command.info.vegans"],
            )
            if (
                guild
                and isinstance(guild, discord.Guild)
                and guild.id != guild_snowflake
            ):
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_vegans.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @commands.command(
        name="video-only-channels",
        help="List video channels.",
    )
    @metadata(permission="command.info.video-channels")
    async def list_video_channels_text_command(
        self,
        ctx: commands.Context,
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(
                warning="This command must be used in a server."
            )
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        if target is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            obj = ctx.guild
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.scope.guild"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=obj.id,
                requested=["command.info.video-channels"],
            )
            if obj.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.id,
                    requested=["other_guilds"],
                )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.scope.channel"],
            )
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=obj.id,
                guild_snowflake=obj.guild.id,
                requested=["command.info.video-channels"],
            )
            if obj.guild.id != guild_snowflake:
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=obj.guild.id,
                    requested=["other_guilds"],
                )
        pages = await list_video_channels.build_pages(obj=obj)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(InfoTextCommands(bot=bot))
