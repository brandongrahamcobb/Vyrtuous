"""!/bin/python3
admin_text_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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
from vyrtuous.inc.helpers import PATH_LOG, at_home
from vyrtuous.listing import (
    list_administrator_roles,
    list_automute_channels,
    list_caps,
    list_permissions,
    list_streams,
    list_video_channels,
)
from vyrtuous.models.category import Category, CategoryObject
from vyrtuous.models.duration import Duration, DurationObject, DurationWrapper
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.utils.channels import video_channel_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import cap_service
from vyrtuous.utils.tracking import stream_service
from vyrtuous.utils.users import moderator_service


class HiddenAdministratorTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Administrator"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx) -> bool:
        if ctx.guild is None:
            raise commands.CheckFailure(
                "This command must be executed inside a server."
            )
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @commands.command(name="aroles", help="Administrator roles.")
    @skip_text_command_help_discovery()
    async def list_administrator_roles_text_command(
        self,
        ctx: commands.Context,
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None:
            guild_snowflake = ctx.guild.id
        else:
            guild_snowflake = guild.id
        if ctx.guild.id != guild_snowflake:
            await moderator_service.check_minimum_role(
                guild_snowflake=guild_snowflake,
                member_snowflake=ctx.author.id,
                lowest_role="Developer",
            )
        pages = await list_administrator_roles.build_pages(
            guild_snowflake=guild_snowflake,
        )
        return await tick.end(success=pages)

    @commands.command(name="cap", help="Cap alias duration for mods.")
    @skip_text_command_help_discovery()
    async def cap_text_command(
        self,
        ctx: commands.Context,
        category: CategoryObject = commands.parameter(
            converter=Category, description="One of: `mute`, `ban`, `tmute`"
        ),
        channel: discord.abc.GuildChannel | None = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
        limit: DurationWrapper = commands.parameter(
            converter=Duration, default=None, description="m/h/d."
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if channel is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = ctx.channel.id
            guild_snowflake = ctx.guild.id
        elif isinstance(
            channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.id
            guild_snowflake = channel.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        if limit is None:
            duration = DurationObject(number=24, prefix="", sign=1, unit="h")
        else:
            duration = limit.duration
        msg = await cap_service.toggle_cap(
            category=category.category,
            channel_snowflake=channel_snowflake,
            duration=duration,
            guild_snowflake=guild_snowflake,
        )
        return await tick.end(success=msg)

    @commands.command(name="caps", help="List caps.")
    @skip_text_command_help_discovery()
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention or server ID.",
        ),
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None:
            guild_snowflake = ctx.guild.id
        else:
            guild_snowflake = guild.id
        obj = target or ctx.channel
        pages = await list_caps.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="debug", help="Shows the last `n` number of logging.")
    @skip_text_command_help_discovery()
    async def debug_text_command(
        self,
        ctx,
        *,
        lines: int = commands.parameter(
            default=3, description="Specify the number of lines"
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if lines <= 0:
            return await tick.end(warning="Lines must be greater than 0")
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

    @commands.command(name="pc", help="View permissions.")
    @skip_text_command_help_discovery()
    async def list_permissions_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: channel ID/mention or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            obj = ctx.channel
        else:
            obj = target
        pages = await list_permissions.build_pages(
            obj=obj,
        )
        return await tick.end(success=pages)

    @commands.command(name="roleid", help="Get role by name.")
    @skip_text_command_help_discovery()
    async def get_role_id_text_command(
        self,
        ctx: commands.Context,
        role_name: str,
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None:
            role = discord.utils.get(ctx.guild.roles, name=role_name)
        else:
            role = discord.utils.get(guild.roles, name=role_name)
        if role:
            return await tick.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await tick.end(
                warning=f"No role named `{role_name}` found in this server."
            )

    @commands.command(name="automutes", help="List automute channels.")
    @skip_text_command_help_discovery()
    async def list_automute_channels_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            discord.Guild, discord.abc.GuildChannel, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: channel ID/mention, or server ID.",
        ),
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None:
            guild_snowflake = ctx.guild.id
        else:
            guild_snowflake = guild.id
        obj = target or ctx.guild
        pages = await list_automute_channels.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @commands.command(name="stream", help="Setup streaming.")
    @skip_text_command_help_discovery()
    async def modify_streaming_text_command(
        self,
        ctx: commands.Context,
        target_channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.TextChannelConverter,
            description="Tag a channel or include its ID.",
        ),
        source: (
            str | discord.abc.GuildChannel | discord.Guild | None
        ) = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify `all` a channel ID/mention or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if source is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            source_obj = ctx.guild
        elif isinstance(source, str):
            await moderator_service.check_minimum_role(
                member_snowflake=ctx.author.id,
                lowest_role="Developer",
            )
            source_obj = None
        else:
            source_obj = source
        if isinstance(
            target_channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            target_channel_obj = target_channel
        else:
            return await tick.end(
                warning="This command must target a valid target channel."
            )
        pages = await stream_service.toggle_stream(
            target_channel=target_channel_obj,
            source=source_obj,
        )
        return await tick.end(success=pages)

    @commands.command(name="streams", help="List streaming routes.")
    @skip_text_command_help_discovery()
    async def list_streaming_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            discord.Guild, discord.abc.GuildChannel, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a channel ID/mention or server a ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            obj = ctx.guild
        else:
            obj = target
        pages = await list_streams.build_pages(obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="v", help="Start/stop video-only channel.")
    @skip_text_command_help_discovery()
    async def toggle_video_channel_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include the ID",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        obj = channel or ctx.channel
        msg = await video_channel_service.toggle_video_channel(
            channel_snowflake=obj.id, guild_snowflake=obj.guild.id
        )
        return await tick.end(success=msg)

    @commands.command(
        name="vs",
        help="List video channels.",
    )
    @skip_text_command_help_discovery()
    async def list_video_channels_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a channel ID/mention or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            obj = ctx.guild
        else:
            obj = target
        pages = await list_video_channels.build_pages(obj=obj)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenAdministratorTextCommands(bot=bot))
