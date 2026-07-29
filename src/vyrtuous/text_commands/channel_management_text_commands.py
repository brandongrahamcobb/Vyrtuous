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

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.models.category import Category, CategoryObject
from vyrtuous.models.duration import Duration, DurationObject, DurationWrapper
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.utils.channels import automute_channel_service, video_channel_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import cap_service
from vyrtuous.utils.permissions import permission_service
from vyrtuous.utils.tracking import stream_service


class ChannelManagementTextCommands(commands.Cog):

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    @commands.command(name="automute", help="Start/stop automute")
    async def toggle_automute_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel | None = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a channel ID/mention.",
        ),
        duration: DurationWrapper | None = commands.parameter(
            converter=Duration,
            default=None,
            description="Specify duration as m/h/d.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if channel is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = ctx.channel.id
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = ctx.guild.id
        elif isinstance(
            channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.id
            guild_snowflake = channel.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            guild_snowflake=guild_snowflake,
            channel_snowflake=channel_snowflake,
            requested=["command.channel.automute"],
        )
        if duration is None:
            duration_obj = DurationObject(number=2, prefix="+", sign=1, unit="h")
        else:
            duration_obj = duration.duration
        pages = await automute_channel_service.toggle_automute(
            author_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            duration=duration_obj,
        )
        return await tick.end(success=pages)

    @commands.command(name="cap", help="Cap duration.")
    async def cap_text_command(
        self,
        ctx: commands.Context,
        category: CategoryObject = commands.parameter(
            converter=Category, description="One of: `ban`, `tmute` or `vmute`."
        ),
        channel: discord.abc.GuildChannel | None = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Specify a channel ID/mention.",
        ),
        limit: DurationWrapper | None = commands.parameter(
            converter=Duration, default=None, description="Specify limit in m/h/d."
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if channel is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
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
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            guild_snowflake=guild_snowflake,
            channel_snowflake=channel_snowflake,
            requested=["command.channel.cap"],
        )
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

    @commands.command(name="stream", help="Setup streaming.")
    async def modify_streaming_text_command(
        self,
        ctx: commands.Context,
        target_channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.TextChannelConverter,
            description="Specify a channel ID/mention or server ID where logs will be streamed.",
        ),
        source: discord.abc.GuildChannel | discord.Guild | None = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a channel ID/mention or server ID where logs will come from.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if source is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            source_obj = ctx.guild
        elif isinstance(
            source,
            (
                discord.Guild,
                discord.VoiceChannel,
                discord.TextChannel,
                discord.StageChannel,
            ),
        ):
            source_obj = source
        else:
            return await tick.end(
                warning="This command must target a valid channel or server."
            )
        if isinstance(
            target_channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            target_channel_obj = target_channel
        else:
            return await tick.end(
                warning="This command must target a valid target channel."
            )
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            guild_snowflake=target_channel_obj.guild.id,
            channel_snowflake=target_channel_obj.id,
            requested=["command.channel.stream"],
        )
        pages = await stream_service.toggle_stream(
            target_channel=target_channel_obj,
            source=source_obj,
        )
        return await tick.end(success=pages)

    @commands.command(name="v", help="Start/stop video-only channel.")
    async def toggle_video_channel_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel | None = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Specify a channel ID/mention.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if channel is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = ctx.channel.id
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = ctx.guild.id
        elif isinstance(
            channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.id
            guild_snowflake = channel.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            guild_snowflake=guild_snowflake,
            channel_snowflake=channel_snowflake,
            requested=["command.channel.video"],
        )
        msg = await video_channel_service.toggle_video_channel(
            channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake
        )
        return await tick.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(ChannelManagementTextCommands(bot=bot))
