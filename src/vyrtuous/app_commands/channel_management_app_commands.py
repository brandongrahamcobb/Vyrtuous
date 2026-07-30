"""!/bin/python3

hidden_admin_app_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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
from vyrtuous.models.category import AppCategory, CategoryObject
from vyrtuous.models.duration import AppDuration, DurationObject, DurationWrapper
from vyrtuous.models.metadata import metadata
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.channels import automute_channel_service, video_channel_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import cap_service
from vyrtuous.utils.tracking import stream_service


class ChannelManagementAppCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    @metadata(permission="command.channel.automute")
    @app_commands.command(name="automute", description="Start/stop automute")
    @app_commands.describe(
        channel="Specify a channel ID/mention.",
        duration="Specify duration as m/h/d.",
    )
    async def toggle_automute_app_command(
        self,
        interaction: discord.Interaction,
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
        duration: app_commands.Transform[DurationWrapper | None, AppDuration] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = interaction.channel.id
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            guild_snowflake=guild_snowflake,
            channel_snowflake=channel_snowflake,
            requested=["command.channel.automute"],
        )
        if duration is None:
            duration_obj = DurationObject(number=2, prefix="", sign=1, unit="h")
        else:
            duration_obj = duration.duration
        pages = await automute_channel_service.toggle_automute(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            duration=duration_obj,
        )
        return await tick.end(success=pages)

    @metadata(permission="command.channel.cap")
    @app_commands.command(name="cap", description="Cap duration.")
    @app_commands.describe(
        channel="Specify a channel ID/mention.",
        category="Specify one of: `ban`, `tmute`, or `vmute`.",
        limit="Specify limit in m/d/h.",
    )
    async def cap_app_command(
        self,
        interaction: discord.Interaction,
        category: app_commands.Transform[CategoryObject, AppCategory],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
        limit: app_commands.Transform[DurationWrapper | None, AppDuration] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = interaction.channel.id
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
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
            guild_snowflake=guild_snowflake,
            duration=duration,
        )
        return await tick.end(success=msg)

    @metadata(permission="command.channel.stream")
    @app_commands.command(name="stream", description="Setup streaming.")
    @app_commands.describe(
        target_channel="Specify a channel ID/mention where logs will be streamed.",
        source="Specify a channel ID/mention or server ID where logs will come from.",
    )
    async def modify_streaming_app_command(
        self,
        interaction: discord.Interaction,
        target_channel: app_commands.Transform[TargetObject, AppTarget],
        source: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if source is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            source_obj = interaction.guild
        elif isinstance(
            source.target,
            (
                discord.Guild,
                discord.VoiceChannel,
                discord.TextChannel,
                discord.StageChannel,
            ),
        ):
            source_obj = source.target
        else:
            return await tick.end(
                warning="This command must target a valid channel or server."
            )
        if isinstance(
            target_channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            target_channel_obj = target_channel.target
        else:
            return await tick.end(
                warning="This command must target a valid target channel."
            )
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            guild_snowflake=target_channel_obj.guild.id,
            channel_snowflake=target_channel_obj.id,
            requested=["command.channel.stream"],
        )
        pages = await stream_service.toggle_stream(
            target_channel=target_channel_obj,
            source=source_obj,
        )
        return await tick.end(success=pages)

    @metadata(permission="command.channel.video-channel")
    @app_commands.command(
        name="video-only", description="Start/stop video-only channel."
    )
    @app_commands.describe(channel="Specify a channel ID/mention.")
    async def toggle_video_only_channel_app_command(
        self,
        interaction: discord.Interaction,
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
        duration: app_commands.Transform[DurationWrapper | None, AppDuration] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = interaction.channel.id
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        if duration is None:
            duration_obj = DurationObject(number=2, prefix="", sign=1, unit="h")
        else:
            duration_obj = duration.duration
        if cap_service.exceeds_cap(
            category="vmute",
            channel_snowflake=channel_snowflake,
            duration=duration_obj,
            guild_snowflake=guild_snowflake,
        ):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                guild_snowflake=guild_snowflake,
                channel_snowflake=channel_snowflake,
                requested=["command.moderation.uncapped"],
            )
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            guild_snowflake=guild_snowflake,
            channel_snowflake=channel_snowflake,
            requested=["command.channel.video-channel"],
        )
        msg = await video_channel_service.toggle_video_channel(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            duration=duration_obj,
        )
        return await tick.end(success=msg)

    async def cog_app_command_error(self, interaction, error):
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))


async def setup(bot: DiscordBot):
    await bot.add_cog(ChannelManagementAppCommands(bot=bot))
