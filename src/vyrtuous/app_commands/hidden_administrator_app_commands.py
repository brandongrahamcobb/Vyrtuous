"""!/bin/python3

hidden_admin_app_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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
from discord import app_commands
from discord.ext import commands

from vyrtuous.app_commands.help_app_command import skip_app_command_help_discovery
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import PATH_LOG
from vyrtuous.listing import (
    list_administrator_roles,
    list_automute_channels,
    list_caps,
    list_permissions,
    list_streams,
    list_video_channels,
)
from vyrtuous.models.category import AppCategory, CategoryObject
from vyrtuous.models.duration import AppDuration, DurationObject, DurationWrapper
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.utils.channels import video_channel_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import cap_service
from vyrtuous.utils.tracking import stream_service
from vyrtuous.utils.users import moderator_service


class HiddenAdministratorAppCommands(commands.Cog):
    PERMISSION_LEVEL = "Administrator"

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def interaction_check(self, interaction: discord.Interaction):
        if interaction.guild is None:
            raise commands.CheckFailure(
                "This command must be executed inside a server."
            )
        if interaction.channel is None:
            raise commands.CheckFailure(
                "This command must be executed in a valid channel."
            )
        await moderator_service.check_minimum_role(
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @app_commands.command(name="aroles", description="Administrator roles.")
    @app_commands.describe(guild="Specify a server ID.")
    @skip_app_command_help_discovery()
    async def list_administrator_roles_app_command(
        self,
        interaction: discord.Interaction,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None:
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if interaction.guild.id != guild_snowflake:
            await moderator_service.check_minimum_role(
                member_snowflake=interaction.user.id,
                lowest_role="Developer",
            )
        pages = await list_administrator_roles.build_pages(
            guild_snowflake=guild_snowflake
        )
        return await tick.end(success=pages)

    @app_commands.command(name="cap", description="Cap alias duration for mods.")
    @app_commands.describe(
        channel="Specify a channel ID/mention.",
        category="Specify one of: `ban`, `tmute`, or `vmute`.",
        limit="Specify limit in m/d/h.",
    )
    @skip_app_command_help_discovery()
    async def cap_app_command(
        self,
        interaction: discord.Interaction,
        category: app_commands.Transform[CategoryObject, AppCategory],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
        limit: app_commands.Transform[DurationWrapper | None, AppDuration] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
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

    @app_commands.command(name="caps", description="List caps.")
    @app_commands.describe(
        target="Specify one of: channel ID/mention or server ID.",
        guild="Specify a server ID.",
    )
    @skip_app_command_help_discovery()
    async def list_caps_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            obj = interaction.guild
        else:
            obj = target.target
        pages = await list_caps.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @app_commands.command(
        name="debug", description="Shows the last `n` number of logging."
    )
    @app_commands.describe(lines="Specify the number of lines.")
    @skip_app_command_help_discovery()
    async def debug_app_command(
        self, interaction: discord.Interaction, lines: int | None = 3
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if lines is None:
            return await tick.end(warning="Lines must not be `None`.")
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

    @app_commands.command(name="pc", description="View permissions.")
    @app_commands.describe(
        target="Specify one of: channel ID/mention or server ID.",
    )
    @skip_app_command_help_discovery()
    async def list_permissions_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            obj = interaction.guild
        else:
            obj = target.target
        pages = await list_permissions.build_pages(
            obj=obj,
        )
        return await tick.end(success=pages)

    @app_commands.command(name="roleid", description="Get role by name.")
    @app_commands.describe(
        role_name="Specify a role name.", guild="Specify a server ID."
    )
    @skip_app_command_help_discovery()
    async def get_role_id_app_command(
        self,
        interaction: discord.Interaction,
        role_name: str,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_obj = interaction.guild
        else:
            if isinstance(guild.target, discord.Guild):
                guild_obj = guild.target
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        role = discord.utils.get(guild_obj.roles, name=role_name)
        if role:
            return await tick.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await tick.end(
                warning=f"No role named `{role_name}` found in this server."
            )

    @app_commands.command(name="automutes", description="List automute channels.")
    @app_commands.describe(
        target="Specify one of: channel ID/mention or server ID.",
        guild="Specify a server ID.",
    )
    @skip_app_command_help_discovery()
    async def list_automute_channels_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if interaction.channel is None:
            return await tick.end(warning="This command must target a valid channel.")
        if guild is None:
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if interaction.guild.id != guild_snowflake:
            await moderator_service.check_minimum_role(
                member_snowflake=interaction.user.id,
                lowest_role="Developer",
            )
        if target is None:
            obj = interaction.guild
        else:
            obj = target.target
        pages = await list_automute_channels.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @app_commands.command(name="stream", description="Setup streaming.")
    @app_commands.describe(
        target_channel="Specify a channel ID/mention where logs will be streamed.",
        source="Specify a channel ID/mention or server ID where logs will come from.",
    )
    @skip_app_command_help_discovery()
    async def modify_streaming_app_command(
        self,
        interaction: discord.Interaction,
        target_channel: app_commands.Transform[TargetObject, AppTarget],
        source: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if source is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            source_obj = interaction.guild
        elif isinstance(source.target, str):
            await moderator_service.check_minimum_role(
                member_snowflake=interaction.user.id,
                lowest_role="Developer",
            )
            source_obj = None
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
        pages = await stream_service.toggle_stream(
            target_channel=target_channel_obj,
            source=source_obj,
        )
        return await tick.end(success=pages)

    @app_commands.command(name="streams", description="List streaming routes.")
    @app_commands.describe(
        target="Specify one of: a channel ID/mention or server ID",
    )
    @skip_app_command_help_discovery()
    async def list_streaming_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            obj = interaction.guild
        else:
            obj = target.target
        pages = await list_streams.build_pages(obj=obj)
        return await tick.end(success=pages)

    @app_commands.command(name="v", description="Start/stop video-only channel.")
    @app_commands.describe(channel="Specify a channel ID/mention.")
    @skip_app_command_help_discovery()
    async def toggle_video_channel_app_command(
        self,
        interaction: discord.Interaction,
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
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
        msg = await video_channel_service.toggle_video_channel(
            channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake
        )
        return await tick.end(success=msg)

    @app_commands.command(
        name="vs",
        description="List video channels.",
    )
    @app_commands.describe(target="Specify one of: a channel ID/mention or server ID.")
    @skip_app_command_help_discovery()
    async def list_video_channels_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            obj = interaction.guild
        else:
            obj = target.target
        pages = await list_video_channels.build_pages(obj=obj)
        return await tick.end(success=pages)

    async def cog_app_command_error(self, interaction, error):
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenAdministratorAppCommands(bot=bot))
