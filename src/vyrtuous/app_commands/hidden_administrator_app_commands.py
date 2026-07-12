"""!/bin/python3

admin_app_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.app_commands.help_app_command import skip_app_command_help_discovery
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
from vyrtuous.models import multi_converter
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
            raise commands.CheckFailure("This command must be used inside a server.")
        if interaction.channel is None:
            raise commands.CheckFailure("This command must be used in a valid channel.")
        await moderator_service.check_minimum_role(
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @app_commands.command(name="aroles", description="Administrator roles.")
    @app_commands.describe(target="Specify a server ID")
    @skip_app_command_help_discovery()
    async def list_administrator_roles_app_command(
        self, interaction: discord.Interaction, target: str | None
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.guild
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_administrator_roles.build_pages(
            obj=obj, is_at_home=is_at_home
        )
        return await tick.end(success=pages)

    @app_commands.command(name="cap", description="Cap alias duration for mods.")
    @app_commands.describe(
        channel="Specify a channel ID/mention.",
        category="Specify one of: `ban`, `tmute`, or `vmute`.",
        limit="Specify limit in m/d/h",
    )
    @skip_app_command_help_discovery()
    async def cap_app_command(
        self,
        interaction: discord.Interaction,
        channel: discord.abc.GuildChannel,
        category: str,
        limit: str | None = "24h",
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        msg = await cap_service.toggle_cap(
            category=str(category), channel=channel, duration_str=limit or "24h"
        )
        return await tick.end(success=msg)

    @app_commands.command(name="caps", description="List caps.")
    @app_commands.describe(target="Specify one of: channel ID/mention or server ID.")
    @skip_app_command_help_discovery()
    async def list_caps_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_caps.build_pages(obj=obj, is_at_home=is_at_home)
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
    @app_commands.describe(target="Specify one of: channel ID/mention or server ID.")
    @skip_app_command_help_discovery()
    async def list_permissions_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_permissions.build_pages(
            guild_snowflake=interaction.guild.id,
            obj=obj,
            is_at_home=is_at_home,
        )
        return await tick.end(success=pages)

    @app_commands.command(name="roleid", description="Get role by name.")
    @app_commands.describe(role_name="Specify a role name.")
    @skip_app_command_help_discovery()
    async def get_role_id_app_command(
        self, interaction: discord.Interaction, role_name: str
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used within a server.")
        role = discord.utils.get(interaction.guild.roles, name=role_name)
        if role:
            return await tick.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await tick.end(
                warning=f"No role named `{role_name}` found in this server."
            )

    @app_commands.command(name="ams", description="List automute channels.")
    @app_commands.describe(target="Specify one of: channel ID/mention or server ID.")
    @skip_app_command_help_discovery()
    async def list_automute_channels_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_automute_channels.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @app_commands.command(name="stream", description="Setup streaming.")
    @app_commands.describe(
        target_channel="Specify a channel ID/mention where logs will be streamed",
        source_channel="Specify a channel ID/mention where logs will come from",
    )
    @skip_app_command_help_discovery()
    async def modify_streaming_app_command(
        self,
        interaction: discord.Interaction,
        target_channel: str,
        source_channel: str | None = "all",
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        target = multi_converter.transform(
            interaction=interaction, argument=target_channel
        )
        if not isinstance(
            target, (discord.VoiceChannel, discord.StageChannel, discord.TextChannel)
        ):
            return await tick.end(warning="This command must target a valid channel.")
        if source_channel == "all":
            pages = await stream_service.toggle_stream(
                guild_snowflake=interaction.guild.id,
                source_channel_snowflake=None,
                target_channel_snowflake=target.id,
            )
        else:
            source = multi_converter.transform(
                interaction=interaction, argument=source_channel
            )
            if not isinstance(
                source,
                (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
            ):
                return await tick.end(
                    warning="This command must source from a valid channel."
                )
            pages = await stream_service.toggle_stream(
                guild_snowflake=interaction.guild.id,
                source_channel_snowflake=source.id,
                target_channel_snowflake=target.id,
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
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_streams.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @app_commands.command(name="v", description="Start/stop video-only channel.")
    @app_commands.describe(channel="Specify a channel ID/mention.")
    @skip_app_command_help_discovery()
    async def toggle_video_channel_app_command(
        self, interaction: discord.Interaction, channel: discord.abc.GuildChannel | None
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        obj = channel or interaction.channel
        if not isinstance(obj, (discord.VoiceChannel, discord.StageChannel)):
            return await tick.end(
                warning="This command must be executed for a valid channel."
            )
        msg = await video_channel_service.toggle_video_channel(
            channel_snowflake=obj.id, guild_snowflake=obj.guild.id
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
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_video_channels.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenAdministratorAppCommands(bot=bot))
