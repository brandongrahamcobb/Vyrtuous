"""!/bin/python3
moderation_app_commands.py A discord.py cog containing moderation commands for the Vyrtuous bot.

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
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import voice_mute_service
from vyrtuous.utils.permissions import permission_service
from vyrtuous.view.infraction_view import InfractionView
from vyrtuous.view.modify_infraction_view import ModifyInfractionView
from vyrtuous.view.view_context import ViewContext


class ModerationAppCommands(commands.Cog):

    def __init__(
        self,
        bot: DiscordBot,
    ):
        self.__bot = bot

    async def cog_app_command_error(self, interaction, error):
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))

    @app_commands.command(name="ban", description="Create a ban.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_ban_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if channel.target.guild.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                channel_snowflake = channel.target.id
                guild_snowflake = channel.target.guild.id
            else:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
        if interaction.channel is None:
            return await tick.end(warning="This command must used in a server channel.")
        else:
            channel_snowflake = interaction.channel.id
        ctx = ViewContext(
            interaction=interaction,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        ctx.category = "ban"
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the ban", view=view, ephemeral=True
        )

    @app_commands.command(name="duration", description="Modify a duration.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def change_moderation_duration_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if channel.target.guild.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                channel_snowflake = channel.target.id
                guild_snowflake = channel.target.guild.id
            else:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = interaction.channel.id
        ctx = ViewContext(
            interaction=interaction,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = ModifyInfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            modal="duration",
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a category", view=view, ephemeral=True
        )

    @app_commands.command(name="flag", description="Create a flag.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_flag_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if channel.target.guild.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                channel_snowflake = channel.target.id
                guild_snowflake = channel.target.guild.id
            else:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
        if interaction.channel is None:
            return await tick.end(warning="This command must used in a server channel.")
        else:
            channel_snowflake = interaction.channel.id

        ctx = ViewContext(
            interaction=interaction,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        ctx.category = "flag"
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the flag", view=view, ephemeral=True
        )

    @app_commands.command(name="mute", description="Create a voice mute.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_voice_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message | None:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if channel.target.guild.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                channel_snowflake = channel.target.id
                guild_snowflake = channel.target.guild.id
            else:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = interaction.channel.id
        if await voice_mute_service.is_voice_muted(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            targets=["server"],
        ):
            return
        ctx = ViewContext(
            interaction=interaction,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        ctx.category = "vmute"
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the voice-mute", view=view, ephemeral=True
        )

    @app_commands.command(name="reason", description="Modify a reason.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def change_moderation_reason_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if channel.target.guild.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                channel_snowflake = channel.target.id
                guild_snowflake = channel.target.guild.id
            else:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = interaction.channel.id
        ctx = ViewContext(
            interaction=interaction,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = ModifyInfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            modal="reason",
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a category", view=view, ephemeral=True
        )

    @app_commands.command(name="tmute", description="Create a text-mute.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_app_text_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if channel.target.guild.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                channel_snowflake = channel.target.id
                guild_snowflake = channel.target.guild.id
            else:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = interaction.channel.id
        ctx = ViewContext(
            interaction=interaction,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        ctx.category = "tmute"
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the text-mute", view=view, ephemeral=True
        )
        await interaction.response.defer()


async def setup(bot: DiscordBot):
    await bot.add_cog(ModerationAppCommands(bot=bot))
