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
from vyrtuous.models.category import AppCategory, CategoryObject
from vyrtuous.models.duration import AppDuration, DurationObject, DurationWrapper
from vyrtuous.models.metadata import metadata
from vyrtuous.models.scope import AppScope, ScopeObject
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    clear_service,
    server_mute_service,
    voice_mute_service,
)
from vyrtuous.permissions import permission_service
from vyrtuous.view.cancel_confirm_view import VerifyView
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

    @metadata(permission="command.moderation.ban")
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
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.moderation.ban"],
        )
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

    @metadata(permission="command.moderation.blacklist")
    @app_commands.command(name="blacklist", description="Blacklist a member.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
        channel="Specify a channel ID/mention.",
    )
    async def toggle_blacklist_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            if interaction.guild is None:
                return await tick.end(warning="This command must used in a server.")
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target valid member.")
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            guild_snowflake=guild_snowflake,
            requested=["command.moderation.blacklist"],
        )
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            member_snowflake=member_snowflake,
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
        )
        msg = await ban_service.toggle_blacklist(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        return await tick.end(success=msg)

    @metadata(permission="command.clear")
    @app_commands.command(name="clear", description="Reset records.")
    @app_commands.describe(
        category="Specify one of: `alias`, `all`, `automute`, `ban`, `flag`, `tmute`, `stream` or `vmute`.",
        scope="Specify one of: `auto`, `server`, `user`",
        target="Specify `all`, a channel ID/mention, member ID/mention, or server ID.",
        guild="Specify a server ID.",
    )
    async def clear_channel_access_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject, AppTarget],
        category: app_commands.Transform[CategoryObject, AppCategory],
        scope: app_commands.Transform[ScopeObject, AppScope],
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
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
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = interaction.channel.id
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.clear"],
        )
        view = VerifyView(
            author_snowflake=interaction.user.id,
            category=str(category),
            obj=target,
        )
        embed = view.build_embed()
        await tick.end(success=embed, view=view)
        await view.wait()
        tick = Tick(bot=self.__bot, interaction=interaction)
        message = await interaction.original_response()
        message_snowflake = message.id
        msg = await clear_service.clear(
            author_snowflake=interaction.user.id,
            category=category.category,
            guild_snowflake=guild_snowflake,
            message_snowflake=message_snowflake,
            message_channel_snowflake=message.channel.id,
            obj=target,
            target=scope.scope,
            view=view,
        )
        return await tick.end(success=msg)

    @metadata(permission="command.moderation.duration")
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
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.moderation.duration"],
        )
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

    @metadata(permission="command.moderation.flag")
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
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.moderation.flag"],
        )
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

    @metadata(permission="command.moderation.voice-mute")
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
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.moderation.voice-mute"],
        )
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

    @metadata(permission="command.moderation.reason")
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
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.moderation.reason"],
        )
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

    @metadata(permission="command.moderation.voice-mute.channel_mute")
    @app_commands.command(name="rmute", description="Room mute (except yourself).")
    @app_commands.describe(
        channel="Specify a channel ID/mention.",
        duration="Specify a duration in m/h/d.",
        reason="Specify a reason.",
    )
    async def channel_mute_app_command(
        self,
        interaction: discord.Interaction,
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
        duration: app_commands.Transform[DurationWrapper | None, AppDuration] = None,
        reason: str | None = "No reason provided.",
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
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
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = interaction.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        if duration is None:
            duration_obj = DurationObject(number=1, prefix="", sign=1, unit="h")
        else:
            duration_obj = duration.duration
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.moderation.voice-mute.channel_mute"],
        )
        pages = await voice_mute_service.channel_mute(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            duration=duration_obj,
            excluded=[interaction.user.id],
            guild_snowflake=guild_snowflake,
            reason=reason or "No reason provided.",
        )
        return await tick.end(success=pages)

    @metadata(permission="command.moderation.voice-mute.server")
    @app_commands.command(name="smute", description="Server mute/server unmute.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
        reason="Specify a reason.",
        guild="Specify a server ID.",
    )
    async def toggle_server_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        reason: str | None = "No reason provided.",
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
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
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = interaction.channel.id
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.moderation.voice-mute.server"],
        )
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            author_snowflake=interaction.user.id,
            member_snowflake=member_snowflake,
        )
        msg = await server_mute_service.toggle_server_mute(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            reason=reason or "No reason provided.",
        )
        return await tick.end(success=msg)

    @metadata(permission="command.moderation.text-mute")
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
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.moderation.text-mute"],
        )
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
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

    @metadata(permission="command.moderation.unvoice-mute.channel_unmute")
    @app_commands.command(name="xrmute", description="Unmute all.")
    @app_commands.describe(
        channel="Specify a channel ID/mention.",
    )
    async def channel_unmute_app_command(
        self,
        interaction: discord.Interaction,
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must used in a server channel."
                )
            channel_snowflake = interaction.channel.id
            if interaction.guild is None:
                return await tick.end(warning="This command must used in a server.")
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
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.moderation.unvoice-mute.channel_unmute"],
        )
        pages = await voice_mute_service.channel_unmute(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            target="click",
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(ModerationAppCommands(bot=bot))
