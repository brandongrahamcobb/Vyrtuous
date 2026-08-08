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
from vyrtuous.db.ban import Ban
from vyrtuous.db.flag import Flag
from vyrtuous.db.text_mute import TextMute
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.modal.duration_modal import DurationModal
from vyrtuous.modal.reason_modal import ReasonModal
from vyrtuous.models.duration import AppDuration, DurationObject, DurationWrapper
from vyrtuous.models.metadata import metadata
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import CheckFailure
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    server_mute_service,
    voice_mute_service,
)
from vyrtuous.view.clear_view import ClearView
from vyrtuous.view.infraction_view import InfractionView
from vyrtuous.view.modify_view import ModifyView


class ModerationAppCommands(commands.Cog):

    def __init__(
        self,
        bot: DiscordBot,
    ):
        self.__bot = bot

    async def cog_app_command_error(self, interaction, error):
        self.__bot.logger.info(str(error))
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error), ephemeral=True)

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
                return await tick.end(
                    warning="This command must be used in a server.", ephemeral=True
                )
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel.",
                    ephemeral=True,
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server.", ephemeral=True
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
                    warning="This command must target a valid channel.", ephemeral=True
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(
                warning=f"This command must target a valid member.", ephemeral=True
            )
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.moderation.ban"],
        )
        ctx = SnowflakeContext(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            interaction=interaction,
            model=Ban,
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
                    warning="This command must be used in a server channel.",
                    ephemeral=True,
                )
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must used in a server.", ephemeral=True
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(
                warning="This command must target a valid channel.", ephemeral=True
            )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(
                warning=f"This command must target a valid member.", ephemeral=True
            )
        await tick.defer()
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
        target="Specify `all`, a channel ID/mention, member ID/mention, or server ID.",
    )
    async def clear_channel_access_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject, AppTarget],
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
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
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.clear"],
        )
        if not isinstance(
            target.target,
            (str, int, discord.Guild, discord.Member, discord.abc.GuildChannel),
        ):
            return await tick.end(
                warning="Invalid target. Must be a `channel`, `member`, `server` or `all`."
            )
        if isinstance(target.target, int):
            await permission_service.has_equal_or_lower_role(
                permission_state=permission_state,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                author_snowflake=interaction.user.id,
                member_snowflake=target.target,
            )
        elif isinstance(target.target, discord.Member):
            await permission_service.has_equal_or_lower_role(
                permission_state=permission_state,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                author_snowflake=interaction.user.id,
                member_snowflake=target.target.id,
            )
        ctx = SnowflakeContext(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=None,
        )
        view = ClearView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            interaction=interaction,
            obj=target.target,
            tick=tick,
        )
        try:
            await view.setup()
        except CheckFailure as e:
            return await tick.end(warning=str(e))
        await interaction.response.send_message(
            content="Specify what to clear.", view=view, ephemeral=True
        )

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
                return await tick.end(
                    warning="This command must be used in a server.", ephemeral=True
                )
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel.",
                    ephemeral=True,
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server.", ephemeral=True
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
                    warning="This command must target a valid channel.", ephemeral=True
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(
                warning=f"This command must target a valid member.", ephemeral=True
            )
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.moderation.duration"],
        )
        ctx = SnowflakeContext(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = ModifyView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            interaction=interaction,
            modal=DurationModal,
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
                return await tick.end(
                    warning="This command must be used in a server.", ephemeral=True
                )
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel.",
                    ephemeral=True,
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server.", ephemeral=True
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
                    warning="This command must target a valid channel.", ephemeral=True
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(
                warning=f"This command must target a valid member.", ephemeral=True
            )
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
        ctx = SnowflakeContext(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            interaction=interaction,
            model=Flag,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the flag", view=view, ephemeral=True
        )

    @metadata(permission="command.moderation.voice-mute")
    @app_commands.command(name="mute", description="Create a voice mute.")
    @app_commands.describe(
        member="Specify a member ID/mention.", channel="Specify a channel ID/mention."
    )
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
                return await tick.end(
                    warning="This command must be used in a server.", ephemeral=True
                )
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel.",
                    ephemeral=True,
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server.", ephemeral=True
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
                    warning="This command must target a valid channel.", ephemeral=True
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(
                warning=f"This command must target a valid member.", ephemeral=True
            )
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
        if await voice_mute_service.is_voice_muted(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            targets=["server"],
        ):
            return
        ctx = SnowflakeContext(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            interaction=interaction,
            model=VoiceMute,
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
                return await tick.end(
                    warning="This command must be used in a server.", ephemeral=True
                )
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel.",
                    ephemeral=True,
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server.", ephemeral=True
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
                    warning="This command must target a valid channel.", ephemeral=True
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(
                warning=f"This command must target a valid member.", ephemeral=True
            )
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.moderation.reason"],
        )
        ctx = SnowflakeContext(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = ModifyView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            interaction=interaction,
            modal=ReasonModal,
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
                    warning="This command must be used in a server channel.",
                    ephemeral=True,
                )
            channel_snowflake = interaction.channel.id
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.target.id
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server.", ephemeral=True
                )
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(
                warning="This command must target a valid channel.", ephemeral=True
            )
        if duration is None:
            duration_obj = DurationObject(number=1, prefix="", sign=1, unit="h")
        else:
            duration_obj = duration.duration
        await tick.defer()
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
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(
                warning=f"This command must target a valid member.", ephemeral=True
            )
        await tick.defer()
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
                return await tick.end(
                    warning="This command must be used in a server.", ephemeral=True
                )
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel.",
                    ephemeral=True,
                )
            channel_snowflake = interaction.channel.id
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(channel.target, discord.abc.GuildChannel):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server.", ephemeral=True
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
                    warning="This command must target a valid channel.", ephemeral=True
                )
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(
                warning=f"This command must target a valid member.", ephemeral=True
            )
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.moderation.text-mute"],
        )
        ctx = SnowflakeContext(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            interaction=interaction,
            model=TextMute,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the text-mute.", view=view, ephemeral=True
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
                    warning="This command must used in a server channel.",
                    ephemeral=True,
                )
            channel_snowflake = interaction.channel.id
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must used in a server.", ephemeral=True
                )
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(
                warning="This command must target a valid channel.", ephemeral=True
            )
        await tick.defer()
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
