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
from vyrtuous.models.duration import Duration, DurationObject, DurationWrapper
from vyrtuous.models.metadata import metadata
from vyrtuous.models.target import Target, TargetObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    server_mute_service,
    voice_mute_service,
)


class ModerationTextCommands(commands.Cog):

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    @commands.command(name="blacklist", help="Blacklist a member.")
    @metadata(permission="command.moderation.blacklist")
    async def toggle_blacklist_text_command(
        self,
        ctx: commands.Context,
        member: TargetObject = commands.parameter(
            converter=Target,
            description="Specify a member ID/mention.",
        ),
        channel: TargetObject | None = commands.parameter(
            converter=Target,
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
                    warning="This command must be used in a server channel."
                )
            if ctx.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            channel_snowflake = ctx.channel.id
            guild_snowflake = ctx.guild.id
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
            member_snowflake=ctx.author.id,
            guild_snowflake=guild_snowflake,
            requested=["command.moderation.blacklist"],
        )
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            member_snowflake=member_snowflake,
            author_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
        )
        msg = await ban_service.toggle_blacklist(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        return await tick.end(success=msg)

    @commands.command(name="rmute", help="Room mute (except yourself).")
    @metadata(permission="command.moderation.voice-mute.channel_mute")
    async def channel_mute_text_command(
        self,
        ctx: commands.Context,
        channel: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a channel ID/mention.",
        ),
        duration: DurationWrapper | None = commands.parameter(
            converter=Duration,
            default=None,
            description="Specify a duration m/h/d.",
        ),
        reason: str = commands.parameter(
            default="No reason provided.", description="Specify a reason."
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if channel is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = ctx.channel.id
            if ctx.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = ctx.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        if duration is None:
            duration_obj = DurationObject(number=1, prefix="", sign=1, unit="h")
        else:
            duration_obj = duration.duration
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.moderation.voice-mute.channel_mute"],
        )
        pages = await voice_mute_service.channel_mute(
            author_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            duration=duration_obj,
            excluded=[ctx.author.id],
            guild_snowflake=guild_snowflake,
            reason=reason,
        )
        return await tick.end(success=pages)

    @commands.command(name="smute", help="Server mute/server unmute.")
    @metadata(permission="command.moderation.voice-mute.server")
    async def toggle_server_mute_text_command(
        self,
        ctx: commands.Context,
        member: TargetObject = commands.parameter(
            converter=Target,
            description="Specify a member ID/mention.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
        reason: str = commands.parameter(
            default="No reason provided",
            description="Optional reason (required for 7 days or more)",
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
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.moderation.voice-mute.server"],
        )
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            author_snowflake=ctx.author.id,
            member_snowflake=member_snowflake,
        )
        msg = await server_mute_service.toggle_server_mute(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            reason=reason,
        )
        return await tick.end(success=msg)

    @commands.command(name="xrmute", help="Unmute all.")
    @metadata(permission="command.moderation.unvoice-mute.channel_unmute")
    async def channel_unmute_text_command(
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
        if channel is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = ctx.channel.id
            if ctx.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = ctx.guild.id
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
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.moderation.unvoice-mute.channel_unmute"],
        )
        pages = await voice_mute_service.channel_unmute(
            author_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            target="click",
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(ModerationTextCommands(bot=bot))
