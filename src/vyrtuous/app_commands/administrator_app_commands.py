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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.listing import list_overwrites, list_server_mutes
from vyrtuous.models.category import AppCategory, CategoryObject
from vyrtuous.models.duration import AppDuration, DurationObject, DurationWrapper
from vyrtuous.models.scope import AppScope, ScopeObject
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    clear_service,
    server_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.users import coordinator_service, moderator_service
from vyrtuous.view.cancel_confirm_view import VerifyView


class AdministratorAppCommands(commands.Cog):
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

    @app_commands.command(name="clear", description="Reset records.")
    @app_commands.describe(
        category="Specify one of: `admin`, `alias`, `all`, `automute`, `ban`, `coord`, `flag`, `mod`, `tmute`, `stream` or `vmute`.",
        scope="Specify one of: `auto`, `server`, `user`",
        target="Specify one of: channel ID/mention, member ID/mention, or server ID.",
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

    @app_commands.command(name="coord", description="Grant/revoke coords.")
    @app_commands.describe(
        channel="Specify a channel ID/mention",
        member="Specify a member ID/mention",
    )
    async def toggle_coordinator_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
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
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        await moderator_service.has_equal_or_lower_role(
            target_member_snowflake=member_snowflake,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
        )
        message = await interaction.original_response()
        message_snowflake = message.id
        msg = await coordinator_service.toggle_coordinator(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake,
            message_channel_snowflake=message.channel.id,
        )
        return await tick.end(success=msg)

    @app_commands.command(name="ow", description="Overwrite stats.")
    @app_commands.describe(
        target="Specify one of: a channel ID/mention or server ID.",
    )
    async def list_overwrites_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        embed = list_overwrites.build_embed(obj=obj)
        return await tick.end(success=embed)

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
        pages = await voice_mute_service.channel_mute(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            duration=duration_obj,
            excluded=[interaction.user.id],
            guild_snowflake=guild_snowflake,
            reason=reason or "No reason provided.",
        )
        return await tick.end(success=pages)

    @app_commands.command(name="rmv", description="VC move.")
    @app_commands.describe(
        target_channel="Specify a `to` channel ID/mention.",
        source_channel="Specify a `from` channel ID/mention.",
    )
    async def channel_move_all_app_command(
        self,
        interaction: discord.Interaction,
        target_channel: app_commands.Transform[TargetObject, AppTarget],
        source_channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if isinstance(
            target_channel.target,
            (discord.VoiceChannel, discord.StageChannel),
        ):
            target_guild_snowflake = target_channel.target.guild.id
            target_guild_name = target_channel.target.guild.name
            target_channel_snowflake = target_channel.target.id
            target_channel_name = target_channel.target.name
            target_channel_obj = target_channel.target
            target_channel_mention = target_channel.target.mention
        else:
            return await tick.end(
                warning="This command must target a valid target channel."
            )
        if source_channel is None:
            if interaction.channel is None or not isinstance(
                interaction.channel, (discord.VoiceChannel, discord.StageChannel)
            ):
                return await tick.end(
                    warning="This command must target a valid source channel."
                )
            source_guild_snowflake = interaction.channel.guild.id
            source_guild_name = interaction.channel.guild.name
            source_channel_snowflake = interaction.channel.id
            source_channel_name = interaction.channel.name
            source_channel_members = interaction.channel.members
            source_channel_mention = interaction.channel.mention
        elif isinstance(
            source_channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            source_guild_snowflake = source_channel.target.guild.id
            source_guild_name = source_channel.target.guild.name
            source_channel_snowflake = source_channel.target.id
            source_channel_name = source_channel.target.name
            source_channel_members = source_channel.target.members
            source_channel_mention = source_channel.target.mention
        else:
            return await tick.end(
                warning="This command must target a valid source channel."
            )
        failed, moved = [], []
        for member in source_channel_members:
            try:
                await member.move_to(target_channel_obj)
                moved.append(member)
            except discord.Forbidden as e:
                failed.append(member)
                self.__bot.logger.warning(
                    f"Unable to move member "
                    f"{member.display_name} ({member.id}) from channel "
                    f"{source_channel_name} ({source_channel_snowflake}) in guild "
                    f"{source_guild_name} ({source_guild_snowflake}) to channel "
                    f"{target_channel_name} ({target_channel_snowflake}) in guild "
                    f"{target_guild_name} ({target_guild_snowflake}) to channel "
                    f"{str(e).capitalize()}"
                )
        embed = discord.Embed(
            title=f"{emojis.get_random_emoji()} "
            f"Moved {source_channel_mention} to "
            f"{target_channel_mention}",
            color=discord.Color.green(),
        )
        if moved:
            embed.add_field(
                name=f"Successfully Moved (`{len(moved)}`)",
                value="\n".join(member.mention for member in moved),
                inline=False,
            )
        else:
            embed.add_field(name="Successfully Moved", value="None", inline=False)
        if failed:
            embed.add_field(
                name=f"Failed to Move ({len(failed)})",
                value="\n".join(member.mention for member in failed),
                inline=False,
            )
        embed.set_footer(
            text=f"Moved from {source_channel_name} to {target_channel_name}"
        )
        return await tick.end(success=embed)

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
        if isinstance(member.target, int):
            member_snowflake = member.target
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        msg = await server_mute_service.toggle_server_mute(
            author_snowflake=interaction.user.id,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            reason=reason or "No reason provided.",
        )
        return await tick.end(success=msg)

    @app_commands.command(name="smutes", description="List mutes.")
    @app_commands.describe(
        target="Specify one of: a member ID/mention or a server ID.",
        guild="Specify a server ID.",
    )
    async def list_server_mutes_app_command(
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
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        pages = await list_server_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

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
        if channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = interaction.channel.id
            if interaction.guild is None:
                return await tick.end(warning="This command must target a valid guild.")
            guild_snowflake = interaction.guild.id
        elif isinstance(
            channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.target.id
            guild_snowflake = channel.target.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        pages = await voice_mute_service.channel_unmute(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            target="user",
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(AdministratorAppCommands(bot=bot))
