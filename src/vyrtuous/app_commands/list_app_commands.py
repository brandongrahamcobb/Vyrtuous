"""!/bin/python3
moderator_app_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.listing import list_bans, list_flags, list_text_mutes, list_voice_mutes
from vyrtuous.models.scope import AppScope, ScopeObject
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import voice_mute_service
from vyrtuous.utils.permissions import permission_service
from vyrtuous.view.infraction_view import InfractionView
from vyrtuous.view.modify_infraction_view import ModifyInfractionView
from vyrtuous.view.view_context import ViewContext


class ModeratorAppCommands(commands.Cog):

    def __init__(
        self,
        bot: DiscordBot,
    ):
        self.__bot = bot

    @app_commands.command(name="duration", description="Modify a duration.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def change_moderation_duration_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
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
        if interaction.channel is None:
            return await tick.end(warning="This command must target a valid channel.")
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

    @app_commands.command(name="reason", description="Modify a reason.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def change_moderation_reason_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
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
        if interaction.channel is None:
            return await tick.end(warning="This command must target a valid channel.")
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

    @app_commands.command(name="ban", description="Create a ban.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_ban_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
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
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
        if interaction.channel is None:
            return await tick.end(warning="This command must target a valid channel.")
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

    @app_commands.command(name="flag", description="Create a flag.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_flag_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
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
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
        if interaction.channel is None:
            return await tick.end(warning="This command must target a valid channel.")
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
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message | None:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
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
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
        if interaction.channel is None:
            return await tick.end(warning="This command must target a valid channel.")
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

    @app_commands.command(name="tmute", description="Create a text-mute.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_app_text_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
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
        if invincible := bot.registry.get(MemberState).invincible.get(
            guild_snowflake, None
        ):
            if member_snowflake in invincible:
                return
        if interaction.channel is None:
            return await tick.end(warning="This command must target a valid channel.")
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

    @app_commands.command(name="bans", description="List bans.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_bans_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.scope.guild"],
                )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning=f"This command must be used in a server channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.guild"],
            )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.channel"],
            )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.member"],
            )
        pages = await list_bans.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    # @app_commands.command(name="coords", description="Lists coords.")
    # @app_commands.describe(
    #     target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    #     guild="Specify a server ID.",
    # )
    # async def list_coordinators_app_command(
    #     self,
    #     interaction: discord.Interaction,
    #     target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    #     guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    # ) -> discord.Message:
    #     tick = Tick(bot=self.__bot, interaction=interaction)
    #     bot: DiscordBot = DiscordBot.get_instance()
    #     permission_state = bot.registry.get(PermissionState)
    #     if guild is None:
    #         if interaction.guild is None:
    #             return await tick.end(
    #                 warning="This command must be used in a valid server."
    #             )
    #         guild_snowflake = interaction.guild.id
    #     else:
    #         if isinstance(guild.target, discord.Guild):
    #             if interaction.guild is None:
    #                 return await tick.end(
    #                     warning="This command must be used in a server."
    #                 )
    #             if guild.target.id != interaction.guild.id:
    #                 await permission_service.has_permissions(
    #                     permission_state=permission_state,
    #                     member_snowflake=interaction.user.id,
    #                     requested=["other_guilds"],
    #                 )
    #             await permission_service.has_permissions(
    #                 permission_state=permission_state,
    #                 member_snowflake=interaction.user.id,
    #                 requested=["command.listing.scope.guild"],
    #             )
    #             guild_snowflake = guild.target.id
    #         else:
    #             return await tick.end(warning="This command must be used in a server.")
    #     if target is None:
    #         if interaction.channel is None:
    #             return await tick.end(
    #                 warning=f"This command must target a valid channel."
    #             )
    #         obj = interaction.channel
    #     else:
    #         obj = target.target
    #     if isinstance(obj, discord.Guild):
    #         await permission_service.has_permissions(
    #             permission_state=permission_state,
    #             member_snowflake=interaction.user.id,
    #             requested=["command.listing.scope.guild"],
    #         )
    #     elif isinstance(obj, discord.abc.GuildChannel):
    #         await permission_service.has_permissions(
    #             permission_state=permission_state,
    #             member_snowflake=interaction.user.id,
    #             requested=["command.listing.scope.channel"],
    #         )
    #     elif isinstance(obj, discord.Member):
    #         await permission_service.has_permissions(
    #             permission_state=permission_state,
    #             member_snowflake=interaction.user.id,
    #             requested=["command.listing.scope.member"],
    #         )
    #     pages = await list_coordinators.build_pages(
    #         guild_snowflake=guild_snowflake, obj=obj
    #     )
    #     return await tick.end(success=pages)

    @app_commands.command(name="flags", description="List flags.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_flags_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.scope.guild"],
                )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning=f"This command must be used in a server channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.guild"],
            )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.channel"],
            )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.member"],
            )
        pages = await list_flags.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    # @app_commands.command(name="mods", description="Lists mods.")
    # @app_commands.describe(
    #     target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    # )
    # async def list_moderators_app_command(
    #     self,
    #     interaction: discord.Interaction,
    #     target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    #     guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    # ) -> discord.Message:
    #     tick = Tick(bot=self.__bot, interaction=interaction)
    #     if guild is None:
    #         if interaction.guild is None:
    #             return await tick.end(
    #                 warning="This command must target a valid server."
    #             )
    #         guild_snowflake = interaction.guild.id
    #     else:
    #         if isinstance(guild.target, discord.Guild):
    #             guild_snowflake = guild.target.id
    #         else:
    #             return await tick.end(
    #                 warning="This command must target a valid server."
    #             )
    #     if target is None:
    #         if interaction.channel is None:
    #             return await tick.end(
    #                 warning=f"This command must target a valid channel."
    #             )
    #         obj = interaction.channel
    #     else:
    #         obj = target.target
    #     pages = await list_moderators.build_pages(
    #         guild_snowflake=guild_snowflake,
    #         obj=obj,
    #     )
    #     return await tick.end(success=pages)

    @app_commands.command(name="mutes", description="List mutes.")
    @app_commands.describe(
        target="Specify one of: a channel ID/mention, member ID/mention or server ID.",
        scope="Specify one of: `all`, `auto`, `click`, or `command`.",
        guild="Specify a server ID.",
    )
    async def list_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        scope: app_commands.Transform[ScopeObject | None, AppScope] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.scope.guild"],
                )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning=f"This command must used in a server channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.guild"],
            )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.channel"],
            )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.member"],
            )
        if scope is None:
            mute_type = "all"
        else:
            mute_type = scope.scope
        match mute_type:
            case "all":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=[
                        "command.listing.voice-mutes.auto"
                        "command.listing.voice-mutes.click"
                        "command.listing.voice-mutes.command"
                    ],
                )
            case "auto":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.auto"],
                )
            case "click":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.click"],
                )
            case "command":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.command"],
                )
        pages = await list_voice_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj, mute_type=mute_type
        )
        return await tick.end(success=pages)

    @app_commands.command(name="summary", description="List user moderation.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
    )
    async def list_moderation_summary_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        scope: app_commands.Transform[ScopeObject | None, AppScope] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.listing.scope.member"],
        )
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.scope.guild"],
                )
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
        if scope is None:
            mute_type = "all"
        else:
            mute_type = scope.scope
        match mute_type:
            case "all":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=[
                        "command.listing.voice-mutes.auto"
                        "command.listing.voice-mutes.click"
                        "command.listing.voice-mutes.command"
                    ],
                )
            case "auto":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.auto"],
                )
            case "click":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.click"],
                )
            case "command":
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.voice-mutes.command"],
                )
        pages: list[discord.Embed] = []
        services = []
        services.append(list_bans)
        services.append(list_flags)
        services.append(list_text_mutes)
        for service in services:
            summary_pages = await service.build_pages(
                guild_snowflake=guild_snowflake,
                obj=member_snowflake,
            )
            if isinstance(summary_pages, list):
                for page in summary_pages:
                    if isinstance(page, discord.Embed):
                        pages.append(page)
        summary_pages = await list_voice_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=member_snowflake, mute_type=mute_type
        )
        if isinstance(summary_pages, list):
            for page in summary_pages:
                if isinstance(page, discord.Embed):
                    pages.append(page)
        if not pages:
            return await tick.end(success="No infractions found")
        return await tick.end(success=pages)

    @app_commands.command(name="tmutes", description="List text-mutes.")
    @app_commands.describe(
        target="Specify a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_text_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(warning="This command must be used in a server.")
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be used in a server."
                    )
                if guild.target.id != interaction.guild.id:
                    await permission_service.has_permissions(
                        permission_state=permission_state,
                        member_snowflake=interaction.user.id,
                        requested=["other_guilds"],
                    )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    requested=["command.listing.scope.guild"],
                )
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning=f"This command must be used in a server channel."
                )
            obj = interaction.channel
        else:
            obj = target.target
        if isinstance(obj, discord.Guild):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.guild"],
            )
        elif isinstance(obj, discord.abc.GuildChannel):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.channel"],
            )
        elif isinstance(obj, discord.Member):
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                requested=["command.listing.scope.member"],
            )
        pages = await list_text_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    async def cog_app_command_error(self, interaction, error):
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))


async def setup(bot: DiscordBot):
    await bot.add_cog(ModeratorAppCommands(bot=bot))
