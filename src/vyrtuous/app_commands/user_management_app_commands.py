"""!/bin/python3
permission_manager_app_commands.py A discord.py cog containing permission manager commands for the Vyrtuous bot.

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
from vyrtuous.models.metadata import metadata
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import hero_service, vegan_service
from vyrtuous.view.grant_view import GrantView
from vyrtuous.view.revoke_view import RevokeView
from vyrtuous.view.view_context import ViewContext


class UserManagementAppCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    @metadata(permission="command.users.grant")
    @app_commands.command(name="grant", description="Grant permission levels.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def grant_group_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        else:
            guild_snowflake = interaction.guild.id
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
            requested=["command.users.grant"],
        )
        ctx = ViewContext(
            interaction=interaction,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = GrantView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the group", view=view, ephemeral=True
        )

    @metadata(permission="command.users.hero")
    @app_commands.command(name="hero", description="Grant/revoke invincibility.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
        target="Specify `all` or a server ID.",
    )
    async def toggle_hero_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        singular = True
        if target is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guilds = [interaction.guild]
        else:
            if isinstance(target.target, discord.Guild):
                guilds = [target.target]
            elif target.target == "all":
                guilds = bot.guilds
                singular = False
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
            member_name = "Unknown"
            simplified_member = self.__bot.registry.get(MemberState).active.get(
                member_snowflake
            )
            if simplified_member:
                member_name = simplified_member[0]
        elif isinstance(member.target, discord.Member):
            member_snowflake = member.target.id
            member_name = member.target.mention
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        if not singular:
            if any(
                member_snowflake in member_set
                for member_set in self.__bot.registry.get(
                    MemberState
                ).invincible.values()
            ):
                enabled = False
            else:
                enabled = True
        else:
            guild_snowflake = guilds[0].id
            if member_snowflake in self.__bot.registry.get(MemberState).invincible.get(
                guild_snowflake, set()
            ):
                enabled = False
            else:
                enabled = True
        if enabled:
            header = f"All moderation events have been forgiven and invincibility has been enabled for {member_name} in "
        else:
            header = f"Invincibility has been disabled for {member_name} in "
        guild_names = []
        for g in guilds:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=g.id,
                requested=["command.users.hero"],
            )
            if enabled:
                await hero_service.add_invincible_member(g.id, member_snowflake)
                await hero_service.unrestrict(
                    guild_snowflake=g.id, member_snowflake=member_snowflake
                )
                guild_names.append(g.name)
            else:
                await hero_service.remove_invincible_member(g.id, member_snowflake)
                guild_names.append(g.name)
        msg = header + ", ".join(guild_names) + "."
        return await tick.end(success=msg)

    @metadata(permission="command.users.revoke")
    @app_commands.command(name="revoke", description="Revoke permission levels.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
    )
    async def revoke_group_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        else:
            guild_snowflake = interaction.guild.id
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
            requested=["command.users.revoke"],
        )
        ctx = ViewContext(
            interaction=interaction,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = RevokeView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the group", view=view, ephemeral=True
        )

    # @app_commands.command(name="survey", description="Survey stage members.")
    # @app_commands.describe(
    #     channel="Specify channel ID/mention.",
    # )
    # async def survey_app_command(
    #     self,
    #     interaction: discord.Interaction,
    #     channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    # ) -> discord.Message:
    #     tick = Tick(bot=self.__bot, interaction=interaction)
    #     if channel is None:
    #         if interaction.guild is None:
    #             return await tick.end(
    #                 warning="This command must target a valid server."
    #             )
    #         if interaction.channel is None:
    #             return await tick.end(
    #                 warning="This command must target a valid channel."
    #             )
    #         channel_snowflake = interaction.channel.id
    #         guild_snowflake = interaction.guild.id
    #     elif isinstance(
    #         channel.target,
    #         (discord.VoiceChannel, discord.StageChannel),
    #     ):
    #         channel_snowflake = channel.target.id
    #         guild_snowflake = channel.target.guild.id
    #     else:
    #         return await tick.end(warning="This command must target a valid channel.")
    #     pages = await moderator_service.survey(
    #         channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake
    #     )
    #     return await tick.end(success=pages)

    @metadata(permission="command.users.vegan")
    @app_commands.command(name="vcow", description="Toggle vegan.")
    @app_commands.describe(member="Specify member ID/mention.", notes="Include notes.")
    async def toggle_vegan_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        notes: str | None,
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
        if not vegan_service.is_vegan(
            guild_snowflake=guild_snowflake, member_snowflake=interaction.user.id
        ):
            return await tick.end(warning="Author is not a vegan.")
        if isinstance(member, int):
            member_snowflake = member
        elif isinstance(member, discord.Member):
            member_snowflake = member.id
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.users.vegan"],
        )
        embed = await vegan_service.toggle_vegan(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            notes=notes,
        )
        return await tick.end(success=embed)

    async def cog_app_command_error(self, interaction, error):
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))


async def setup(bot: DiscordBot):
    await bot.add_cog(UserManagementAppCommands(bot=bot))
