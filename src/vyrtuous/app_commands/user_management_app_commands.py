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
from vyrtuous.cache.permissions import PermissionGroup
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.models.group import AppGroup, GroupObject
from vyrtuous.models.metadata import metadata
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import autoassign_role_service, hero_service, vegan_service
from vyrtuous.view.grant_view import GrantView
from vyrtuous.view.revoke_view import RevokeView


class UserManagementAppCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def cog_app_command_error(self, interaction, error):
        self.__bot.logger.info(str(error))
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))

    @metadata(permission="command.users.autoassign")
    @app_commands.command(name="autoassign", description="Toggle an autoassign role.")
    @app_commands.describe(
        group="Specify a group alias/name.",
        role="Specify a role ID/mention.",
        channel="Specify a channel ID/mention.",
        guild="Specify a server ID.",
    )
    async def toggle_autoassign_role_app_command(
        self,
        interaction: discord.Interaction,
        group: app_commands.Transform[GroupObject, AppGroup],
        role: app_commands.Transform[TargetObject, AppTarget],
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        if guild is None:
            guild_snowflake = interaction.guild.id
        else:
            if not isinstance(guild.target, discord.Guild):
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = guild.target.id
        if channel is None:
            channel_snowflake = None
        else:
            if not isinstance(channel.target, discord.abc.GuildChannel):
                return await tick.end(
                    warning="This command must be used in a server channel."
                )
            channel_snowflake = channel.target.id
        if not isinstance(role.target, discord.Role):
            return await tick.end(warning="This command must target a valid role.")
        else:
            role_snowflake = role.target.id
        if not isinstance(group.group, PermissionGroup):
            return await tick.end(warning="This command must target a valid group.")
        else:
            group_obj = group.group
            effective_group = await permission_service.resolve_effective_group(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake or interaction.channel.id,
                guild_snowflake=guild_snowflake,
            )
            if effective_group:
                if group_obj.alias not in effective_group.ancestors:
                    return await tick.end(
                        warning="You cannot autoassign a group you do not inherit."
                    )
            else:
                return await tick.end(
                    warning="You cannot autoassign a group you do not inherit."
                )
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake or interaction.channel.id,
            guild_snowflake=guild_snowflake,
            requested=["command.users.autoassign"],
        )
        if interaction.guild.id != guild_snowflake:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=interaction.user.id,
                channel_snowflake=channel_snowflake or interaction.channel.id,
                guild_snowflake=guild_snowflake,
                requested=["other_guilds"],
            )
        pages = await autoassign_role_service.toggle_autoassign_role(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            group=group_obj,
            guild_snowflake=guild_snowflake,
            role_snowflake=role_snowflake,
        )
        return await tick.end(success=pages)

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
        ctx = SnowflakeContext(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = GrantView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            interaction=interaction,
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
        ctx = SnowflakeContext(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        view = RevokeView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            interaction=interaction,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Specify the group", view=view, ephemeral=True
        )

    @metadata(permission="command.users.vegan")
    @app_commands.command(name="vegan", description="Toggle vegan.")
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
            requested=["command.users.vegan"],
        )
        embed = await vegan_service.toggle_vegan(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            notes=notes or "No notes provided.",
        )
        return await tick.end(success=embed)


async def setup(bot: DiscordBot):
    await bot.add_cog(UserManagementAppCommands(bot=bot))
