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
from vyrtuous.cache.permissions import PermissionGroup
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.models.group import Group, GroupObject
from vyrtuous.models.metadata import metadata
from vyrtuous.models.target import Target, TargetObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import autoassign_role_service, hero_service, vegan_service


class UserManagementTextCommands(commands.Cog):

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    @commands.command(name="autoassign", help="Toggle an autoassign role.")
    @metadata(permission="command.users.autoassign")
    async def toggle_autoassign_role_text_command(
        self,
        ctx: commands.Context,
        group: GroupObject = commands.parameter(
            converter=Group,
            description="Specify a group alias or name.",
        ),
        role: TargetObject = commands.parameter(
            converter=Target,
            description="Specify a role ID/mention.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if guild is None:
            guild_snowflake = ctx.guild.id
        else:
            if not isinstance(guild.target, discord.Guild):
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = guild.target.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
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
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
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
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.users.autoassign"],
        )
        if ctx.guild.id != guild_snowflake:
            await permission_service.has_permissions(
                permission_state=permission_state,
                member_snowflake=ctx.author.id,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                requested=["other_guilds"],
            )
        pages = await autoassign_role_service.toggle_autoassign_role(
            author_snowflake=ctx.author.id,
            group=group.group,
            guild_snowflake=guild_snowflake,
            role_snowflake=role_snowflake,
        )
        return await tick.end(success=pages)

    @commands.command(name="hero", help="Grant/revoke invincibility.")
    @metadata(permission="command.users.hero")
    async def toggle_hero_text_command(
        self,
        ctx: commands.Context,
        member: TargetObject = commands.parameter(
            converter=Target,
            description="Specify a member ID/mention.",
        ),
        target: TargetObject | None = commands.parameter(
            converter=Target,
            default=None,
            description="Specify a `all` or a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        singular = True
        if target is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guilds = [ctx.guild]
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
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
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
                member_snowflake=ctx.author.id,
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

    @commands.command(name="vegan", help="Toggle vegan.")
    @metadata(permission="command.users.vegan")
    async def toggle_vegan_text_command(
        self,
        ctx: commands.Context,
        member: TargetObject = commands.parameter(
            converter=Target,
            description="Specify a member ID/mention.",
        ),
        guild: TargetObject | None = commands.parameter(
            converter=TargetObject,
            default=None,
            description="Specify a server ID.",
        ),
        notes: str | None = commands.parameter(
            default=None,
            description="Include notes.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
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
        if not vegan_service.is_vegan(
            guild_snowflake=guild_snowflake, member_snowflake=ctx.author.id
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
            member_snowflake=ctx.author.id,
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
    await bot.add_cog(UserManagementTextCommands(bot=bot))
