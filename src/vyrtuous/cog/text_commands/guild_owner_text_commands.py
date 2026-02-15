"""!/bin/python3
guild_owner_text_commands.py A discord.py cog containing guild owner commands for the Vyrtuous bot.

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

from typing import Literal, Optional

import discord
from discord.ext import commands

from vyrtuous.administrator.administrator_service import AdministratorRoleService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cog.help_command import skip_help_discovery
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.field.snowflake import MemberSnowflake, RoleSnowflake
from vyrtuous.owner.guild_owner import GuildOwner
from vyrtuous.owner.guild_owner_service import guild_owner_predicator
from vyrtuous.utils.discord_object_service import DiscordObject
from vyrtuous.utils.message_service import MessageService
from vyrtuous.utils.permission_service import PermissionService
from vyrtuous.utils.state_service import StateService


class GuildOwnerTextCommands(commands.Cog):
    ROLE = GuildOwner

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    @commands.command(name="admin", help="Toggle administrator role.")
    @guild_owner_predicator()
    async def toggle_administrator_by_role_text_command(
        self, ctx: commands.Context, role: RoleSnowflake
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        default_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        do = DiscordObject(ctx=ctx)
        role_dict = await do.determine_from_target(target=role)
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(role_dict.get("columns", None))
        pages = await AdministratorRoleService.toggle_administrator_role(
            role_dict=role_dict, updated_kwargs=updated_kwargs
        )
        await StateService.send_pages(title="Administrators", pages=pages, state=state)

    @commands.command(name="hero", help="Grant/revoke invincibility.")
    @guild_owner_predicator()
    @skip_help_discovery()
    async def invincibility_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            description="Tag a member or include their ID"
        ),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        do = DiscordObject(ctx=ctx)
        member_dict = await do.determine_from_target(target=member)
        where_kwargs = member_dict.get("columns", None)
        enabled = PermissionService.toggle_enabled()
        if enabled:
            PermissionService.add_invincible_member(**where_kwargs)
            await PermissionService.unrestrict(**where_kwargs)
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member_dict.get('mention', None)}."
            )
        else:
            PermissionService.remove_invincible_member(**where_kwargs)
            msg = f"Invincibility has been disabled for {member_dict.get('mention', None)}"
        return await state.end(success=msg)

    @commands.command(name="devs", help="List devs.")
    @guild_owner_predicator()
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        *,
        target: str | None = commands.parameter(
            default=None, description="'all', a specific server or user mention/ID"
        ),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        do = DiscordObject(ctx=ctx)
        target = target or "all"
        object_dict = await do.determine_from_target(target=target)
        pages = await DeveloperService.build_pages(object_dict=object_dict)
        await StateService.send_pages(title="Developers", pages=pages, state=state)

    @commands.command(name="sync", help="Sync app commands.")
    @guild_owner_predicator()
    async def sync_text_command(
        self,
        ctx: commands.Context,
        spec: Optional[Literal["~", "*", "^"]] = None,
        *,
        guilds: commands.Greedy[discord.Object] = None,
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        synced = []
        if not guilds:
            if spec == "~":
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == "*":
                ctx.bot.tree.copy_global_to(guild=ctx.guild)
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == "^":
                ctx.bot.tree.clear_commands(guild=ctx.guild)
                await ctx.bot.tree.sync(guild=ctx.guild)
            else:
                synced = await ctx.bot.tree.sync()
            try:
                if spec is None:
                    msg = f"Synced {len(synced)} commands globally."
                else:
                    msg = f"Synced {len(synced)} commands to the current server."
                return await state.end(success=msg)
            except Exception as e:
                return await state.end(warning=str(e).capitalize())
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        return await state.end(success=f"Synced the tree to {ret}/{len(guilds)}.")


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildOwnerTextCommands(bot))
