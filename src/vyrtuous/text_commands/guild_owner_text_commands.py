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

from typing import Literal, Optional, Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.guild_owner import NotGuildOwner
from vyrtuous.cache.registry import MemberState
from vyrtuous.inc.helpers import at_home
from vyrtuous.listing import list_developers, list_heroes
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import (
    administrator_role_service,
    developer_service,
    guild_owner_service,
    hero_service,
    sysadmin_service,
)


class GuildOwnerTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Guild Owner"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx):
        if ctx.guild is None:
            raise commands.CheckFailure("This command must be used inside a server.")
        context = SnowflakeContext(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
        )
        for verify in (
            sysadmin_service.is_sysadmin_wrapper,
            developer_service.is_developer_wrapper,
            guild_owner_service.is_guild_owner_wrapper,
        ):
            try:
                if await verify(context=context):
                    return True
            except commands.CheckFailure:
                continue
        raise NotGuildOwner

    @commands.command(name="admin", help="Toggle administrator role.")
    async def toggle_administrator_by_role_text_command(
        self,
        ctx: commands.Context,
        role: discord.Role = commands.parameter(
            converter=commands.RoleConverter, description="Tag a role or its ID"
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        pages = await administrator_role_service.toggle_administrator_role(
            role=role,
        )
        return await tick.end(success=pages)

    @commands.command(name="hero", help="Grant/revoke invincibility.")
    @skip_text_command_help_discovery()
    async def invincibility_text_command(
        self,
        ctx: commands.Context,
        member: discord.Member = commands.parameter(
            converter=commands.MemberConverter,
            description="Tag a member or include their ID",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if member.id in (self.__bot.registry.get(MemberState).invincible[ctx.guild.id]):
            await hero_service.add_invincible_member(member.guild.id, member.id)
            await hero_service.unrestrict(
                guild_snowflake=member.guild.id, member_snowflake=member.id
            )
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member.mention}."
            )
        else:
            await hero_service.remove_invincible_member(member.guild.id, member.id)
            self.__bot.registry.get(MemberState).invincible[ctx.guild.id].remove(
                member.id
            )
            msg = f"Invincibility has been disabled for {member.mention}"
        return await tick.end(success=msg)

    @commands.command(name="devs", help="List devs.")
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[str, discord.Member, None] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="'all', or user mention/ID",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        obj = target or "all"
        pages = await list_developers.build_pages(obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="heroes", help="List heroes.")
    async def list_heroes_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[str, discord.Member] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="'all', or user mention/ID",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        is_at_home = at_home(source=ctx)
        pages = await list_heroes.build_pages(
            is_at_home=is_at_home,
            obj=target,
        )
        return await tick.end(success=pages)

    @commands.command(name="sync", help="Sync app commands.")
    async def sync_text_command(
        self,
        ctx: commands.Context,
        spec: Optional[Literal["~", "*", "^"]] = None,
        *,
        guilds: Union[commands.Greedy[discord.Object], None] = None,
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
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
                return await tick.end(success=msg)
            except Exception as e:
                return await tick.end(warning=str(e).capitalize())
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        return await tick.end(success=f"Synced the tree to {ret}/{len(guilds)}.")


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildOwnerTextCommands(bot=bot))
