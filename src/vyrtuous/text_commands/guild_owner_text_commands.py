"""!/bin/python3
guild_owner_text_commands.py A discord.py cog containing guild owner commands for the Vyrtuous bot.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from typing import Literal, Optional, Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.listing import list_developers
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import administrator_role_service, moderator_service


class GuildOwnerTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Guild Owner"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx) -> bool:
        if ctx.guild is None:
            raise commands.CheckFailure(
                "This command must be executed inside a server."
            )
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @commands.command(name="admin", help="Toggle administrator role.")
    async def toggle_administrator_by_role_text_command(
        self,
        ctx: commands.Context,
        role: discord.Role = commands.parameter(
            converter=commands.RoleConverter, description="Tag a role or its ID"
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        pages = await administrator_role_service.toggle_administrator_role(
            author_snowflake=ctx.author.id,
            guild_snowflake=ctx.guild.id,
            message_snowflake=ctx.message.id,
            message_channel_snowflake=ctx.message.channel.id,
            role_snowflake=role.id,
        )
        return await tick.end(success=pages)

    @commands.command(name="devs", help="List devs.")
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        member: Union[int, discord.Member, None] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a member ID/mention.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        pages = await list_developers.build_pages(obj=member)
        return await tick.end(success=pages)

    @commands.command(name="sync", help="Sync app commands.")
    async def sync_text_command(
        self,
        ctx: commands.Context,
        spec: Optional[Literal["~", "*", "^"]] = None,
        *,
        guilds: Union[commands.Greedy[discord.Object], None] = None,
    ) -> discord.Message:
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
