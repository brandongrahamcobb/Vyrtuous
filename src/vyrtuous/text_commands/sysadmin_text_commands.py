"""!/bin/python3
sysadmin_text_commands.py A discord.py cog containing sysadmin commands for the Vyrtuous bot.

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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.cache.sysadmin import NotSysadmin
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.upload import upload_service
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import developer_service, sysadmin_service


class SysadminTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Sysadmin"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx: commands.Context):
        if ctx.guild is None:
            raise commands.CheckFailure("This command must be used inside a server.")
        context = SnowflakeContext(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
        )
        if await sysadmin_service.is_sysadmin_wrapper(context=context):
            return True
        raise NotSysadmin

    @commands.command(name="assign", help="Assign developer.")
    async def assign_bug_to_developer_text_command(
        self,
        ctx: commands.Context,
        reference: str = commands.parameter(
            description="Include an issue reference ID"
        ),
        member: discord.Member = commands.parameter(
            converter=commands.MemberConverter,
            description="Tag a member or include their ID",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        embed = await developer_service.handle_developer_assignment(
            member=member,
            reference=reference,
        )
        return await tick.end(success=embed)

    @commands.command(name="dev", help="Grant/revoke devs.")
    async def toggle_developer_text_command(
        self,
        ctx: commands.Context,
        member: Union[int, discord.Member] = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if isinstance(member, discord.Member):
            member_snowflake = member.id
        else:
            simplified_member = self.__bot.registry.get(MemberState).active.get(
                member, None
            )
            if not simplified_member:
                raise commands.MemberNotFound(str(member))
            else:
                member_snowflake = member
        msg = await developer_service.toggle_developer(
            member_snowflake=member_snowflake,
        )
        return await tick.end(success=msg)

    @commands.command(name="upload", help="Create the upload document.")
    async def uploads_text_command(
        self,
        ctx: commands.Context,
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        await upload_service.build_latex_document()
        return await tick.end(success="Success!")


async def setup(bot: DiscordBot):
    await bot.add_cog(SysadminTextCommands(bot=bot))
