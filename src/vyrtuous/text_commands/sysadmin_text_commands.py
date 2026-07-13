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
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import developer_service, moderator_service


class SysadminTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Sysadmin"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx: commands.Context) -> bool:
        if ctx.guild is None:
            raise commands.CheckFailure("This command must be executed inside a server.")
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @commands.command(name="dev", help="Grant/revoke devs.")
    async def toggle_developer_text_command(
        self,
        ctx: commands.Context,
        member: Union[int, discord.Member] = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID",
        ),
    ) -> discord.Message:
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
            author_snowflake=ctx.author.id,
            member_snowflake=member_snowflake,
            message_snowflake=ctx.message.id,
            message_channel_snowflake=ctx.message.channel.id,
        )
        return await tick.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(SysadminTextCommands(bot=bot))
