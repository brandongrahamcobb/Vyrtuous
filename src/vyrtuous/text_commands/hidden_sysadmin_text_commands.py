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

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.upload import upload_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import developer_service, moderator_service


class SysadminTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Sysadmin"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx: commands.Context) -> bool:
        if ctx.guild is None:
            raise commands.CheckFailure("This command must be used inside a server.")
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @commands.command(name="assign", help="Assign developer.")
    @skip_text_command_help_discovery()
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
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        embed = await developer_service.handle_developer_assignment(
            member=member,
            reference=reference,
        )
        return await tick.end(success=embed)

    @commands.command(name="upload", help="Create the upload document.")
    @skip_text_command_help_discovery()
    async def uploads_text_command(
        self,
        ctx: commands.Context,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        await upload_service.build_latex_document()
        return await tick.end(success="Success!")


async def setup(bot: DiscordBot):
    await bot.add_cog(SysadminTextCommands(bot=bot))
