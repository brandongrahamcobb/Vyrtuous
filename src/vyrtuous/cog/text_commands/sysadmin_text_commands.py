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

from typing import Any, Coroutine, Union

from discord.ext import commands
import discord

from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.sysadmin.sysadmin import Sysadmin
from vyrtuous.sysadmin.sysadmin_service import SysadminService, NotSysadmin
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.utils.discord_object_service import DiscordObjectService
from vyrtuous.utils.message_service import MessageService
from vyrtuous.utils.state_service import StateService


class SysadminTextCommands(commands.Cog):
    ROLE = Sysadmin

    def __init__(self, bot: DiscordBot):
        self.__bot = bot
        self.message_service = MessageService(self.__bot)
        self.__bug_service = BugService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__discord_object_service = DiscordObjectService()
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__developer_service = DeveloperService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            emoji=self.__emoji,
        )
        self.__sysadmin_service = SysadminService(
            author_service=self.__author_service,
            bot=self.__bot,
        )

    async def cog_check(self, ctx) -> Coroutine[Any, Any, bool]:
        async def predicate(
            source: Union[commands.Context, discord.Interaction, discord.Message],
        ):
            if await self.__sysadmin_service.is_sysadmin_wrapper(source):
                return True
            raise NotSysadmin

        predicate._permission_level = "Sysadmin"
        return predicate(ctx)

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
        state = StateService(ctx=ctx)
        member_dict = self.__discord_object_service.translate(obj=member)
        embed = await self.__bug_service.assign_bug_to_developer(
            reference=reference, member_dict=member_dict
        )
        return await state.end(success=embed)

    @commands.command(name="dev", help="Grant/revoke devs.")
    async def toggle_developer_text_command(
        self,
        ctx: commands.Context,
        member: discord.Member = commands.parameter(
            converter=commands.MemberConverter,
            description="Tag a member or include their ID",
        ),
    ):
        state = StateService(ctx=ctx)
        default_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        member_dict = self.__discord_object_service.translate(obj=member)
        msg = await self.__developer_service.toggle_developer(
            member_dict=member_dict,
        )
        return await state.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(SysadminTextCommands(bot))
