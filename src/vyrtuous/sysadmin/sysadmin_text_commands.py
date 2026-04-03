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

from typing import Any, Coroutine

import discord
from discord.ext import commands

from vyrtuous.active_members.active_member_service import ActiveMemberService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_builder import DurationBuilder
from vyrtuous.sysadmin.sysadmin import Sysadmin
from vyrtuous.sysadmin.sysadmin_service import NotSysadmin, SysadminService
from vyrtuous.upload.upload_service import UploadService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.default_context import DefaultContext
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.discord_object_service import MultiConverter
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.state_service import StateService


class SysadminTextCommands(commands.Cog):
    ROLE = Sysadmin

    def __init__(self, *, bot: DiscordBot | None = None):
        self.__bot = bot
        self.__emoji = Emojis()
        self.__author_service = AuthorService()
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__active_member_service = ActiveMemberService(
            bot=self.__bot, database_factory=self.__database_factory
        )
        self.__duration_builder = DurationBuilder()
        self.__bug_service = BugService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__developer_service = DeveloperService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            database_factory=self.__database_factory,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
        )
        self.__sysadmin_service = SysadminService(
            active_member_service=self.__active_member_service,
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__upload_service = UploadService(
            bot=self.__bot, database_factory=self.__database_factory
        )

    async def cog_check(self, ctx: commands.Context) -> Coroutine[Any, Any, bool]:
        context = DefaultContext(ctx=ctx)
        if await self.__sysadmin_service.is_sysadmin_wrapper(context=context):
            return True
        raise NotSysadmin

    cog_check._permission_level = "Sysadmin"

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
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            upload_service=self.__upload_service,
        )
        embed = await self.__developer_service.assign_bug_to_developer(
            member=member,
            reference=reference,
        )
        return await state.end(success=embed)

    @commands.command(name="dev", help="Grant/revoke devs.")
    async def toggle_developer_text_command(
        self,
        ctx: commands.Context,
        member: int | discord.Member = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID",
        ),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            upload_service=self.__upload_service,
        )
        if isinstance(member, discord.Member):
            display_name = member.display_name
            member_snowflake = member.id
        else:
            display_name = None
            member_snowflake = member
        try:
            msg = await self.__developer_service.toggle_developer(
                display_name=display_name,
                member_snowflake=member_snowflake,
            )
        except:
            import traceback

            traceback.print_exc()
        return await state.end(success=msg)

    @commands.command(name="upload", help="Create the upload document.")
    async def uploads_text_command(
        self,
        ctx: commands.Context,
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            upload_service=self.__upload_service,
        )
        await self.__upload_service.build_latex_document()
        return await state.end(success="Success!")


async def setup(bot: DiscordBot):
    await bot.add_cog(SysadminTextCommands(bot=bot))
