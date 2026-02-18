"""!/bin/python3
coordinator_text_commands.py A discord.py cog containing coordinator commands for the Vyrtuous bot.

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

import discord
from discord.ext import commands

from vyrtuous.administrator.administrator_service import AdministratorService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.cog.help_command import skip_text_command_help_discovery
from vyrtuous.coordinator.coordinator import Coordinator
from vyrtuous.coordinator.coordinator_service import CoordinatorService, NotCoordinator
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration import Duration, DurationObject
from vyrtuous.moderator.moderator_service import ModeratorService
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.stage_room.stage_service import StageService

# from vyrtuous.field.snowflake import ChannelSnowflake, MemberSnowflake
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.discord_object_service import DiscordObjectService
from vyrtuous.utils.emojis import Emojis

# from vyrtuous.utils.permission_service import PermissionService
from vyrtuous.utils.state_service import StateService


class CoordinatorTextCommands(commands.Cog):
    ROLE = Coordinator

    def __init__(self, bot: DiscordBot):
        self.__author_service = AuthorService()
        self.__bot = bot
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__moderator_service = ModeratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__stage_service = StageService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
        )
        self.__bug_service = BugService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__developer_service = DeveloperService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            database_factory=self.__database_factory,
            emoji=self.__emoji,
        )
        self.__guild_owner_service = GuildOwnerService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__administrator_service = AdministratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__coordinator_service = CoordinatorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__sysadmin_service = SysadminService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__discord_object_service = DiscordObjectService()

    async def cog_check(self, ctx) -> Coroutine[Any, Any, bool]:
        async def predicate(
            source: Union[commands.Context, discord.Interaction, discord.Message],
        ):
            for verify in (
                self.__sysadmin_service.is_sysadmin_wrapper,
                self.__developer_service.is_developer_wrapper,
                self.__guild_owner_service.is_guild_owner_wrapper,
                self.__administrator_service.is_administrator_wrapper,
                self.__coordinator_service.is_coordinator_at_all_wrapper,
            ):
                try:
                    if await verify(source):
                        return True
                except commands.CheckFailure:
                    continue
            raise NotCoordinator

        predicate._permission_level = "Coordinator"
        return await predicate(ctx)

    @commands.command(name="mod", help="Grant/revoke mods.")
    async def toggle_moderator_text_command(
        self,
        ctx: commands.Context,
        member: discord.Member = commands.parameter(
            converter=commands.MemberConverter,
            description="Tag a member or include their ID",
        ),
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            description="Tag a channel or include its ID.",
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
        default_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        updated_kwargs = default_kwargs.copy()
        channel_dict = self.__discord_object_service.to_dict(obj=channel)
        member_dict = self.__discord_object_service.to_dict(obj=member)
        updated_kwargs.update(channel_dict.get("columns", None))
        # await PermissionService.has_equal_or_lower_role(
        #     target_member_snowflake=int(member.get("id", None)),
        #     **updated_kwargs,
        # )
        msg = await self.__moderator_service.toggle_moderator(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            member_dict=member_dict,
        )
        return await state.end(success=msg)

    @commands.command(name="stage", help="Start/stop stage")
    @skip_text_command_help_discovery()
    async def toggle_stage_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            description="Tag a channel or include its ID.",
        ),
        *,
        duration: Duration = commands.parameter(
            default=DurationObject("1h"),
            description="Options: (+|-)duration(m|h|d) 0 - permanent / 24h - default",
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
        default_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        obj = channel or ctx.channel
        channel_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__stage_service.toggle_stage(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            duration=duration,
        )
        return await state.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorTextCommands(bot))
