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

from typing import Any, Coroutine

import discord
from discord.ext import commands

from vyrtuous.active_members.active_member_service import ActiveMemberService
from vyrtuous.administrator.administrator_service import AdministratorService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.coordinator.coordinator import Coordinator
from vyrtuous.coordinator.coordinator_service import CoordinatorService, NotCoordinator
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_builder import DurationBuilder
from vyrtuous.moderator.help_text_command import skip_text_command_help_discovery
from vyrtuous.moderator.moderator_service import ModeratorService
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.stage_room.stage_service import StageService
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.upload.upload_service import UploadService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.data_service import DataService
from vyrtuous.utils.default_context import DefaultContext
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.message_service import PaginatorService
from vyrtuous.utils.state_service import StateService
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService
from vyrtuous.utils.discord_object_service import MultiConverter
from vyrtuous.ban.ban_service import BanService


class CoordinatorTextCommands(commands.Cog):
    ROLE = Coordinator

    def __init__(self, *, bot: DiscordBot | None = None):
        self.__author_service = AuthorService()
        self.__bot = bot
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__active_member_service = ActiveMemberService(
            bot=self.__bot, database_factory=self.__database_factory
        )
        self.__emoji = Emojis()
        self.__duration_builder = DurationBuilder()
        self.__paginator_service = PaginatorService(bot=self.__bot)
        self.__bug_service = BugService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__administrator_service = AdministratorService(
            active_member_service=self.__active_member_service,
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__coordinator_service = CoordinatorService(
            active_member_service=self.__active_member_service,
            author_service=self.__author_service,
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
        self.__guild_owner_service = GuildOwnerService(
            active_member_service=self.__active_member_service,
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__sysadmin_service = SysadminService(
            active_member_service=self.__active_member_service,
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__moderator_service = ModeratorService(
            active_member_service=self.__active_member_service,
            administrator_service=self.__administrator_service,
            author_service=self.__author_service,
            bot=self.__bot,
            coordinator_service=self.__coordinator_service,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            guild_owner_service=self.__guild_owner_service,
            sysadmin_service=self.__sysadmin_service,
        )
        self.__stream_service = StreamService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            paginator_service=self.__paginator_service,
        )
        self.__data_service = DataService(
            database_factory=self.__database_factory,
            duration_builder=self.__duration_builder,
            moderator_service=self.__moderator_service,
        )
        self.__voice_mute_service = VoiceMuteService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            stream_service=self.__stream_service,
        )
        self.__stage_service = StageService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            voice_mute_service=self.__voice_mute_service,
        )
        self.__upload_service = UploadService(
            bot=self.__bot, database_factory=self.__database_factory
        )
        self.__ban_service = BanService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )

    async def cog_check(self, ctx) -> Coroutine[Any, Any, bool]:
        context = DefaultContext(ctx=ctx)
        for verify in (
            self.__sysadmin_service.is_sysadmin_wrapper,
            self.__developer_service.is_developer_wrapper,
            self.__guild_owner_service.is_guild_owner_wrapper,
            self.__administrator_service.is_administrator_wrapper,
            self.__coordinator_service.is_coordinator_at_all_wrapper,
        ):
            try:
                if await verify(context=context):
                    return True
            except commands.CheckFailure:
                continue
        raise NotCoordinator

    cog_check._permission_level = "Coordinator"

    @commands.command(name="blacklist", help="Blacklist overwrite cleanup.")
    async def toggle_blacklist_text_command(
        self,
        ctx: commands.Context,
        member: int | discord.Member = commands.parameter(
            converter=MultiConverter,
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
            upload_service=self.__upload_service,
        )
        context = DefaultContext(ctx=ctx)
        await self.__moderator_service.check_minimum_role(
            channel_snowflake=channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role="Coordinator",
        )
        await self.__moderator_service.has_equal_or_lower_role(
            target_member_snowflake=int(member.id),
            member_snowflake=context.author.id,
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
        )
        if isinstance(member, discord.Member):
            member_snowflake = member.id
        else:
            member_snowflake = member
        msg = await self.__ban_service.toggle_blacklist(
            channel=channel,
            member_snowflake=member_snowflake,
        )
        return await state.end(success=msg)

    @commands.command(name="mod", help="Grant/revoke mods.")
    async def toggle_moderator_text_command(
        self,
        ctx: commands.Context,
        member: int | discord.Member = commands.parameter(
            converter=MultiConverter,
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
            upload_service=self.__upload_service,
        )
        context = DefaultContext(ctx=ctx)
        await self.__moderator_service.check_minimum_role(
            channel_snowflake=channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role="Coordinator",
        )
        await self.__moderator_service.has_equal_or_lower_role(
            target_member_snowflake=int(member.id),
            member_snowflake=context.author.id,
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
        )
        if isinstance(member, discord.Member):
            display_name = member.display_name
            member_snowflake = member.id
        else:
            display_name = None
            member_snowflake = member
        msg = await self.__moderator_service.toggle_moderator(
            channel=channel,
            display_name=display_name,
            member_snowflake=member_snowflake,
        )
        return await state.end(success=msg)

    @commands.command(name="purge", help="Delete messages.")
    async def purge_text_command(
        self,
        ctx: commands.Context,
        member: int | discord.Member = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID",
        ),
        amount: int = commands.parameter(
            default=25, description="Number of messages to delete"
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
        context = DefaultContext(ctx=ctx)
        await self.__moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role="Coordinator",
        )
        await self.__moderator_service.has_equal_or_lower_role(
            target_member_snowflake=int(member.id),
            member_snowflake=context.author.id,
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.channel.guild.id,
        )
        if isinstance(member, discord.Member):
            member_snowflake = int(member.id)
            display_name = str(member.mention)
        else:
            member_snowflake = int(member)
            member = self.__active_member_service.active_members.get(
                member_snowflake, None
            )
            if member:
                display_name = member.get("name", None)
            else:
                display_name = member_snowflake
        count = int(0)
        async for msg in ctx.channel.history():
            if amount == count:
                break
            if msg.author.id == member_snowflake:
                await msg.delete()
                count += 1
        return await state.end(
            success=f"Successfully deleted {count} messages from {display_name} in {ctx.channel.mention}."
        )

    @commands.command(name="stage", help="Start/stop stage")
    @skip_text_command_help_discovery()
    async def toggle_stage_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
        *,
        duration: str = commands.parameter(
            default="1h",
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
            upload_service=self.__upload_service,
        )
        context = DefaultContext(ctx=ctx)
        resolved_channel = channel or ctx.channel
        await self.__moderator_service.check_minimum_role(
            channel_snowflake=resolved_channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role="Coordinator",
        )
        pages = await self.__stage_service.toggle_stage(
            channel=resolved_channel, context=context, duration_value=duration
        )
        return await state.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorTextCommands(bot=bot))
