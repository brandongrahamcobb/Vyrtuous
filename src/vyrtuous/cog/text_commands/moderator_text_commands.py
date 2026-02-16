"""!/bin/python3
moderator_text_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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
from vyrtuous.alias.alias_service import AliasService
from vyrtuous.ban.ban_service import BanService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService

# from vyrtuous.cog.help_command import skip_help_discovery
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_service import DurationService
from vyrtuous.flag.flag_service import FlagService
from vyrtuous.moderator.moderator import Moderator
from vyrtuous.moderator.moderator_service import ModeratorService, NotModerator
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.stage_room.stage_service import StageService
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.temporary_room.temporary_room_service import TemporaryRoomService
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.data_service import DataService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.discord_object_service import DiscordObjectService, MultiConverter
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.home import at_home
from vyrtuous.utils.state_service import StateService
from vyrtuous.vegan.vegan_service import VeganService
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService


class ModeratorTextCommands(commands.Cog):
    ROLE = Moderator

    def __init__(self, bot: DiscordBot):
        self.config = bot.config
        self.__author_service = AuthorService()
        self.__bot = bot
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__duration_service = DurationService()
        self.__sysadmin_service = SysadminService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__bug_service = BugService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__alias_service = AliasService(
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
        )
        self.__developer_service = DeveloperService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            emoji=self.__emoji,
        )
        self.__coordinator_service = CoordinatorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__administrator_service = AdministratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__discord_object_service = DiscordObjectService()
        self.__temporary_room_service = TemporaryRoomService(
            alias_service=self.__alias_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__moderator_service = ModeratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__data_service = DataService(
            duration_service=self.__duration_service,
            moderator_service=self.__moderator_service,
        )
        self.__stream_service = StreamService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
        )
        self.__voice_mute_service = VoiceMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__moderator_service = ModeratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__ban_service = BanService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
        )
        self.__flag_service = FlagService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__vegan_service = VeganService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__text_mute_service = TextMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
        )
        self.__guild_owner_service = GuildOwnerService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )

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
                self.__moderator_service.is_moderator_at_all_wrapper,
            ):
                try:
                    if await verify(source):
                        return True
                except commands.CheckFailure:
                    continue
            raise NotModerator

        predicate._permission_level = "Moderator"
        return await predicate(ctx)

    @commands.command(name="admins", help="Lists admins.")
    # @skip_help_discovery()
    async def list_administrators_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: `all`, channel ID/mention or server ID.",
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
        obj = target or ctx.guild
        is_at_home = at_home(source=ctx)
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__administrator_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="bans", help="List bans.")
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention or server ID.",
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
        obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__ban_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="cmds", help="List aliases.")
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
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
        obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__alias_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="coords", help="Lists coords.")
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: `all`, channel ID/mention, or server ID.",
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
        obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__coordinator_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="del", help="Delete message.")
    # @skip_help_discovery()
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        message: discord.Message = commands.parameter(
            converter=commands.MessageConverter, description="Message snowflake"
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
        try:
            msg = await ctx.channel.fetch_message(message)
        except discord.NotFound:
            return await state.end(warning=f"Message `{message}` not found.")
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())
        return await state.end(success=f"Message `{message}` deleted successfully.")

    @commands.command(name="flags", help="List flags.")
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
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
        obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__flag_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="ls", help="List new vegans.")
    # @skip_help_discovery()
    async def list_new_vegans_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
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
        obj = target or ctx.guild
        is_at_home = at_home(source=ctx)
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__vegan_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(
        name="migrate",
        help="Migrate a temporary room to a new channel by snowflake.",
    )
    # @skip_help_discovery()
    async def migrate_temp_room_text_command(
        self,
        ctx: commands.Context,
        old_name: str = commands.parameter(description="Provide a channel name"),
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
        channel_dict = self.__discord_object_service.to_dict(obj=channel)
        msg = await self.__temporary_room_service.migrate_temporary_room(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            old_name=old_name,
        )
        return await state.end(success=msg)

    @commands.command(name="mods", help="Lists mods.")
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
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
        obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__moderator_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="mutes", help="List mutes.")
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
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
        obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__voice_mute_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="mstage", help="Toggle stage mute/unmute.")
    # @skip_help_discovery()
    async def stage_mute_text_command(
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
        channel = channel or ctx.channel
        channel_dict = self.__discord_object_service.to_dict(obj=channel)
        member_dict = self.__discord_object_service.to_dict(obj=member)
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(channel_dict.get("columns", None))
        updated_kwargs.update(member_dict.get("columns", None))
        msg = await self.__stage_service.toggle_stage_mute(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            member_dict=member_dict,
        )
        return await state.end(success=msg)

    # @commands.command(name="summary", help="List user moderation.")
    # async def list_moderation_summary_text_command(
    #     self,
    #     ctx: commands.Context,
    #     member: discord.Member = commands.parameter(
    #         converter=commands.MemberConverter,
    #         description="Specify a member ID/mention."
    #     ),
    # ):
    #     state = StateService(
    #         author_service=self.__author_service,
    #         bot=self.__bot,
    #         bug_service=self.__bug_service,
    #         ctx=ctx,
    #         developer_service=self.__developer_service,
    #         emoji=self.__emoji,
    #     )
    #     pages = []
    #     dir_paths = []
    #     dir_paths.append(Path("/app/vyrtuous/db/infractions"))
    #     obj = member or int(ctx.member)
    #     is_at_home = at_home(source=ctx)
    #     member_dict = self.__discord_object_service.to_dict(obj=obj)
    #     for obj in dir_to_classes(dir_paths=dir_paths, parent=AliasService):
    #         object_pages = await obj.build_pages(
    #             object_dict=member_dict, is_at_home=is_at_home
    #         )
    #         if object_pages:
    #             pages.extend(object_pages)
    #     await StateService.send_pages(title="Infractions", pages=pages, state=state)

    @commands.command(name="survey", help="Survey stage members.")
    # @skip_help_discovery()
    async def stage_survey_text_command(
        self,
        ctx: commands.Context,
        channel: Union[discord.abc.GuildChannel, None] = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
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
        obj = channel or ctx.channel
        channel_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__moderator_service.survey(
            channel_dict=channel_dict, guild_snowflake=ctx.guild.id
        )
        return await state.end(success=pages)

    @commands.command(name="tmutes", help="List text-mutes.")
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
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
        obj = target or int(ctx.channel.id)
        is_at_home = at_home(source=ctx)
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__text_mute_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)


async def setup(bot: DiscordBot):
    cog = ModeratorTextCommands(bot)
    await bot.add_cog(cog)
