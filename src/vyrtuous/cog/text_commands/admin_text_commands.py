"""!/bin/python3
admin_text_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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

from vyrtuous.administrator.administrator import Administrator
from vyrtuous.administrator.administrator_service import (
    AdministratorRoleService, AdministratorService, NotAdministrator)
from vyrtuous.alias.alias_service import AliasService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.cap.cap_service import CapService
# from vyrtuous.cog.help_command import skip_help_discovery
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_service import DurationService
from vyrtuous.field.category import Category
from vyrtuous.inc.helpers import PATH_LOG
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.server_mute.server_mute_service import ServerMuteService
from vyrtuous.stage_room.stage_service import StageService
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.temporary_room.temporary_room_service import TemporaryRoomService
from vyrtuous.utils.author_service import AuthorService
# from vyrtuous.utils.clear_service import ClearService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.discord_object_service import (DiscordObjectService,
                                                   MultiConverter)
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.home import at_home
from vyrtuous.utils.logger import logger
from vyrtuous.utils.message_service import MessageService
# from vyrtuous.utils.permission_service import PermissionService
from vyrtuous.utils.state_service import StateService
from vyrtuous.video_room.video_room_service import VideoRoomService
# from vyrtuous.view.cancel_confirm_view import VerifyView
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService


class AdminTextCommands(commands.Cog):
    ROLE = Administrator

    def __init__(self, bot: DiscordBot):
        self.__author_service = AuthorService()
        self.__bot = bot
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__duration_service = DurationService()
        self.__alias_service = AliasService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__bug_service = BugService(
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
        self.__guild_owner_service = GuildOwnerService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
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
        self.__sysadmin_service = SysadminService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.message_service = MessageService(self.__bot)
        self.__administrator_role_service = AdministratorRoleService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__cap_service = CapService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
        )
        self.__discord_object_service = DiscordObjectService()
        self.__server_mute_service = ServerMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__temporary_room_service = TemporaryRoomService(
            alias_service=self.__alias_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__stream_service = StreamService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__video_room_service = VideoRoomService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__voice_mute_service = VoiceMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
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
            ):
                try:
                    if await verify(source):
                        return True
                except commands.CheckFailure:
                    continue
            raise NotAdministrator

        predicate._permission_level = "Administrator"
        return await predicate(ctx)

    @commands.command(
        name="alias",
        help="Alias creation.",
    )
    async def create_alias_text_command(
        self,
        ctx: commands.Context,
        category: Category = commands.parameter(
            description="Specify a category for a `ban`, `flag`, `role`, `tmute`, `vegan` or `vmute` action."
        ),
        alias_name: str = commands.parameter(description="Alias/Pseudonym"),
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            description="Tag a channel or include the ID",
        ),
        *,
        role: discord.Role = commands.parameter(
            converter=commands.RoleConverter,
            default=None,
            description="Tag a role or include the ID.",
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
        kwargs = {"alias_name": alias_name, "category": category, "channel": channel}
        if role:
            kwargs.update({"role": role})
        msg = await self.__alias_service.create_alias(**kwargs)
        return await state.end(success=msg)

    @commands.command(name="aroles", help="Administrator roles.")
    # @skip_help_discovery()
    async def list_administrator_roles_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel
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
        obj = target or ctx.guild
        object_dict = self.__discord_object_service.translate(obj=obj)
        is_at_home = at_home(source=ctx)
        pages = await self.__administrator_role_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="cap", help="Cap alias duration for mods.")
    # @skip_help_discovery()
    async def cap_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            description="Tag a channel or include its ID.",
        ),
        category: Category = commands.parameter(
            description="One of: `mute`, `ban`, `tmute`"
        ),
        *,
        hours: int = commands.parameter(default=24, description="# of hours"),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        channel_dict = self.__discord_object_service.translate(obj=channel)
        msg = await self.__cap_service.toggle_cap(
            category=category, channel_dict=channel_dict, hours=hours
        )
        return await state.end(success=msg)

    @commands.command(name="caps", help="List caps.")
    # @skip_help_discovery()
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel, None
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
        object_dict = self.__discord_object_service.translate(obj=obj)
        pages = await self.__cap_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="clear", help="Reset database.")
    async def clear_channel_access_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel, discord.Member
        ] = commands.parameter(
            converter=MultiConverter,
            description="Specify 'all', tag a channel/guild/member or include its ID",
        ),
        *,
        category: Category = commands.parameter(
            default="all",
            description="Specify one of: `alias`, `arole`, `all`, `ban`, `coord`, "
            "flag`, `mod`, `troom`, `tmute`, `stage`, `stream`, `vegan`, `vmute` or `vroom`.",
        ),
    ):
        default_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        object_dict = self.__discord_object_service.translate(obj=target)
        where_kwargs = object_dict.get("columns", None)
        # view = VerifyView(
        #     category=str(category),
        #     mention=object_dict.get("mention", "All"),
        #     author_snowflake=ctx.author.id,
        #     **where_kwargs,
        # )
        # embed = view.build_embed()
        # await ctx.reply(embed=embed, view=view)
        # await view.wait()
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        # msg = await ClearService.clear(
        #     category=category,
        #     default_kwargs=default_kwargs,
        #     object_dict=object_dict,
        #     where_kwargs=where_kwargs,
        #     target=target,
        #     view=view,
        # )
        return await state.end(success=msg)

    @commands.command(name="coord", help="Grant/revoke coords.")
    async def toggle_coordinator_text_command(
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
        channel_dict = self.__discord_object_service.translate(obj=channel)
        member_dict = self.__discord_object_service.translate(obj=member)
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(channel_dict.get("columns", None))
        # await PermissionService.has_equal_or_lower_role(
        #     target_member_snowflake=int(member_dict.get("id", None)),
        #     **updated_kwargs,
        # )
        msg = await self.__coordinator_service.toggle_coordinator(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            member_dict=member_dict,
        )
        return await state.end(success=msg)

    @commands.command(name="debug", help="Shows the last `n` number of logging.")
    async def debug_text_command(
        self,
        ctx,
        *,
        lines: int = commands.parameter(
            default=3, description="Specify the number of lines"
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
        if lines <= 0:
            return await state.end(warning="Lines must be greater than 0")
        try:
            with open(PATH_LOG, "r") as f:
                content = f.readlines()[-lines:]
                content = [line.split(" - ", 3)[-1] for line in content]
        except FileNotFoundError:
            return await state.end(warning="Log file not found")
        output = "".join(content)
        if len(output) > 1900:
            output = output[-1900:]
        return await state.end(success=f"```log\n{output}\n```")

    # @commands.command(name="pc", help="View permissions.")
    # @skip_help_discovery()
    # async def list_permissions_text_command(
    #     self,
    #     ctx: commands.Context,
    #     target: str | None = commands.parameter(
    #         default=None,
    #         description="Specify one of: `all`, channel ID/mention or server ID.",
    #     ),
    # ):
    #     state = StateService(ctx=ctx)
    #     default_kwargs = {
    #         "channel_snowflake": int(ctx.channel.id),
    #         "guild_snowflake": int(ctx.guild.id),
    #         "member_snowflake": int(ctx.author.id),
    #     }
    #     obj = target or int(ctx.channel.id)
    #     object_dict = self.__discord_object_service.translate(obj=obj)
    #     is_at_home = at_home(source=ctx)
    #     updated_kwargs = default_kwargs.copy()
    #     updated_kwargs.update(object_dict.get("columns", None))
    #     if target and str(target).lower() == "all":
    #         # await PermissionService.check(**updated_kwargs, lowest_role="Guild Owner")
    #         channel_objs = [
    #             channel_obj
    #             for guild in self.bot.guilds
    #             for channel_obj in guild.channels
    #         ]
    #     elif hasattr(object_dict.get("object", None), "channels"):
    #         channel_objs = object_dict.get("object", None).channels
    #     else:
    #         channel_objs = [object_dict.get("object", None)]
    #     pages = await PermissionService.build_pages(
    #         channel_objs=channel_objs,
    #         default_kwargs=default_kwargs,
    #         is_at_home=is_at_home,
    #     )
    #     if pages:
    #         return await state.end(warning=pages)
    #     else:
    #         return await state.end(
    #             success=f"{self.bot.user.display_name} has all permissions for `{target}`."
    #         )

    @commands.command(name="rmute", help="Room mute (except yourself).")
    async def room_mute_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel | None = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
        *,
        reason: str = commands.parameter(
            default="No reason provided.", description="Specify a reason."
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
        channel_dict = self.__discord_object_service.translate(obj=obj)
        pages = await self.__voice_mute_service.room_mute(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            reason=reason,
        )
        return await state.end(success=pages)

    @commands.command(name="rmv", help="VC move.")
    async def room_move_all_text_command(
        self,
        ctx: commands.Context,
        source_channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            description="Tag a channel or include its ID.",
        ),
        target_channel: discord.abc.GuildChannel = commands.parameter(
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
        failed, moved = [], []
        source_channel_dict = self.__discord_object_service.translate(
            obj=source_channel
        )
        target_channel_dict = self.__discord_object_service.translate(
            obj=target_channel
        )
        for member in source_channel_dict.get("object", None).members:
            try:
                await member.move_to(target_channel_dict.get("object", None))
                moved.append(member)
            except discord.Forbidden as e:
                failed.append(member)
                logger.warning(
                    f"Unable to move member "
                    f"{member.display_name} ({member.id}) from channel "
                    f"{source_channel_dict.get('name', None)} ({source_channel}) to channel "
                    f"{target_channel_dict.get('name', None)} ({target_channel}) in guild "
                    f"{ctx.guild.name} ({ctx.guild.id}). "
                    f"{str(e).capitalize()}"
                )
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} "
            f"Moved {source_channel_dict.get('mention', None)} to "
            f"{target_channel_dict.get('mention', None)}",
            color=discord.Color.green(),
        )
        if moved:
            embed.add_field(
                name=f"Successfully Moved (`{len(moved)}`)",
                value="\n".join(member.mention for member in moved),
                inline=False,
            )
        else:
            embed.add_field(name="Successfully Moved", value="None", inline=False)
        if failed:
            embed.add_field(
                name=f"Failed to Move ({len(failed)})",
                value="\n".join(member.mention for member in failed),
                inline=False,
            )
        embed.set_footer(
            text=f"Moved from {source_channel_dict.get('name', None)} "
            f"to {target_channel_dict.get('name', None)}"
        )
        return await state.end(success=embed)

    @commands.command(name="roleid", help="Get role by name.")
    # @skip_help_discovery()
    async def get_role_id_text_command(self, ctx: commands.Context, *, role_name: str):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            return await state.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await state.end(
                warning=f"No role named `{role_name}` found in this server."
            )

    @commands.command(name="smute", help="Server mute/server unmute.")
    async def toggle_server_mute_text_command(
        self,
        ctx: commands.Context,
        member: discord.Member = commands.parameter(
            converter=commands.MemberConverter,
            description="Tag a member or include their ID",
        ),
        *,
        reason: str = commands.parameter(
            default="No reason provided",
            description="Optional reason (required for 7 days or more)",
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
        member_dict = self.__discord_object_service.translate(obj=member)
        msg = await self.__server_mute_service.toggle_server_mute(
            default_kwargs=default_kwargs, member_dict=member_dict, reason=reason
        )
        return await state.end(success=msg)

    @commands.command(name="smutes", help="List mutes.")
    async def list_server_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel, discord.Member, None
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
        object_dict = self.__discord_object_service.translate(obj=obj)
        is_at_home = at_home(source=ctx)
        pages = await self.__server_mute_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="stages", help="List stages.")
    # @skip_help_discovery()
    async def list_stages_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel, discord.Member, None
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
        object_dict = self.__discord_object_service.translate(obj=obj)
        is_at_home = at_home(source=ctx)
        pages = await self.__stage_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return state.end(success=pages)

    @commands.command(
        name="temp", help="Toggle a temporary room and assign an owner.", hidden=True
    )
    # @skip_help_discovery()
    async def toggle_temp_room_text_command(
        self,
        ctx: commands.Context,
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
        channel_dict = self.__discord_object_service.translate(obj=channel)
        msg = await self.__temporary_room_service.toggle_temporary_room(
            channel_dict=channel_dict
        )
        return await state.end(success=msg)

    @commands.command(
        name="temps",
        help="List temporary rooms with matching command aliases.",
    )
    # @skip_help_discovery()
    async def list_temp_rooms_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel, None
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
        object_dict = self.__discord_object_service.translate(obj=obj)
        pages = await self.__temporary_room_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="stream", help="Setup streaming.")
    # @skip_help_discovery()
    async def modify_streaming_text_command(
        self,
        ctx: commands.Context,
        target_channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            description="Tag a channel or include its ID.",
        ),
        source_channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None
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
        # default_kwargs = {
        #     "channel_snowflake": int(ctx.channel.id),
        #     "guild_snowflake": int(ctx.guild.id),
        #     "member_snowflake": int(ctx.author.id),
        # }
        target_channel_dict = self.__discord_object_service.translate(obj=channel)
        source_channel_dict = self.__discord_object_service.translate(obj=channel)
        pages = await self.__stream_service.toggle_stream(
            default_kwargs=default_kwargs,
            source_channel_dict=source_channel_dict,
            target_channel_dict=target_channel_dict,
        )
        return await state.end(success=pages)

    @commands.command(name="streams", help="List streaming routes.")
    # @skip_help_discovery()
    async def list_streaming_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel
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
        object_dict = self.__discord_object_service.translate(obj=obj)
        pages = await self.__stream_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="vr", help="Start/stop video-only room.")
    # @skip_help_discovery()
    async def toggle_video_room_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include the ID",
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
        channel_dict = self.__discord_object_service.translate(obj=obj)
        msg = await self.__video_room_service.toggle_video_room(
            channel_dict=channel_dict
        )
        return await state.end(success=msg)

    @commands.command(
        name="vrs",
        help="List video rooms.",
    )
    # @skip_help_discovery()
    async def list_video_rooms_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Include `all`, channel or server ID.",
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
        object_dict = self.__discord_object_service.translate(obj=obj)
        is_at_home = at_home(source=ctx)
        pages = await self.__video_room_service.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(name="xalias", help="Delete alias.")
    async def delete_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: str = commands.parameter(description="Include an alias name"),
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
        msg = await self.__alias_service.delete_alias(
            alias_name=alias_name, default_kwargs=default_kwargs
        )
        return await state.end(success=msg)

    @commands.command(name="xrmute", help="Unmute all.")
    async def room_unmute_text_command(
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
        channel_dict = self.__discord_object_service.translate(obj=obj)
        pages = await self.__voice_mute_service.room_unmute(
            channel_dict=channel_dict, guild_snowflake=ctx.guild.id
        )
        return state.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(AdminTextCommands(bot))
