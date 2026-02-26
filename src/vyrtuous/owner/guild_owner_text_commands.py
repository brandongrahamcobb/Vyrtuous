"""!/bin/python3
guild_owner_text_commands.py A discord.py cog containing guild owner commands for the Vyrtuous bot.

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

from typing import Any, Coroutine, Literal, Optional, Union

import discord
from discord.ext import commands

from vyrtuous.administrator.administrator_service import (
    AdministratorRoleService,
    AdministratorService,
)
from vyrtuous.ban.ban_service import BanService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_builder import DurationBuilder
from vyrtuous.flag.flag_service import FlagService
from vyrtuous.moderator.help_text_command import skip_text_command_help_discovery
from vyrtuous.moderator.moderator_service import ModeratorService
from vyrtuous.owner.guild_owner import GuildOwner
from vyrtuous.owner.guild_owner_service import GuildOwnerService, NotGuildOwner
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.upload.upload_service import UploadService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.data_service import DataService
from vyrtuous.utils.default_context import DefaultContext
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.discord_object_service import DiscordObjectService, MultiConverter
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.hero_service import HeroService
from vyrtuous.utils.home import at_home
from vyrtuous.utils.message_service import PaginatorService
from vyrtuous.utils.state_service import StateService
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService


class GuildOwnerTextCommands(commands.Cog):
    ROLE = GuildOwner

    def __init__(self, *, bot: DiscordBot | None = None):
        self.__bot = bot
        self.__author_service = AuthorService()
        self.__emoji = Emojis()
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__duration_builder = DurationBuilder()
        self.__bug_service = BugService(
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
        self.__coordinator_service = CoordinatorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__developer_service = DeveloperService(
            bot=self.__bot,
            bug_service=self.__bug_service,
            database_factory=self.__database_factory,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
        )
        self.__guild_owner_service = GuildOwnerService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__sysadmin_service = SysadminService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__moderator_service = ModeratorService(
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
        self.__data_service = DataService(
            database_factory=self.__database_factory,
            duration_builder=self.__duration_builder,
            moderator_service=self.__moderator_service,
        )
        self.__guild_owner_service = GuildOwnerService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__paginator_service = PaginatorService(bot=self.__bot)
        self.__stream_service = StreamService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            paginator_service=self.__paginator_service,
        )
        self.__discord_object_service = DiscordObjectService()
        self.__voice_mute_service = VoiceMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            stream_service=self.__stream_service,
        )
        self.__ban_service = BanService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
        )
        self.__flag_service = FlagService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__text_mute_service = TextMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
        )
        self.__hero_service = HeroService(
            ban_service=self.__ban_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            flag_service=self.__flag_service,
            text_mute_service=self.__text_mute_service,
            voice_mute_service=self.__voice_mute_service,
        )
        self.__upload_service = UploadService(
            bot=self.__bot, database_factory=self.__database_factory
        )
        self.__administrator_role_service = AdministratorRoleService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )

    async def cog_check(self, ctx) -> Coroutine[Any, Any, bool]:
        context = DefaultContext(ctx=ctx)
        for verify in (
            self.__sysadmin_service.is_sysadmin_wrapper,
            self.__developer_service.is_developer_wrapper,
            self.__guild_owner_service.is_guild_owner_wrapper,
        ):
            try:
                if await verify(context=context):
                    return True
            except commands.CheckFailure:
                continue
        raise NotGuildOwner

    cog_check._permission_level = "Guild Owner"

    @commands.command(name="admin", help="Toggle administrator role.")
    async def toggle_administrator_by_role_text_command(
        self,
        ctx: commands.Context,
        role: discord.Role = commands.parameter(
            converter=commands.RoleConverter, description="Tag a role or its ID"
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
        try:
            pages = await self.__administrator_role_service.toggle_administrator_role(
                role=role,
            )
        except:
            import traceback

            traceback.print_exc()
        return await state.end(success=pages)

    @commands.command(name="hero", help="Grant/revoke invincibility.")
    @skip_text_command_help_discovery()
    async def invincibility_text_command(
        self,
        ctx: commands.Context,
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
        enabled = self.__hero_service.toggle_enabled()
        if enabled:
            self.__hero_service + (member.guild.id, member.id)
            await self.__hero_service.unrestrict(
                guild_snowflake=member.guild.id, member_snowflake=member.id
            )
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member.mention}."
            )
        else:
            self.__hero_service - (member.guild.id, member.id)
            msg = f"Invincibility has been disabled for {member.mention}"
        return await state.end(success=msg)

    @commands.command(name="devs", help="List devs.")
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[str, discord.Member, None] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="'all', or user mention/ID",
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
        if not target:
            obj = None
        else:
            obj = target or ctx.guild
        try:
            pages = await self.__developer_service.build_pages(obj=obj)
        except:
            import traceback

            traceback.print_exc()
        return await state.end(success=pages)

    @commands.command(name="heroes", help="List heroes.")
    async def list_heroes_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[str, discord.Member] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="'all', or user mention/ID",
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
        is_at_home = at_home(source=ctx)
        pages = await self.__hero_service.build_pages(
            is_at_home=is_at_home,
            obj=target,
        )
        return await state.end(success=pages)

    @commands.command(name="sync", help="Sync app commands.")
    async def sync_text_command(
        self,
        ctx: commands.Context,
        spec: Optional[Literal["~", "*", "^"]] = None,
        *,
        guilds: commands.Greedy[discord.Object] = None,
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
                return await state.end(success=msg)
            except Exception as e:
                return await state.end(warning=str(e).capitalize())
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        return await state.end(success=f"Synced the tree to {ret}/{len(guilds)}.")


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildOwnerTextCommands(bot=bot))
