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

from vyrtuous.administrator.administrator_service import AdministratorRoleService
from vyrtuous.ban.ban_service import BanService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.cog.text_commands.help_text_command import (
    skip_text_command_help_discovery,
)
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_service import DurationService
from vyrtuous.flag.flag_service import FlagService
from vyrtuous.moderator.moderator_service import ModeratorService
from vyrtuous.owner.guild_owner import GuildOwner
from vyrtuous.owner.guild_owner_service import GuildOwnerService, NotGuildOwner
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.data_service import DataService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.discord_object_service import DiscordObjectService, MultiConverter
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.hero_service import HeroService
from vyrtuous.utils.message_service import PaginatorService
from vyrtuous.utils.state_service import StateService
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService
from vyrtuous.utils.home import at_home


class GuildOwnerTextCommands(commands.Cog):
    ROLE = GuildOwner

    def __init__(self, *, bot: DiscordBot | None = None):
        self.__bot = bot
        self.__author_service = AuthorService()
        self.__emoji = Emojis()
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__duration_service = DurationService()
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
        self.__moderator_service = ModeratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            emoji=self.__emoji,
        )
        self.__data_service = DataService(
            duration_service=self.__duration_service,
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
        self.__administrator_role_service = AdministratorRoleService(
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
        self.__voice_mute_service = VoiceMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            stream_service=self.__stream_service,
        )
        self.__ban_service = BanService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
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
            duration_service=self.__duration_service,
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

    async def cog_check(self, ctx) -> Coroutine[Any, Any, bool]:
        async def predicate(
            source: Union[commands.Context, discord.Interaction, discord.Message],
        ):
            for verify in (
                self.__sysadmin_service.is_sysadmin_wrapper,
                self.__developer_service.is_developer_wrapper,
                self.__guild_owner_service.is_guild_owner_wrapper,
            ):
                try:
                    if await verify(source):
                        return True
                except commands.CheckFailure:
                    continue
            raise NotGuildOwner

        predicate._permission_level = "Guild Owner"
        return await predicate(ctx)

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
        )
        default_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        role_dict = self.__discord_object_service.to_dict(obj=role)
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(role_dict.get("columns", None))
        pages = await self.__administrator_role_service.toggle_administrator_role(
            role_dict=role_dict, updated_kwargs=updated_kwargs
        )
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
        )
        member_dict = self.__discord_object_service.to_dict(obj=member)
        where_kwargs = member_dict.get("columns", None)
        enabled = self.__hero_service.toggle_enabled()
        if enabled:
            self.__hero_service + (
                where_kwargs.get("guild_snowflake", None),
                where_kwargs.get("member_snowflake", None),
            )
            await self.__hero_service.unrestrict(**where_kwargs)
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member_dict.get('mention', None)}."
            )
        else:
            self.__hero_service - (
                where_kwargs.get("guild_snowflake", None),
                where_kwargs.get("member_snowflake", None),
            )
            msg = f"Invincibility has been disabled for {member_dict.get('mention', None)}"
        return await state.end(success=msg)

    @commands.command(name="devs", help="List devs.")
    async def list_developers_text_command(
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
        )
        obj = target or "all"
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__developer_service.build_pages(object_dict=object_dict)
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
        )
        default_kwargs = {
            "channel_snowflake": ctx.channel.id,
            "guild_snowflake": ctx.guild.id,
            "member_snowflake": ctx.author.id,
        }
        is_at_home = at_home(source=ctx)
        obj = target or "all"
        object_dict = self.__discord_object_service.to_dict(obj=obj)
        pages = await self.__hero_service.build_pages(
            is_at_home=is_at_home,
            default_kwargs=default_kwargs,
            object_dict=object_dict,
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
