"""!/bin/python3
generic_event_listeners.py A discord.py cog containing generic event listeners for the Vyrtuous bot.

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

import traceback

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.active_members import active_member_service
from vyrtuous.administrator.administrator_service import AdministratorService
from vyrtuous.alias.alias_context import AliasContext
from vyrtuous.ban.ban_service import BanService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.cap.cap_service import CapService
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_builder import DurationBuilder
from vyrtuous.moderator.moderator_service import ModeratorService
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.upload.upload_service import UploadService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.data_service import DataService
from vyrtuous.utils.default_context import DefaultContext
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.message_service import PaginatorService
from vyrtuous.utils.state_service import StateService
from vyrtuous.active_members.active_member_service import ActiveMemberService


class GenericEventListeners(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self._ready_done = False
        self.__bot = bot
        self.__author_service = AuthorService()
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__duration_builder = DurationBuilder()
        self.__paginator_service = PaginatorService(bot=self.__bot)
        self.__active_member_service = ActiveMemberService(
            bot=self.__bot, database_factory=self.__database_factory
        )
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
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            paginator_service=self.__paginator_service,
        )
        self.__ban_service = BanService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__text_mute_service = TextMuteService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__cap_service = CapService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
        )
        self.__data_service = DataService(
            database_factory=self.__database_factory,
            duration_builder=self.__duration_builder,
            moderator_service=self.__moderator_service,
        )
        self.__upload_service = UploadService(
            bot=self.__bot, database_factory=self.__database_factory
        )

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if after.author.bot:
            return
        if before.content != after.content:
            ctx = await self.__bot.get_context(after)
            if ctx.command:
                await self.__bot.invoke(ctx)
            else:
                await self.on_message(after)

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message):
        if isinstance(message.channel, discord.VoiceChannel) or isinstance(
            message.channel, discord.StageChannel
        ):
            await self.__active_member_service.update_active_member(
                guild_snowflake=message.guild.id,
                member_snowflake=message.author.id,
                name=message.author.display_name,
            )
        if (
            not message.guild
            or self.__bot.config["release_mode"]
            and message.author.id == self.__bot.user.id
        ):
            return
        await self.__ban_service.is_banned_then_kick_and_reset_cooldown(
            channel=message.channel, member=message.author
        )
        await self.__text_mute_service.is_text_muted_then_mute_and_reset_cooldown(
            channel=message.channel, member=message.author
        )
        if not message.content.startswith(self.__bot.config["discord_command_prefix"]):
            return
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            message=message,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            upload_service=self.__upload_service,
        )
        try:
            d_ctx = DefaultContext(message=message)
            ctx = AliasContext(
                active_member_service=self.__active_member_service,
                bot=self.__bot,
                cap_service=self.__cap_service,
                content=message.content,
                database_factory=self.__database_factory,
                default_ctx=d_ctx,
                dictionary_service=self.__dictionary_service,
                emoji=self.__emoji,
                moderator_service=self.__moderator_service,
            )
            if not await ctx.setup():
                return
            await self.__moderator_service.has_equal_or_lower_role(
                **d_ctx.to_dict(),
                target_member_snowflake=ctx.member_snowflake,
            )
            service = ctx.alias.service(
                active_member_service=self.__active_member_service,
                bot=self.__bot,
                database_factory=self.__database_factory,
                data_service=self.__data_service,
                dictionary_service=self.__dictionary_service,
                duration_builder=self.__duration_builder,
                emoji=self.__emoji,
                stream_service=self.__stream_service,
            )
            await service.enforce_or_undo(
                ctx=ctx, default_ctx=d_ctx, source=message, state=state
            )
        except (
            commands.BadArgument,
            commands.CheckFailure,
            commands.MissingRequiredArgument,
        ) as e:
            if isinstance(e, commands.MissingRequiredArgument):
                missing = e.param.name
                return await state.end(error=f"Missing required argument: `{missing}`")
            else:
                return await state.end(error=str(e))

    @commands.Cog.listener()
    async def on_command_error(self, ctx, error):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            upload_service=self.__upload_service,
        )
        if isinstance(error, commands.BadArgument):
            return await state.end(error=str(error))
        elif isinstance(error, commands.CheckFailure):
            return await state.end(error=str(error))
        elif isinstance(error, commands.CommandInvokeError):
            return await state.end(error=str(error))
        elif isinstance(error, commands.MissingRequiredArgument):
            missing = error.param.name
            return await state.end(error=f"Missing required argument: `{missing}`")

    @commands.Cog.listener()
    async def on_app_command_error(self, interaction, error):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            interaction=interaction,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        if isinstance(error, app_commands.CheckFailure):
            return await state.end(error=str(error))

    @commands.Cog.listener()
    async def on_ready(self):
        if getattr(self, "_ready_done", False):
            return
        self._ready_done = True
        method_names = [cmd.callback.__name__ for cmd in self.__bot.commands]
        self.__bot.logger.info(method_names)


async def setup(bot: DiscordBot):
    await bot.add_cog(GenericEventListeners(bot))
