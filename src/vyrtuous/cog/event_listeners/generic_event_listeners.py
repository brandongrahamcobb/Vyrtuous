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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.alias.alias_context import AliasContext
from vyrtuous.ban.ban_service import BanService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.cap.cap_service import CapService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_service import DurationService
from vyrtuous.moderator.moderator_service import ModeratorService
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.discord_object_service import DiscordObjectService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.logger import logger
from vyrtuous.utils.message_service import PaginatorService
from vyrtuous.utils.state_service import StateService


class GenericEventListeners(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self._ready_done = False
        self.__bot = bot
        self.__author_service = AuthorService()
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__duration_service = DurationService()
        self.__paginator_service = PaginatorService(bot=self.__bot)
        self.__moderator_service = ModeratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__stream_service = StreamService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            moderator_service=self.__moderator_service,
            paginator_service=self.__paginator_service,
        )
        self.__ban_service = BanService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__text_mute_service = TextMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
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
        self.__cap_service = CapService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
        )
        self.__discord_object_service = DiscordObjectService()

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if after.author.bot:
            return
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)
            else:
                await self.on_message(after)

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message):
        try:
            if (
                not message.guild
                or self.config["release_mode"]
                and message.author.id == self.bot.user.id
            ):
                return
            await self.__ban_service.is_banned_then_kick_and_reset_cooldown(
                channel=message.channel, member=message.author
            )
            await self.__text_mute_service.is_text_muted_then_mute_and_reset_cooldown(
                channel=message.channel, member=message.author
            )
            if not message.content.startswith(self.config["discord_command_prefix"]):
                return
            state = StateService(
                author_service=self.__author_service,
                bot=self.__bot,
                bug_service=self.__bug_service,
                message=message,
                developer_service=self.__developer_service,
                emoji=self.__emoji,
            )
            ctx = AliasContext(
                bot=self.__bot,
                cap_service=self.__cap_service,
                database_factory=self.__database_factory,
                dictionary_service=self.__dictionary_service,
                discord_object_service=self.__discord_object_service,
                duration_service=self.__duration_service,
                emoji=self.__emoji,
                message=message,
            )
            await ctx.setup()
            await self.__moderator_service.has_equal_or_lower_role(
                channel_snowflake=ctx.target_channel_snowflake,
                guild_snowflake=ctx.source_guild_snowflake,
                member_snowflake=ctx.source_member_snowflake,
                target_member_snowflake=ctx.target_member_snowflake,
            )
            await ctx.alias.service.enforce_or_undo(
                ctx=ctx, source=message, state=state
            )
        except:
            import traceback

            traceback.print_exc()

    @commands.Cog.listener()
    async def on_command_error(self, ctx, error):
        state = StateService(ctx=ctx)
        logger.error(str(error))
        if isinstance(error, commands.BadArgument):
            return await state.end(error=str(error))
        elif isinstance(error, commands.CheckFailure):
            return await state.end(error=str(error))
        elif isinstance(error, commands.MissingRequiredArgument):
            missing = error.param.name
            return await state.end(error=f"Missing required argument: `{missing}`")

    @commands.Cog.listener()
    async def on_app_command_error(self, interaction, error):
        state = StateService(interaction=interaction)
        logger.error(str(error))
        if isinstance(error, app_commands.CheckFailure):
            return await state.end(error=str(error))

    @commands.Cog.listener()
    async def on_ready(self):
        if getattr(self, "_ready_done", False):
            return
        self._ready_done = True
        method_names = [cmd.callback.__name__ for cmd in self.bot.commands]
        logger.info(method_names)


async def setup(bot: DiscordBot):
    await bot.add_cog(GenericEventListeners(bot))
