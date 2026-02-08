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

from vyrtuous.db.alias.alias_context import AliasContext
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.discord_object_service import DiscordObject
from vyrtuous.commands.discord_object_service import DiscordObjectNotFound
from vyrtuous.db.alias.alias_service import AliasService
from vyrtuous.commands.messaging.state_service import StateService
from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.infractions.ban.ban_service import BanService
from vyrtuous.db.infractions.tmute.text_mute_service import TextMuteService
from vyrtuous.utils.logger import logger


class GenericEventListeners(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self._ready_done = False

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
            await BanService.ban_overwrite(
                channel=message.channel, member=message.author
            )
            await TextMuteService.text_mute_overwrite(
                channel=message.channel, member=message.author
            )
            if not message.content.startswith(self.config["discord_command_prefix"]):
                return
            state = StateService(message=message)
            ctx = AliasContext(message=message)
            await ctx.setup()
            await PermissionService.has_equal_or_lower_role(
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
        elif isinstance(error, DiscordObjectNotFound):
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
