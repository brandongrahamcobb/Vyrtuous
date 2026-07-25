"""!/bin/python3
generic_event_listeners.py A discord.py cog containing generic event listeners for the Vyrtuous bot.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from datetime import datetime, timezone

import discord
from discord.ext import commands

from vyrtuous.aliases.alias_context import AliasContext
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.alias import NotAlias
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.users import moderator_service, role_service


class GenericEventListeners(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self._ready_done = False
        self.__bot = bot

    @commands.Cog.listener()
    async def on_message_edit(self, before, after) -> None:
        if after.author.bot:
            return
        if before.content != after.content:
            ctx = await self.__bot.get_context(after)
            if ctx.command:
                await self.__bot.invoke(ctx)
            else:
                await self.on_message(after)

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message) -> discord.Message | None:
        if isinstance(message.channel, discord.VoiceChannel) or isinstance(
            message.channel, discord.StageChannel
        ):
            self.__bot.registry.get(MemberState).active[message.author.id] = (
                message.author.display_name,
                datetime.now(timezone.utc),
            )
        if (
            not message.guild
            or self.__bot.config["release_mode"]
            and self.__bot.user is not None
            and message.author.id == self.__bot.user.id
        ):
            return None
        if not isinstance(message.channel, discord.abc.GuildChannel):
            return None
        if not isinstance(message.author, discord.Member):
            return None
        await ban_service.is_banned_then_kick_and_reset_cooldown(
            channel=message.channel, member=message.author
        )
        await text_mute_service.is_text_muted_then_mute_and_reset_cooldown(
            channel=message.channel, member=message.author
        )
        if not message.content.startswith(self.__bot.config["discord_command_prefix"]):
            return None

        tick = Tick(bot=self.__bot, message=message)
        try:
            alias_ctx = AliasContext(
                content=message.content,
                guild_snowflake=message.guild.id,
            )
            try:
                if not await alias_ctx.setup():
                    return None
            except NotAlias:
                return None
            invincible_members_in_guild = self.__bot.registry.get(
                MemberState
            ).invincible.get(message.guild.id)
            if invincible_members_in_guild is not None:
                if alias_ctx.member_snowflake in invincible_members_in_guild:
                    return None
            await moderator_service.has_equal_or_lower_role(
                channel_snowflake=alias_ctx.channel_snowflake,
                guild_snowflake=message.guild.id,
                member_snowflake=message.author.id,
                target_member_snowflake=alias_ctx.member_snowflake,
            )
            match alias_ctx.category:
                case "ban":
                    embed = await ban_service.enforce_or_undo(
                        alias_ctx=alias_ctx, message=message
                    )
                    return await tick.end(success=embed)
                case "flag":
                    embed = await flag_service.enforce_or_undo(
                        alias_ctx=alias_ctx, message=message
                    )
                    return await tick.end(success=embed)
                case "role":
                    embed = await role_service.enforce_or_undo(
                        alias_ctx=alias_ctx, message=message
                    )
                    return await tick.end(success=embed)
                case "tmute":
                    embed = await text_mute_service.enforce_or_undo(
                        alias_ctx=alias_ctx, message=message
                    )
                    return await tick.end(success=embed)
                case "vmute":
                    if voice_mute_service.is_voice_muted(
                        guild_snowflake=alias_ctx.guild_snowflake,
                        member_snowflake=alias_ctx.member_snowflake,
                        targets=["server"],
                    ):
                        return None
                    embed = await voice_mute_service.enforce_or_undo(
                        alias_ctx=alias_ctx, message=message
                    )
                    return await tick.end(success=embed)
        except (
            commands.BadArgument,
            commands.CheckFailure,
            commands.MissingRequiredArgument,
        ) as e:
            if isinstance(e, commands.MissingRequiredArgument):
                missing = e.param.name
                return await tick.end(error=f"Missing required argument: `{missing}`")
            else:
                import traceback

                traceback.print_exc()
                return await tick.end(error=str(e))
        return None

    @commands.Cog.listener()
    async def on_command_error(self, ctx, error) -> discord.Message | None:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if isinstance(error, commands.BadArgument):
            return await tick.end(error=str(error))
        elif isinstance(error, commands.CheckFailure):
            return await tick.end(error=str(error))
        elif isinstance(error, commands.CommandInvokeError):
            return await tick.end(error=str(error))
        elif isinstance(error, commands.MissingRequiredArgument):
            missing = error.param.name
            return await tick.end(error=f"Missing required argument: `{missing}`")
        return None

    @commands.Cog.listener()
    async def on_ready(self) -> None:
        if getattr(self, "_ready_done", False):
            return
        self._ready_done = True
        method_names = [cmd.callback.__name__ for cmd in self.__bot.commands]
        self.__bot.logger.info(method_names)


async def setup(bot: DiscordBot):
    await bot.add_cog(GenericEventListeners(bot))
