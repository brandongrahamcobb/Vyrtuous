"""!/bin/python3
sysadmin_text_commands.py A discord.py cog containing sysadmin commands for the Vyrtuous bot.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import discord
from discord.ext import commands

from vyrtuous.aliases import alias_service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.models.category import Category, CategoryObject
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.permissions import permission_service


class AliasManagementTextCommands(commands.Cog):

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    @commands.command(
        name="alias",
        help="Alias creation.",
    )
    async def create_alias_text_command(
        self,
        ctx: commands.Context,
        category: CategoryObject = commands.parameter(
            converter=Category,
            description="Specify a category for a `ban`, `flag`, `role`, `tmute`, or `vmute` action.",
        ),
        alias_name: str = commands.parameter(description="Alias/Pseudonym"),
        channel: discord.abc.GuildChannel | None = commands.parameter(
            converter=MultiConverter,
            description="Specify a channel ID/mention",
        ),
        role: discord.Role | None = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a role ID/mention.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if channel is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            channel_snowflake = ctx.channel.id
            guild_snowflake = ctx.guild.id
        elif isinstance(
            channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = channel.id
            guild_snowflake = channel.guild.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            guild_snowflake=guild_snowflake,
            channel_snowflake=channel_snowflake,
            requested=["command.alias.create"],
        )
        if role is None:
            msg = await alias_service.create_alias(
                alias_name=alias_name,
                category=category.category,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
            )
            return await tick.end(success=msg)
        else:
            if isinstance(role, discord.Role):
                role_snowflake = role.id
            else:
                return await tick.end(warning="This command must target a valid role.")
            msg = await alias_service.create_alias(
                alias_name=alias_name,
                category=category.category,
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                role_snowflake=role_snowflake,
            )
            return await tick.end(success=msg)

    @commands.command(name="xalias", help="Delete alias.")
    async def delete_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: str = commands.parameter(description="Include an alias name"),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state = bot.registry.get(PermissionState)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        else:
            guild_snowflake = ctx.guild.id
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be executed in a server channel."
            )
        else:
            channel_snowflake = ctx.channel.id
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            guild_snowflake=guild_snowflake,
            channel_snowflake=channel_snowflake,
            requested=["command.alias.delete"],
        )
        msg = await alias_service.delete_alias(
            alias_name=alias_name, guild_snowflake=ctx.guild.id
        )
        return await tick.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(AliasManagementTextCommands(bot=bot))
