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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.moderator import NotModerator
from vyrtuous.inc.helpers import at_home
from vyrtuous.listing import (
    list_administrators,
    list_vegans,
)
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import (
    administrator_service,
    coordinator_service,
    developer_service,
    guild_owner_service,
    moderator_service,
    sysadmin_service,
    vegan_service,
)


class HiddenModeratorTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Moderator"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx: commands.Context) -> bool:
        if ctx.guild is None:
            raise commands.CheckFailure("This command must be used inside a server.")
        context = SnowflakeContext(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
        )
        for verify in (
            sysadmin_service.is_sysadmin_wrapper,
            developer_service.is_developer_wrapper,
            guild_owner_service.is_guild_owner_wrapper,
            administrator_service.is_administrator_wrapper,
            coordinator_service.is_coordinator_at_all_wrapper,
            moderator_service.is_moderator_at_all_wrapper,
        ):
            try:
                if await verify(context=context):
                    return True
            except commands.CheckFailure:
                continue
        raise NotModerator

    @commands.command(name="admins", help="Lists admins.")
    @skip_text_command_help_discovery()
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
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.guild
        is_at_home = at_home(source=ctx)
        pages = await list_administrators.build_pages(is_at_home=is_at_home, obj=obj)
        return await tick.end(success=pages)

    # @commands.command(name="mstage", help="Toggle stage mute/unmute.")
    # @skip_text_command_help_discovery()
    # async def stage_mute_text_command(
    #     self,
    #     ctx: commands.Context,
    #     member: discord.Member = commands.parameter(
    #         converter=commands.MemberConverter,
    #         description="Tag a member or include their ID",
    #     ),
    #     channel: discord.abc.GuildChannel = commands.parameter(
    #         converter=commands.VoiceChannelConverter,
    #         description="Tag a channel or include its ID.",
    #     ),
    # ):
    #     tick = Tick(bot=self.__bot, ctx=ctx)
    #     if ctx.guild is None:
    #         return await tick.end(warning="This command must be used in a server.")
    #     context = SnowflakeContext(
    #         channel_snowflake=ctx.channel.id,
    #         guild_snowflake=ctx.guild.id,
    #         member_snowflake=ctx.author.id,
    #     )
    #     obj = channel or ctx.channel
    #     await moderator_service.check_minimum_role(
    #         channel_snowflake=obj.id,
    #         guild_snowflake=ctx.guild.id,
    #         member_snowflake=ctx.author.id,
    #         lowest_role="Moderator",
    #     )
    #     msg = await automute_channel_service.toggle_stage_mute(
    #         channel=obj,
    #         context=context,
    #         member=member,
    #     )
    #     return await tick.end(success=msg)
    @commands.command(name="survey", help="Survey stage members.")
    @skip_text_command_help_discovery()
    async def stage_survey_text_command(
        self,
        ctx: commands.Context,
        channel: Union[discord.abc.GuildChannel, None] = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        obj = channel or ctx.channel
        pages = await moderator_service.survey(
            channel=obj,
        )
        return await tick.end(success=pages)

    @commands.command(name="vcow", help="Toggle vegan.")
    async def toggle_vegan_text_command(
        self,
        ctx: commands.Context,
        member: int | discord.Member = commands.parameter(
            converter=commands.MemberConverter,
            description="Tag a member or include their ID.",
        ),
        *,
        notes: str | None = commands.parameter(
            default=None,
            description="Include notes.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server")
        if not vegan_service.is_vegan(
            guild_snowflake=ctx.guild.id, member_snowflake=ctx.author.id
        ):
            return await tick.end(warning="Author is not a vegan.")
        if isinstance(member, discord.Member):
            embed = await vegan_service.toggle_vegan(
                guild_snowflake=ctx.guild.id, member_snowflake=member.id, notes=notes
            )
            return await tick.end(success=embed)

    @commands.command(name="vegans", help="List new vegans.")
    @skip_text_command_help_discovery()
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
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.guild
        is_at_home = at_home(source=ctx)
        pages = await list_vegans.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenModeratorTextCommands(bot=bot))
