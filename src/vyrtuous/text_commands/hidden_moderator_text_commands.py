"""!/bin/python3
hidden_moderator_text_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.listing import list_administrators, list_vegans
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import moderator_service, vegan_service


class HiddenModeratorTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Moderator"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx: commands.Context) -> bool:
        if ctx.guild is None:
            raise commands.CheckFailure(
                "This command must be executed inside a server."
            )
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @commands.command(name="admins", help="Lists admins.")
    @skip_text_command_help_discovery()
    async def list_administrators_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            int, discord.Guild, discord.Role, discord.Member, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: member ID/mention or server ID.",
        ),
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify one a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None or not isinstance(guild, discord.Guild):
            guild_snowflake = ctx.guild.id
        else:
            guild_snowflake = guild.id
        if ctx.guild.id != guild_snowflake:
            await moderator_service.check_minimum_role(
                member_snowflake=ctx.author.id,
                lowest_role="Developer",
            )
        obj = target or ctx.guild
        pages = await list_administrators.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @commands.command(name="survey", help="Survey stage members.")
    @skip_text_command_help_discovery()
    async def survey_text_command(
        self,
        ctx: commands.Context,
        channel: Union[discord.abc.GuildChannel, None] = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if channel is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = ctx.channel.id
            guild_snowflake = ctx.guild.id
        else:
            channel_snowflake = channel.id
            guild_snowflake = channel.guild.id
        pages = await moderator_service.survey(
            channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake
        )
        return await tick.end(success=pages)

    @commands.command(name="vcow", help="Toggle vegan.")
    async def toggle_vegan_text_command(
        self,
        ctx: commands.Context,
        member: int | discord.Member = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID.",
        ),
        notes: str | None = commands.parameter(
            default=None,
            description="Include notes.",
        ),
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify one a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None or not isinstance(guild, discord.Guild):
            guild_snowflake = ctx.guild.id
        else:
            guild_snowflake = guild.id
        if not vegan_service.is_vegan(
            guild_snowflake=guild_snowflake, member_snowflake=ctx.author.id
        ):
            return await tick.end(warning="Author is not a vegan.")
        if isinstance(member, discord.Member):
            member_snowflake = member.id
        else:
            member_snowflake = member
        embed = await vegan_service.toggle_vegan(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            notes=notes,
        )
        return await tick.end(success=embed)

    @commands.command(name="vegans", help="List new vegans.")
    @skip_text_command_help_discovery()
    async def list_new_vegans_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            int, discord.Member, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a channel ID/mention, member ID/mention, or server ID.",
        ),
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify one a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None or not isinstance(guild, discord.Guild):
            guild_snowflake = ctx.guild.id
        else:
            guild_snowflake = guild.id
        if target is None:
            obj = ctx.guild
        else:
            obj = target
        pages = await list_vegans.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenModeratorTextCommands(bot=bot))
