"""!/bin/python3
coordinator_text_commands.py A discord.py cog containing coordinator commands for the Vyrtuous bot.

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
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import at_home
from vyrtuous.listing import list_bans
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import ban_service
from vyrtuous.utils.users import (
    moderator_service,
)


class HiddenCoordinatorTextCommands(commands.Cog):
    PERMISSION_LEVEL = "Coordinator"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx):
        if ctx.guild is None:
            raise commands.CheckFailure("This command must be used inside a server.")
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @commands.command(name="blacklist", help="Blacklist overwrite cleanup.")
    @skip_text_command_help_discovery()
    async def toggle_blacklist_text_command(
        self,
        ctx: commands.Context,
        member: discord.Member = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID",
        ),
        *,
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        await moderator_service.check_minimum_role(
            channel_snowflake=channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role="Coordinator",
        )
        await moderator_service.has_equal_or_lower_role(
            target_member_snowflake=int(member.id),
            member_snowflake=ctx.author.id,
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
        )
        target = channel or ctx.channel
        if not isinstance(
            target, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
        ):
            return await tick.end(warning="This command must target a valid channel.")
        msg = await ban_service.toggle_blacklist(
            channel=target,
            member_snowflake=member.id,
        )
        return await tick.end(success=msg)

    @commands.command(name="blacklists", help="List blacklists.")
    @skip_text_command_help_discovery()
    async def list_blacklists_text_command(
        self,
        ctx: commands.Context,
        *,
        target: (
            str | discord.Member | discord.Guild | discord.abc.GuildChannel
        ) = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Tag a channel, guild, member, all or include their ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        pages = await list_bans.build_blacklist_pages(is_at_home=is_at_home, obj=obj)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenCoordinatorTextCommands(bot=bot))
