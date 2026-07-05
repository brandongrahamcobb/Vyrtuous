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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.inc.helpers import at_home
from vyrtuous.listing import list_heroes
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import (
    hero_service,
    moderator_service,
)


class HiddenGuildOwnerTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Guild Owner"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx) -> bool:
        if ctx.guild is None:
            raise commands.CheckFailure("This command must be used inside a server.")
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @commands.command(name="hero", help="Grant/revoke invincibility.")
    @skip_text_command_help_discovery()
    async def invincibility_text_command(
        self,
        ctx: commands.Context,
        member: discord.Member = commands.parameter(
            converter=commands.MemberConverter,
            description="Tag a member or include their ID",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if member.id in (self.__bot.registry.get(MemberState).invincible[ctx.guild.id]):
            await hero_service.add_invincible_member(member.guild.id, member.id)
            await hero_service.unrestrict(
                guild_snowflake=member.guild.id, member_snowflake=member.id
            )
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member.mention}."
            )
        else:
            await hero_service.remove_invincible_member(member.guild.id, member.id)
            self.__bot.registry.get(MemberState).invincible[ctx.guild.id].remove(
                member.id
            )
            msg = f"Invincibility has been disabled for {member.mention}"
        return await tick.end(success=msg)

    @commands.command(name="heroes", help="List heroes.")
    @skip_text_command_help_discovery()
    async def list_heroes_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[str, discord.Member] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="'all', or user mention/ID",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        is_at_home = at_home(source=ctx)
        pages = await list_heroes.build_pages(
            is_at_home=is_at_home,
            obj=target,
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenGuildOwnerTextCommands(bot=bot))
