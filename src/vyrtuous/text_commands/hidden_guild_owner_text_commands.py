"""!/bin/python3
hidden_guild_owner_text_commands.py A discord.py cog containing guild owner commands for the Vyrtuous bot.

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
from vyrtuous.cache.registry import MemberState
from vyrtuous.listing import list_heroes
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import hero_service, moderator_service


class HiddenGuildOwnerTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Guild Owner"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx) -> bool:
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

    @commands.command(name="hero", help="Grant/revoke invincibility.")
    @skip_text_command_help_discovery()
    async def toggle_invincibility_text_command(
        self,
        ctx: commands.Context,
        member: int | discord.Member = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID.",
        ),
        target: Union[str, discord.Guild, None] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        singular = True
        if target is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guilds = [ctx.guild]
        else:
            if isinstance(target, discord.Guild):
                guilds = [target]
            elif target == "all":
                guilds = bot.guilds
                singular = False
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if isinstance(member, int):
            member_snowflake = member
            member_name = "Unknown"
            simplified_member = self.__bot.registry.get(MemberState).active.get(
                member_snowflake
            )
            if simplified_member:
                member_name = simplified_member[0]
        else:
            member_snowflake = member.id
            member_name = member.mention
        if not singular:
            if any(
                member_snowflake in member_set
                for member_set in self.__bot.registry.get(
                    MemberState
                ).invincible.values()
            ):
                enabled = False
            else:
                enabled = True
        else:
            guild_snowflake = guilds[0].id
            if member_snowflake in self.__bot.registry.get(MemberState).invincible.get(
                guild_snowflake, set()
            ):
                enabled = False
            else:
                enabled = True
        if enabled:
            header = f"All moderation events have been forgiven and invincibility has been enabled for {member_name} in "
        else:
            header = f"Invincibility has been disabled for {member_name} in "
        guild_names = []
        for g in guilds:
            if enabled:
                await hero_service.add_invincible_member(g.id, member_snowflake)
                await hero_service.unrestrict(
                    guild_snowflake=g.id, member_snowflake=member_snowflake
                )
                guild_names.append(g.name)
            else:
                header = f"Invincibility has been disabled for {member_name} in "
                await hero_service.remove_invincible_member(g.id, member_snowflake)
                guild_names.append(g.name)
        msg = header + ", ".join(guild_names) + "."
        return await tick.end(success=msg)

    @commands.command(name="heroes", help="List heroes.")
    @skip_text_command_help_discovery()
    async def list_heroes_text_command(
        self,
        ctx: commands.Context,
        target: Union[int, discord.Member, discord.Guild, None] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a member mention/ID or server ID.",
        ),
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None:
            guild_snowflake = ctx.guild.id
        else:
            guild_snowflake = guild.id
        obj = target or ctx.guild
        pages = await list_heroes.build_pages(
            guild_snowflake=guild_snowflake,
            obj=obj,
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenGuildOwnerTextCommands(bot=bot))
