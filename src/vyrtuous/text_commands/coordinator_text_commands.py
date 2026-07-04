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
from vyrtuous.db.coordinator import NotCoordinator
from vyrtuous.inc.helpers import at_home
from vyrtuous.listing import list_bans
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import ban_service
from vyrtuous.utils.channels import automute_channel_service
from vyrtuous.utils.users import (
    administrator_service,
    coordinator_service,
    developer_service,
    guild_owner_service,
    moderator_service,
    sysadmin_service,
)


class CoordinatorTextCommands(commands.Cog):
    PERMISSION_LEVEL = "Coordinator"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx):
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
        ):
            try:
                if await verify(context=context):
                    return True
            except commands.CheckFailure:
                continue
        raise NotCoordinator

    @commands.command(name="blacklist", help="Blacklist overwrite cleanup.")
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
    ):
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
        channel = channel or ctx.channel
        msg = await ban_service.toggle_blacklist(
            channel=channel,
            member_snowflake=member.id,
        )
        return await tick.end(success=msg)

    @commands.command(name="blacklists", help="List blacklists.")
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
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        pages = await list_bans.build_blacklist_pages(is_at_home=is_at_home, obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="mod", help="Grant/revoke mods.")
    async def toggle_moderator_text_command(
        self,
        ctx: commands.Context,
        member: discord.Member = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID",
        ),
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            description="Tag a channel or include its ID.",
        ),
    ):
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
        msg = await moderator_service.toggle_moderator(
            channel=channel,
            member_snowflake=member.id,
        )
        return await tick.end(success=msg)

    @commands.command(name="roles", help="List role members.")
    async def list_roles_text_command(
        self,
        ctx: commands.Context,
        role: discord.Role = commands.parameter(
            converter=commands.RoleConverter,
            default=None,
            description="Tag a role or include its ID.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        embeds = []
        members = [member for member in ctx.guild.members if role in member.roles]
        chunk_size = 12
        for index in range(0, len(members), chunk_size):
            chunk = members[index : index + chunk_size]
            description = "\n".join(
                f"{position + 1}. {member.mention} ({member.id})"
                for position, member in enumerate(chunk, start=index)
            )
            embed = discord.Embed(
                title=f"{role.name} Members",
                description=description or "No members found.",
                color=role.color if role.color.value else discord.Color.blurple(),
            )

            embeds.append(embed)
        if not embeds:
            embed = discord.Embed(
                title=f"{role.name} Members",
                description="No members found.",
                color=role.color if role.color.value else discord.Color.blurple(),
            )
            embeds.append(embed)
        return await tick.end(success=embeds)

    # @commands.command(name="stage", help="Start/stop stage")
    # @skip_text_command_help_discovery()
    # async def toggle_stage_text_command(
    #     self,
    #     ctx: commands.Context,
    #     channel: discord.abc.GuildChannel = commands.parameter(
    #         converter=commands.VoiceChannelConverter,
    #         default=None,
    #         description="Tag a channel or include its ID.",
    #     ),
    #     *,
    #     duration: str = commands.parameter(
    #         default="1h",
    #         description="Options: (+|-)duration(m|h|d) 0 - permanent / 24h - default",
    #     ),
    # ):
    #     tick = Tick(bot=self.__bot, ctx=ctx)
    #     if ctx.guild is None:
    #         return await tick.end(warning="This command must be executed in a server.")
    #     context = SnowflakeContext(
    #         channel_snowflake=ctx.channel.id,
    #         guild_snowflake=ctx.guild.id,
    #         member_snowflake=ctx.author.id,
    #     )
    #     resolved_channel = channel or ctx.channel
    #     await moderator_service.check_minimum_role(
    #         channel_snowflake=resolved_channel.id,
    #         guild_snowflake=ctx.guild.id,
    #         member_snowflake=ctx.author.id,
    #         lowest_role="Coordinator",
    #     )
    #     pages = await automute_channel_service.toggle_stage(
    #         channel=resolved_channel, context=context, duration_value=duration
    #     )
    #     return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorTextCommands(bot=bot))
