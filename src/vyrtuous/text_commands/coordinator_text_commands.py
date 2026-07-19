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
from vyrtuous.models.duration import Duration, DurationObject, DurationWrapper
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.utils.channels import automute_channel_service
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import moderator_service


class CoordinatorTextCommands(commands.Cog):
    PERMISSION_LEVEL = "Coordinator"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx):
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
    ) -> discord.Message:
        try:
            tick = Tick(bot=self.__bot, ctx=ctx)
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must be executed in a server."
                )
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
                author_snowflake=ctx.author.id,
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                member_snowflake=member.id,
                message_snowflake=ctx.message.id,
                message_channel_snowflake=ctx.message.channel.id,
            )
        except:
            import traceback

            traceback.print_exc()
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
    ) -> discord.Message:
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

    @commands.command(name="automute", help="Start/stop automute")
    async def toggle_automute_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
        duration: DurationWrapper | None = commands.parameter(
            converter=Duration,
            default=None,
            description="m/h/d",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        if duration is None:
            duration_obj = DurationObject(number=1, prefix="+", sign=1, unit="h")
        else:
            duration_obj = duration.duration
        resolved_channel = channel or ctx.channel
        pages = await automute_channel_service.toggle_automute(
            author_snowflake=ctx.author.id,
            channel_snowflake=resolved_channel.id,
            guild_snowflake=ctx.guild.id,
            duration=duration_obj,
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorTextCommands(bot=bot))
