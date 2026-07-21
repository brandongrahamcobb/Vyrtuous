"""!/bin/python3
admin_text_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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

from vyrtuous.aliases import alias_service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.listing import list_overwrites, list_server_mutes
from vyrtuous.models.category import Category, CategoryObject
from vyrtuous.models.duration import Duration, DurationObject, DurationWrapper
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.models.scope import Scope, ScopeObject
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    clear_service,
    server_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.users import coordinator_service, moderator_service
from vyrtuous.view.cancel_confirm_view import VerifyView


class AdministratorTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Administrator"

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
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            description="Tag a channel or include the ID",
        ),
        role: discord.Role | None = commands.parameter(
            converter=commands.RoleConverter,
            default=None,
            description="Tag a role or include the ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if role is None:
            msg = await alias_service.create_alias(
                alias_name=alias_name,
                category=category.category,
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
            )
            return await tick.end(success=msg)
        else:
            msg = await alias_service.create_alias(
                alias_name=alias_name,
                category=category.category,
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                role_snowflake=role.id,
            )
            return await tick.end(success=msg)

    @commands.command(name="clear", help="Reset records.")
    async def clear_channel_access_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            int, discord.Guild, discord.abc.GuildChannel, discord.Member, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify 'all', a channel ID/mention, a member ID/mention or server ID.",
        ),
        *,
        category: CategoryObject = commands.parameter(
            converter=Category,
            default="all",
            description="Specify one of: `admin`, `alias`, `all`, `automute`, `ban`, `coord`, "
            "flag`, `mod`, `tmute`, `stream` or `vmute`.",
        ),
        scope: ScopeObject | None = commands.parameter(
            converter=Scope,
            default=None,
            description="Specify one of: `auto`, `click`, `command` or `server`.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        obj = target
        view = VerifyView(
            author_snowflake=ctx.author.id,
            category=str(category),
            obj=obj,
        )
        embed = view.build_embed()
        await tick.end(success=embed, view=view)
        await view.wait()
        tick = Tick(bot=self.__bot, ctx=ctx)
        if scope is None:
            mute_type = "click"
        else:
            mute_type = scope.scope
        msg = await clear_service.clear(
            author_snowflake=ctx.author.id,
            category=str(category),
            guild_snowflake=ctx.guild.id,
            message_snowflake=ctx.message.id,
            message_channel_snowflake=ctx.channel.id,
            obj=obj,
            target=mute_type,
            view=view,
        )
        return await tick.end(success=msg)

    @commands.command(name="coord", help="Grant/revoke coords.")
    async def toggle_coordinator_text_command(
        self,
        ctx: commands.Context,
        member: int | discord.Member = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID",
        ),
        channel: discord.abc.GuildChannel | None = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        if channel is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
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
        if isinstance(member, int):
            member_snowflake = member
        elif isinstance(member, discord.Member):
            member_snowflake = member.id
        await moderator_service.has_equal_or_lower_role(
            target_member_snowflake=member_snowflake,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
        )
        msg = await coordinator_service.toggle_coordinator(
            author_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            message_snowflake=ctx.message.id,
            message_channel_snowflake=ctx.message.channel.id,
        )
        return await tick.end(success=msg)

    @commands.command(name="ow", help="Overwrite stats.")
    async def list_overwrites_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a channel ID/mention or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target is None:
            if ctx.channel is None:
                return await tick.end(
                    warning=f"This command must target a valid channel."
                )
            obj = ctx.channel
        else:
            obj = target
        embed = list_overwrites.build_embed(obj=obj)
        return await tick.end(success=embed)

    @commands.command(name="rmute", help="Room mute (except yourself).")
    async def channel_mute_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel | None = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
        duration: DurationWrapper | None = commands.parameter(
            converter=Duration,
            default=None,
            description="Specify a duration m/h/d.",
        ),
        reason: str = commands.parameter(
            default="No reason provided.", description="Specify a reason."
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        if channel is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
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
        if duration is None:
            duration_obj = DurationObject(number=1, prefix="", sign=1, unit="h")
        else:
            duration_obj = duration.duration
        pages = await voice_mute_service.channel_mute(
            author_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            duration=duration_obj,
            guild_snowflake=guild_snowflake,
            reason=reason,
        )
        return await tick.end(success=pages)

    @commands.command(name="rmv", help="VC move.")
    async def channel_move_all_text_command(
        self,
        ctx: commands.Context,
        target_channel: (
            discord.VoiceChannel | discord.StageChannel
        ) = commands.parameter(
            converter=MultiConverter,
            description="Tag a channel or include its ID.",
        ),
        source_channel: (
            discord.VoiceChannel | discord.StageChannel | None
        ) = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if isinstance(
            target_channel,
            (discord.VoiceChannel, discord.StageChannel),
        ):
            target_guild_snowflake = target_channel.guild.id
            target_guild_name = target_channel.guild.name
            target_channel_snowflake = target_channel.id
            target_channel_name = target_channel.name
            target_channel_obj = target_channel
            target_channel_mention = target_channel.mention
        if source_channel is None:
            if ctx.channel is None or not isinstance(
                ctx.channel, (discord.VoiceChannel, discord.StageChannel)
            ):
                return await tick.end(
                    warning="This command must target a valid source channel."
                )
            source_guild_snowflake = ctx.channel.guild.id
            source_guild_name = ctx.channel.guild.name
            source_channel_snowflake = ctx.channel.id
            source_channel_name = ctx.channel.name
            source_channel_members = ctx.channel.members
            source_channel_mention = ctx.channel.mention
        elif isinstance(
            source_channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            source_guild_snowflake = source_channel.guild.id
            source_guild_name = source_channel.guild.name
            source_channel_snowflake = source_channel.id
            source_channel_name = source_channel.name
            source_channel_members = source_channel.members
            source_channel_mention = source_channel.mention
        failed, moved = [], []
        for member in source_channel_members:
            try:
                await member.move_to(target_channel_obj)
                moved.append(member)
            except discord.Forbidden as e:
                failed.append(member)
                self.__bot.logger.warning(
                    f"Unable to move member "
                    f"{member.display_name} ({member.id}) from channel "
                    f"{source_channel_name} ({source_channel_snowflake}) in guild "
                    f"{source_guild_name} ({source_guild_snowflake}) to channel "
                    f"{target_channel_name} ({target_channel_snowflake}) in guild "
                    f"{target_guild_name} ({target_guild_snowflake}) to channel "
                    f"{str(e).capitalize()}"
                )
        embed = discord.Embed(
            title=f"{emojis.get_random_emoji()} "
            f"Moved {source_channel_mention} to "
            f"{target_channel_mention}",
            color=discord.Color.green(),
        )
        if moved:
            embed.add_field(
                name=f"Successfully Moved (`{len(moved)}`)",
                value="\n".join(member.mention for member in moved),
                inline=False,
            )
        else:
            embed.add_field(name="Successfully Moved", value="None", inline=False)
        if failed:
            embed.add_field(
                name=f"Failed to Move ({len(failed)})",
                value="\n".join(member.mention for member in failed),
                inline=False,
            )
        embed.set_footer(
            text=f"Moved from {source_channel_name} to {target_channel_name}"
        )
        return await tick.end(success=embed)

    @commands.command(name="smute", help="Server mute/server unmute.")
    async def toggle_server_mute_text_command(
        self,
        ctx: commands.Context,
        member: Union[int, discord.Member] = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID",
        ),
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify a server ID.",
        ),
        reason: str = commands.parameter(
            default="No reason provided",
            description="Optional reason (required for 7 days or more)",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None:
            guild_snowflake = ctx.guild.id
        else:
            guild_snowflake = guild.id
        if isinstance(member, int):
            member_snowflake = member
        elif isinstance(member, discord.Member):
            member_snowflake = member.id
        msg = await server_mute_service.toggle_server_mute(
            author_snowflake=ctx.author.id,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            reason=reason,
        )
        return await tick.end(success=msg)

    @commands.command(name="smutes", help="List mutes.")
    async def list_server_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            int, discord.Guild, discord.abc.GuildChannel, discord.Member, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
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
        if target is None:
            if ctx.channel is None:
                return await tick.end(
                    warning=f"This command must target a valid channel."
                )
            obj = ctx.channel
        else:
            obj = target
        pages = await list_server_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @commands.command(name="xalias", help="Delete alias.")
    async def delete_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: str = commands.parameter(description="Include an alias name"),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a guild.")
        msg = await alias_service.delete_alias(
            alias_name=alias_name, guild_snowflake=ctx.guild.id
        )
        return await tick.end(success=msg)

    @commands.command(name="xrmute", help="Unmute all.")
    async def channel_unmute_text_command(
        self,
        ctx: commands.Context,
        channel: Union[discord.abc.GuildChannel, None] = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a guild.")
        if channel is None:
            if ctx.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
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
        pages = await voice_mute_service.channel_unmute(
            author_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            target="user",
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(AdministratorTextCommands(bot=bot))
