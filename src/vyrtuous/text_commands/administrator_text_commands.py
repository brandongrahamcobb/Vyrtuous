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
from vyrtuous.db.administrator import NotAdministrator
from vyrtuous.inc.helpers import PATH_LOG, at_home
from vyrtuous.listing import (
    list_administrator_roles,
    list_caps,
    list_permissions,
    list_server_mutes,
    list_streams,
    list_video_channels,
)
from vyrtuous.models.category import Category
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    cap_service,
    server_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.channels import automute_channel_service, video_channel_service
from vyrtuous.utils.tracking import stream_service
from vyrtuous.utils.users import (
    administrator_service,
    coordinator_service,
    developer_service,
    guild_owner_service,
    moderator_service,
    sysadmin_service,
)


class AdminTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Administrator"

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
        ):
            try:
                if await verify(context=context):
                    return True
            except commands.CheckFailure:
                continue
        raise NotAdministrator

    @commands.command(
        name="alias",
        help="Alias creation.",
    )
    async def create_alias_text_command(
        self,
        ctx: commands.Context,
        category: Category = commands.parameter(
            description="Specify a category for a `ban`, `flag`, `role`, `tmute`, `vegan` or `vmute` action."
        ),
        alias_name: str = commands.parameter(description="Alias/Pseudonym"),
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            description="Tag a channel or include the ID",
        ),
        *,
        role: discord.Role = commands.parameter(
            converter=commands.RoleConverter,
            default=None,
            description="Tag a role or include the ID.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        kwargs = {"alias_name": alias_name, "category": category, "channel": channel}
        if role:
            kwargs.update({"role": role})
        msg = await alias_service.create_alias(**kwargs)
        return await tick.end(success=msg)

    @commands.command(name="aroles", help="Administrator roles.")
    @skip_text_command_help_discovery()
    async def list_administrator_roles_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.guild
        is_at_home = at_home(source=ctx)
        pages = await list_administrator_roles.build_pages(
            obj=obj, is_at_home=is_at_home
        )
        return await tick.end(success=pages)

    @commands.command(name="cap", help="Cap alias duration for mods.")
    @skip_text_command_help_discovery()
    async def cap_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            description="Tag a channel or include its ID.",
        ),
        category: Category = commands.parameter(
            description="One of: `mute`, `ban`, `tmute`"
        ),
        *,
        hours: int = commands.parameter(default=24, description="# of hours"),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        msg = await cap_service.toggle_cap(
            category=category, channel=channel, hours=hours
        )
        return await tick.end(success=msg)

    @commands.command(name="caps", help="List caps.")
    @skip_text_command_help_discovery()
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention or server ID.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        pages = await list_caps.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    # @commands.command(name="clear", help="Reset database.")
    # @skip_text_command_help_discovery()
    # async def clear_channel_access_text_command(
    #     self,
    #     ctx: commands.Context,
    #     target: Union[
    #         str, discord.Guild, discord.abc.GuildChannel, discord.Member
    #     ] = commands.parameter(
    #         converter=MultiConverter,
    #         default=None,
    #         description="Specify 'all', tag a channel/guild/member or include its ID",
    #     ),
    #     *,
    #     category: Category = commands.parameter(
    #         default="all",
    #         description="Specify one of: `alias`, `arole`, `all`, `ban`, `coord`, "
    #         "flag`, `mod`, `troom`, `tmute`, `stage`, `stream`, `vegan`, `vmute` or `vroom`.",
    #     ),
    # ):
    #     tick = Tick(bot=self.__bot, ctx=ctx)
    #     context = SnowflakeContext(
    #         channel_snowflake=ctx.channel.id,
    #         guild_snowflake=ctx.guild.id,
    #         member_snowflake=ctx.author.id,
    #     )
    #     if target == "all":
    #         obj = None
    #     else:
    #         obj = target
    #     view = VerifyView(
    #         category=str(category),
    #         obj=obj,
    #         author_snowflake=ctx.author.id,
    #     )
    #     embed = view.build_embed()
    #     await tick.end(success=embed, view=view)
    #     await view.wait()
    #     tick = Tick(bot=self.__bot, ctx=ctx)
    #     msg = await self.__clear_service.clear(
    #         category=category,
    #         default_ctx=context,
    #         obj=obj,
    #         target=target,
    #         view=view,
    #         source=ctx,
    #     )
    #     return await tick.end(success=msg)

    @commands.command(name="coord", help="Grant/revoke coords.")
    async def toggle_coordinator_text_command(
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
        await moderator_service.has_equal_or_lower_role(
            target_member_snowflake=int(member.id),
            member_snowflake=ctx.author.id,
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
        )
        msg = await coordinator_service.toggle_coordinator(
            channel=channel,
            member_snowflake=member.id,
        )
        return await tick.end(success=msg)

    @commands.command(name="debug", help="Shows the last `n` number of logging.")
    async def debug_text_command(
        self,
        ctx,
        *,
        lines: int = commands.parameter(
            default=3, description="Specify the number of lines"
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if lines <= 0:
            return await tick.end(warning="Lines must be greater than 0")
        try:
            with open(PATH_LOG, "r") as f:
                content = f.readlines()[-lines:]
                content = [line.split(" - ", 3)[-1] for line in content]
        except FileNotFoundError:
            return await tick.end(warning="Log file not found")
        output = "".join(content)
        if len(output) > 1900:
            output = output[-1900:]
        return await tick.end(success=f"```log\n{output}\n```")

    @commands.command(name="ow", help="Overwrite stats.")
    async def list_overwrites_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify the ID or mention.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        member_count, role_count, total_count = 0, 0, 0
        if isinstance(target, str):
            for guild in self.__bot.guilds:
                for channel in guild.channels:
                    for target_obj, overwrite in channel.overwrites.items():
                        if any(v is not None for v in overwrite):
                            total_count += 1
                            if isinstance(target_obj, discord.Member):
                                member_count += 1
                            else:
                                role_count += 1
            embed = discord.Embed(
                title="All Servers Overwrites Summary", color=discord.Color.blue()
            )
            embed.add_field(name="Role overwrites", value=str(role_count), inline=False)
            embed.add_field(
                name="Member overwrites", value=str(member_count), inline=False
            )
            embed.add_field(
                name="Total overwrites", value=str(total_count), inline=False
            )
            return await tick.end(success=embed)
        elif isinstance(target, discord.Guild):
            for channel in target.channels:
                for target_obj, overwrite in channel.overwrites.items():
                    if any(v is not None for v in overwrite):
                        total_count += 1
                        if isinstance(target_obj, discord.Member):
                            member_count += 1
                        else:
                            role_count += 1
            embed = discord.Embed(
                title=f"Server {target.name} Overwrites Summary",
                color=discord.Color.blue(),
            )
            embed.add_field(name="Role overwrites", value=str(role_count), inline=False)
            embed.add_field(
                name="Member overwrites", value=str(member_count), inline=False
            )
            embed.add_field(
                name="Total overwrites", value=str(total_count), inline=False
            )
            return await tick.end(success=embed)
        else:
            obj = target or ctx.channel
            if isinstance(obj, (discord.TextChannel, discord.VoiceChannel)):
                for target_obj, overwrite in obj.overwrites.items():
                    if any(v is not None for v in overwrite):
                        total_count += 1
                        if isinstance(target_obj, discord.Member):
                            member_count += 1
                        else:
                            role_count += 1
                embed = discord.Embed(
                    title="Channel Overwrites", color=discord.Color.blue()
                )
                embed.add_field(name="Channel", value=obj.name, inline=False)
                embed.add_field(
                    name="Role overwrites", value=str(role_count), inline=False
                )
                embed.add_field(
                    name="Member overwrites", value=str(member_count), inline=False
                )
                embed.add_field(
                    name="Total overwrites", value=str(total_count), inline=False
                )
                return await tick.end(success=embed)

    @commands.command(name="pc", help="View permissions.")
    @skip_text_command_help_discovery()
    async def list_permissions_text_command(
        self,
        ctx: commands.Context,
        target: Union[str, discord.Guild, None] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: `all`, channel ID/mention or server ID.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        context = SnowflakeContext(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
        )
        if target == "all":
            obj = None
        else:
            obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        pages = await list_permissions.build_pages(
            obj=obj,
            context=context,
            is_at_home=is_at_home,
        )
        return await tick.end(success=pages)

    @commands.command(name="rmute", help="Room mute (except yourself).")
    async def channel_mute_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel | None = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include its ID.",
        ),
        *,
        reason: str = commands.parameter(
            default="No reason provided.", description="Specify a reason."
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        obj = channel or ctx.channel
        pages = await voice_mute_service.channel_mute(
            author=ctx.author,
            channel=obj,
            reason=reason,
        )
        return await tick.end(success=pages)

    @commands.command(name="rmv", help="VC move.")
    async def channel_move_all_text_command(
        self,
        ctx: commands.Context,
        source_channel: (
            discord.VoiceChannel | discord.StageChannel
        ) = commands.parameter(
            converter=MultiConverter,
            description="Tag a channel or include its ID.",
        ),
        target_channel: (
            discord.VoiceChannel | discord.StageChannel
        ) = commands.parameter(
            converter=MultiConverter,
            description="Tag a channel or include its ID.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used within a server.")
        failed, moved = [], []
        for member in source_channel.members:
            try:
                await member.move_to(target_channel)
                moved.append(member)
            except discord.Forbidden as e:
                failed.append(member)
                self.__bot.logger.warning(
                    f"Unable to move member "
                    f"{member.display_name} ({member.id}) from channel "
                    f"{source_channel.name} ({source_channel.id}) to channel "
                    f"{target_channel.name} ({target_channel.id}) in guild "
                    f"{ctx.guild.name} ({ctx.guild.id}). "
                    f"{str(e).capitalize()}"
                )
        embed = discord.Embed(
            title=f"{emojis.get_random_emoji()} "
            f"Moved {source_channel.mention} to "
            f"{target_channel.mention}",
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
            text=f"Moved from {source_channel.name} to {target_channel.name}"
        )
        return await tick.end(success=embed)

    @commands.command(name="roleid", help="Get role by name.")
    @skip_text_command_help_discovery()
    async def get_role_id_text_command(self, ctx: commands.Context, *, role_name: str):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used within a server.")
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            return await tick.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await tick.end(
                warning=f"No role named `{role_name}` found in this server."
            )

    @commands.command(name="smute", help="Server mute/server unmute.")
    async def toggle_server_mute_text_command(
        self,
        ctx: commands.Context,
        member: discord.Member = commands.parameter(
            converter=commands.MemberConverter,
            description="Tag a member or include their ID",
        ),
        *,
        reason: str = commands.parameter(
            default="No reason provided",
            description="Optional reason (required for 7 days or more)",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return tick.end(warning="This command must be used in a server.")
        context = SnowflakeContext(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
        )
        msg = await server_mute_service.toggle_server_mute(
            context=context, member=member, reason=reason
        )
        return await tick.end(success=msg)

    @commands.command(name="smutes", help="List mutes.")
    async def list_server_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel, discord.Member, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):

        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.guild
        is_at_home = at_home(source=ctx)
        pages = await list_server_mutes.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @commands.command(name="stages", help="List stages.")
    @skip_text_command_help_discovery()
    async def list_stages_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel, discord.Member, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        pages = await automute_channel_service.build_pages(
            obj=obj, is_at_home=is_at_home
        )
        return await tick.end(success=pages)

    @commands.command(name="stream", help="Setup streaming.")
    @skip_text_command_help_discovery()
    async def modify_streaming_text_command(
        self,
        ctx: commands.Context,
        target_channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.TextChannelConverter,
            description="Tag a channel or include its ID.",
        ),
        source: Union[str, discord.abc.GuildChannel] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="`All` or tag a channel/include its ID.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        pages = await stream_service.toggle_stream(
            source=source,
            target_channel=target_channel,
        )
        return await tick.end(success=pages)

    @commands.command(name="streams", help="List streaming routes.")
    @skip_text_command_help_discovery()
    async def list_streaming_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[
            str, discord.Guild, discord.abc.GuildChannel
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: `all`, channel ID/mention, or server ID.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        pages = await list_streams.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @commands.command(name="vo", help="Start/stop video-only channel.")
    @skip_text_command_help_discovery()
    async def toggle_video_channel_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Tag a channel or include the ID",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        obj = channel or ctx.channel
        msg = await video_channel_service.toggle_video_channel(channel=obj)
        return await tick.end(success=msg)

    @commands.command(
        name="vos",
        help="List video channels.",
    )
    @skip_text_command_help_discovery()
    async def list_video_channels_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Include `all`, channel or server ID.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        pages = await list_video_channels.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @commands.command(name="xalias", help="Delete alias.")
    async def delete_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: str = commands.parameter(description="Include an alias name"),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return tick.end(warning="This command must be used in a guild.")
        context = SnowflakeContext(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
        )
        msg = await alias_service.delete_alias(alias_name=alias_name, context=context)
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
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        obj = channel or ctx.channel
        pages = await voice_mute_service.channel_unmute(channel=obj)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(AdminTextCommands(bot=bot))
