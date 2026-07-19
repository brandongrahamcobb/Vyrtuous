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
from vyrtuous.cache.registry import MemberState
from vyrtuous.inc.helpers import at_home
from vyrtuous.listing import (
    list_aliases,
    list_bans,
    list_coordinators,
    list_flags,
    list_moderators,
    list_text_mutes,
    list_voice_mutes,
)
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import moderator_service


class ModeratorTextCommands(commands.Cog):

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

    @commands.command(name="bans", help="List bans.")
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, discord.Member, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: channel ID/mention, member ID/mention or server ID.",
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
        obj = target or ctx.channel
        pages = await list_bans.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="cmds", help="List aliases.")
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must target a valid server.")
        if guild is None:
            guild_snowflake = ctx.guild.id
        else:
            guild_snowflake = guild.id
        obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        pages = await list_aliases.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @commands.command(name="coords", help="Lists coords.")
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: `all`, channel ID/mention, or server ID.",
        ),
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
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
        pages = await list_coordinators.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @commands.command(name="del", help="Delete message.")
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        msg: discord.Message = commands.parameter(
            converter=commands.MessageConverter, description="Message snowflake"
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if msg.channel.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        await moderator_service.has_equal_or_lower_role(
            target_member_snowflake=int(msg.author.id),
            member_snowflake=ctx.author.id,
            channel_snowflake=msg.channel.id,
            guild_snowflake=msg.channel.guild.id,
        )
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await tick.end(error=str(e).capitalize())
        return await tick.end(success=f"Message `{msg.id}` deleted successfully.")

    @commands.command(name="flags", help="List flags.")
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            int, discord.abc.GuildChannel, discord.Guild, discord.Member, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: channel ID/mention, member ID/mention, or server ID.",
        ),
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
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
        pages = await list_flags.build_pages(guild_snowflake=guild_snowflake, obj=obj)
        return await tick.end(success=pages)

    @commands.command(name="mods", help="Lists mods.")
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
        guild: Union[discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify a server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if guild is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
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
        pages = await list_moderators.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @commands.command(name="mutes", help="List mutes.")
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
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
        if guild is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
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
        pages = await list_voice_mutes.build_pages(
            guild_snowflake=guild_snowflake, obj=obj
        )
        return await tick.end(success=pages)

    @commands.command(name="purge", help="Delete messages.")
    async def purge_text_command(
        self,
        ctx: commands.Context,
        member: Union[int, discord.Member] = commands.parameter(
            converter=MultiConverter,
            description="Tag a member or include their ID",
        ),
        amount: int = commands.parameter(
            default=25, description="Number of messages to delete"
        ),
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Specify the channel or ID",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role="Coordinator",
        )
        if isinstance(member, discord.Member):
            member_snowflake = int(member.id)
            display_name = str(member.mention)
        else:
            member_snowflake = int(member)
            simplified_member = self.__bot.registry.get(MemberState).active.get(
                member_snowflake, None
            )
            if simplified_member:
                display_name = simplified_member[0]
            else:
                raise commands.MemberNotFound(str(member))
        await moderator_service.has_equal_or_lower_role(
            target_member_snowflake=int(member_snowflake),
            member_snowflake=ctx.author.id,
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
        )
        target = channel or ctx.channel
        if not isinstance(
            target, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
        ):
            return await tick.end(warning="This command must target a valid channel.")
        count = int(0)
        async for msg in target.history():
            if amount == count:
                break
            if msg.author.id == member_snowflake:
                await msg.delete()
                count += 1
        return await tick.end(
            success=f"Successfully deleted {count} messages from {display_name} in {target.mention}."
        )

    @commands.command(name="summary", help="List user moderation.")
    async def list_moderation_summary_text_command(
        self,
        ctx: commands.Context,
        member: discord.Member = commands.parameter(
            converter=commands.MemberConverter,
            description="Specify a member ID/mention.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        pages: list[discord.Embed] = []
        obj = member
        is_at_home = at_home(source=ctx)
        services = []
        services.append(list_bans)
        services.append(list_flags)
        services.append(list_text_mutes)
        services.append(list_voice_mutes)
        for service in services:
            summary_pages = await service.build_pages(obj=obj, is_at_home=is_at_home)
            if isinstance(summary_pages, list):
                for page in summary_pages:
                    if isinstance(page, discord.Embed):
                        pages.append(page)
        if not pages:
            return await tick.end(success="No infractions found")
        return await tick.end(success=pages)

    @commands.command(name="tmutes", help="List text-mutes.")
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Union[
            str, discord.abc.GuildChannel, discord.Guild, None
        ] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
        else:
            obj = target or ctx.channel
        is_at_home = at_home(source=ctx)
        pages = await list_text_mutes.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(ModeratorTextCommands(bot=bot))
