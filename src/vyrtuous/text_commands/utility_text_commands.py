"""!/bin/python3
sysadmin_text_commands.py A discord.py cog containing sysadmin commands for the Vyrtuous bot.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.utils.errors.error import MemberNotFound
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.permissions import permission_service


class UtilityTextCommands(commands.Cog):

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

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
        else:
            channel_snowflake = msg.channel.id
            guild_snowflake = msg.channel.guild.id
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.utility.delete"],
        )
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            member_snowflake=int(msg.author.id),
            author_snowflake=ctx.author.id,
            channel_snowflake=msg.channel.id,
            guild_snowflake=msg.channel.guild.id,
        )
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await tick.end(error=str(e).capitalize())
        return await tick.end(success=f"Message `{msg.id}` deleted successfully.")

    @commands.command(name="ping", help="Ping me!")
    async def ping_text_command(self, ctx: commands.Context) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be executed in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            requested=["command.utility.ping"],
        )
        return await tick.end(success="Pong!")

    @commands.command(name="purge", help="Delete messages.")
    async def purge_text_command(
        self,
        ctx: commands.Context,
        member: Union[int, discord.Member, None] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify the member ID/mention.",
        ),
        amount: int = commands.parameter(
            default=25, description="Number of messages to delete."
        ),
        channel: discord.abc.GuildChannel | None = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify the channel ID/mention.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if channel is None:
            if ctx.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            if ctx.channel is None or not isinstance(
                ctx.channel,
                (discord.TextChannel, discord.VoiceChannel, discord.StageChannel),
            ):
                return await tick.end(
                    warning="This command must target a valid server channel."
                )
            channel_obj = ctx.channel
            guild_snowflake = ctx.guild.id
        else:
            if not isinstance(
                channel,
                (discord.TextChannel, discord.VoiceChannel, discord.StageChannel),
            ):
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_obj = channel
            guild_snowflake = channel.guild.id
        if member is None:
            member_snowflake = None
            display_name = None
        else:
            if isinstance(member, int):
                member_snowflake = member
                simplified_member = self.__bot.registry.get(MemberState).active.get(
                    member_snowflake, None
                )
                if simplified_member:
                    display_name = simplified_member[0]
                else:
                    raise MemberNotFound(str(member))
            elif isinstance(member, discord.Member):
                member_snowflake = member.id
                display_name = str(member.mention)
        count = 0
        skipped = 0
        async for msg in channel_obj.history():
            if amount == count:
                break
            try:
                await permission_service.has_equal_or_lower_role(
                    permission_state=permission_state,
                    member_snowflake=int(msg.author.id),
                    author_snowflake=ctx.author.id,
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=guild_snowflake,
                )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=ctx.author.id,
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=guild_snowflake,
                    requested=["command.utility.purge"],
                )
                if member_snowflake is not None:
                    if msg.author.id == member_snowflake:
                        await msg.delete()
                        count += 1
                else:
                    await msg.delete()
                    count += 1
            except Exception:
                skipped += 1
                continue
        message = f"Deleted {count} messages. Skipped {skipped} messages"
        if member_snowflake is not None and channel_obj is not None:
            message = f"Deleted {count} and skipped {skipped} messages from {display_name} in {channel_obj.mention}."
        elif member_snowflake is not None:
            message = (
                f"Deleted {count} and skipped {skipped} messages from {display_name}."
            )
        elif channel_obj is not None:
            message = f"Deleted {count} and skipped {skipped} messages in {channel_obj.mention}."
        return await tick.end(success=message)

    @commands.command(name="rmv", help="VC move.")
    async def channel_move_all_text_command(
        self,
        ctx: commands.Context,
        target_channel: (
            discord.VoiceChannel | discord.StageChannel | None
        ) = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a `to` channel ID/mention.",
        ),
        source_channel: (
            discord.VoiceChannel | discord.StageChannel | None
        ) = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify a `from` channel ID/mention.",
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
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
        else:
            return await tick.end(
                warning="This command must target a valid target channel."
            )
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
        else:
            return await tick.end(
                warning="This command must target a valid source channel."
            )
        await permission_service.has_permissions(
            permission_state=permission_state,
            channel_snowflake=target_channel_snowflake,
            member_snowflake=ctx.author.id,
            guild_snowflake=target_guild_snowflake,
            requested=["command.movement.channel_move"],
        )
        await permission_service.has_permissions(
            permission_state=permission_state,
            channel_snowflake=source_channel_snowflake,
            member_snowflake=ctx.author.id,
            guild_snowflake=source_guild_snowflake,
            requested=["command.movement.channel_move"],
        )
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


async def setup(bot: DiscordBot):
    await bot.add_cog(UtilityTextCommands(bot=bot))
