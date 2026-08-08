"""!/bin/python3

hidden_admin_app_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.models.message import AppMessage, MessageObject
from vyrtuous.models.metadata import metadata
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import MemberNotFound
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.tick import Tick


class UtilityAppCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def cog_app_command_error(self, interaction, error):
        self.__bot.logger.info(str(error))
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error), ephemeral=True)

    @metadata(permission="command.utility.delete")
    @app_commands.command(name="del", description="Delete a message.")
    @app_commands.describe(msg="Specify a message ID.")
    async def delete_message_app_command(
        self,
        interaction: discord.Interaction,
        msg: app_commands.Transform[MessageObject, AppMessage],
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.defer()
        if msg.message.channel.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        else:
            channel_snowflake = msg.message.channel.id
            guild_snowflake = msg.message.channel.guild.id
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            requested=["command.utility.delete"],
        )
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            member_snowflake=int(msg.message.author.id),
            author_snowflake=interaction.user.id,
            channel_snowflake=msg.message.channel.id,
            guild_snowflake=msg.message.channel.guild.id,
        )
        try:
            await msg.message.delete()
        except discord.Forbidden as e:
            return await tick.end(error=str(e).capitalize())
        return await tick.end(
            success=f"Message `{msg.message.id}` deleted successfully."
        )

    @metadata(permission="command.utility.ping")
    @app_commands.command(name="ping", description="Ping me!")
    async def ping_app_command(
        self, interaction: discord.Interaction
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.defer()
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            requested=["command.utility.ping"],
        )
        return await tick.end(success="Pong!")

    @metadata(permission="command.utility.purge")
    @app_commands.command(name="purge", description="Delete messages.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
        amount="Specify the number of messages.",
        channel="Specify a channel ID/mention.",
    )
    async def purge_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject | None, AppTarget] = None,
        amount: int = 25,
        channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.defer()
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if channel is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            if interaction.channel is None or not isinstance(
                interaction.channel,
                (discord.TextChannel, discord.VoiceChannel, discord.StageChannel),
            ):
                return await tick.end(
                    warning="This command must target a valid server channel."
                )
            channel_obj = interaction.channel
            guild_snowflake = interaction.guild.id
        else:
            if not isinstance(
                channel.target,
                (discord.TextChannel, discord.VoiceChannel, discord.StageChannel),
            ):
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_obj = channel.target
            guild_snowflake = channel.target.guild.id
        if member is None:
            member_snowflake = None
            display_name = None
        else:
            if isinstance(member.target, int):
                member_snowflake = member.target
                simplified_member = self.__bot.registry.get(MemberState).active.get(
                    member_snowflake, None
                )
                if simplified_member:
                    display_name = simplified_member[0]
                else:
                    raise MemberNotFound(str(member))
            elif isinstance(member.target, discord.Member):
                member_snowflake = member.target.id
                display_name = str(member.target.mention)
            else:
                return await tick.end(
                    warning=f"This command must target a valid member."
                )
        count = 0
        skipped = 0
        async for msg in channel_obj.history():
            if count == 0:
                count += 1
                continue
            if amount == count:
                break
            try:
                await permission_service.has_equal_or_lower_role(
                    permission_state=permission_state,
                    member_snowflake=int(msg.author.id),
                    author_snowflake=interaction.user.id,
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=guild_snowflake,
                )
                await permission_service.has_permissions(
                    permission_state=permission_state,
                    member_snowflake=interaction.user.id,
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=guild_snowflake,
                    requested=["command.utility.purge"],
                )
                if member_snowflake is not None:
                    if msg.author.id == member_snowflake:
                        if interaction.message:
                            if msg.id == interaction.message.id:
                                continue
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

    @metadata(permission="command.utility.move")
    @app_commands.command(name="rmove", description="VC move.")
    @app_commands.describe(
        target_channel="Specify a `to` channel ID/mention.",
        source_channel="Specify a `from` channel ID/mention.",
    )
    async def channel_move_all_app_command(
        self,
        interaction: discord.Interaction,
        target_channel: app_commands.Transform[TargetObject, AppTarget],
        source_channel: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.defer()
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        if isinstance(
            target_channel.target,
            (discord.VoiceChannel, discord.StageChannel),
        ):
            target_guild_snowflake = target_channel.target.guild.id
            target_guild_name = target_channel.target.guild.name
            target_channel_snowflake = target_channel.target.id
            target_channel_name = target_channel.target.name
            target_channel_obj = target_channel.target
            target_channel_mention = target_channel.target.mention
        else:
            return await tick.end(
                warning="This command must target a valid target channel."
            )
        if source_channel is None:
            if interaction.channel is None or not isinstance(
                interaction.channel, (discord.VoiceChannel, discord.StageChannel)
            ):
                return await tick.end(
                    warning="This command must target a valid source channel."
                )
            source_guild_snowflake = interaction.channel.guild.id
            source_guild_name = interaction.channel.guild.name
            source_channel_snowflake = interaction.channel.id
            source_channel_name = interaction.channel.name
            source_channel_members = interaction.channel.members
            source_channel_mention = interaction.channel.mention
        elif isinstance(
            source_channel.target,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            source_guild_snowflake = source_channel.target.guild.id
            source_guild_name = source_channel.target.guild.name
            source_channel_snowflake = source_channel.target.id
            source_channel_name = source_channel.target.name
            source_channel_members = source_channel.target.members
            source_channel_mention = source_channel.target.mention
        else:
            return await tick.end(
                warning="This command must target a valid source channel."
            )
        await permission_service.has_permissions(
            permission_state=permission_state,
            channel_snowflake=target_channel_snowflake,
            member_snowflake=interaction.user.id,
            guild_snowflake=target_guild_snowflake,
            requested=["command.utility.move"],
        )
        await permission_service.has_permissions(
            permission_state=permission_state,
            channel_snowflake=source_channel_snowflake,
            member_snowflake=interaction.user.id,
            guild_snowflake=source_guild_snowflake,
            requested=["command.utility.move"],
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
    await bot.add_cog(UtilityAppCommands(bot=bot))
