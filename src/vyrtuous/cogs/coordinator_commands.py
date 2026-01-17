"""coordinator_commands.py A discord.py cog containing coordinator commands for the Vyrtuous bot.

Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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

from datetime import datetime, timedelta, timezone
from typing import Optional

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.roles.moderator import Moderator
from vyrtuous.database.roles.vegan import Vegan
from vyrtuous.database.actions.voice_mute import VoiceMute
from vyrtuous.properties.snowflake import (
    AppChannelSnowflake,
    AppMemberSnowflake,
    ChannelSnowflake,
    MemberSnowflake,
)
from vyrtuous.service.logging_service import logger
from vyrtuous.service.check_service import (
    coordinator_predicator,
    has_equal_or_higher_role,
    not_bot,
)
from vyrtuous.service.messaging.message_service import MessageService
from vyrtuous.service.messaging.state_service import StateService
from vyrtuous.service.resolution.channel_service import resolve_channel
from vyrtuous.service.resolution.member_service import resolve_member
from vyrtuous.utils.emojis import get_random_emoji


class CoordinatorCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot, self.bot.db_pool)

    # DONE
    @app_commands.command(name="mod", description="Grant/revoke mods.")
    @app_commands.describe(
        member="Tag a member or include their ID",
        channel="Tag a channel or include its ID",
    )
    @coordinator_predicator()
    async def create_moderator_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake,
    ):
        state = StateService(interaction)
        try:
            channel_obj = await resolve_channel(
                ctx_interaction_or_message=interaction, channel_str=channel
            )
        except Exception as e:
            logger.warning(str(e).capitalize())
            channel_obj = interaction.channel
            await self.message_service.send_message(
                interaction,
                content=f"\U000026a0\U0000fe0f Defaulting to {channel_obj.mention}.",
            )
        try:
            member_obj = await resolve_member(
                ctx_interaction_or_message=interaction, member_str=member
            )
            not_bot(
                ctx_interaction_or_message=interaction, member_snowflake=member_obj.id
            )
            await has_equal_or_higher_role(
                ctx_interaction_or_message=interaction,
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id,
                sender_snowflake=interaction.user.id,
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
            
        action = None
        moderator = await Moderator.select(
            channel_snowflake=channel_obj.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=member_obj.id,
        )

        if moderator:
            await Moderator.delete(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id,
            )
            action = "revoked"
        else:
            moderator = Moderator(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id,
            )
            await moderator.create()
            action = "granted"

        try:
            return await state.end(
                success=f"{get_random_emoji()} "
                f"Moderator access for {member_obj.mention} has been "
                f"{action} in {channel_obj.mention}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @commands.command(name="mod", help="Grant/revoke mods.")
    @coordinator_predicator()
    async def create_moderator_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            default=None, description="Tag a member or include their ID"
        ),
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
    ):
        state = StateService(ctx)

        try:
            channel_obj = await resolve_channel(
                ctx_interaction_or_message=ctx, channel_str=channel
            )
        except Exception as e:
            logger.warning(str(e).capitalize())
            channel_obj = ctx.channel
            await self.message_service.send_message(
                ctx_interaction_or_message=ctx,
                content=f"\U000026a0\U0000fe0f Defaulting to {channel_obj.mention}.",
            )
        try:
            member_obj = await resolve_member(
                ctx_interaction_or_message=ctx, member_str=member
            )
            not_bot(ctx, member_snowflake=member_obj.id)
            await has_equal_or_higher_role(
                ctx_interaction_or_message=ctx,
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id,
                sender_snowflake=ctx.author.id,
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        action = None
        moderator = await Moderator.select(
            channel_snowflake=channel_obj.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=member_obj.id,
        )

        if moderator:
            await Moderator.delete(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id,
            )
            action = "revoked"
        else:
            moderator = Moderator(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id,
            )
            try:
                await moderator.create()
            except Exception as e:
                print(e)
            action = "granted"

        try:
            return await state.end(
                success=f"{get_random_emoji()} "
                f"Moderator access for {member_obj.mention} has been "
                f"{action} in {channel_obj.mention}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @app_commands.command(name="rmute", description="Room mute (except yourself).")
    @app_commands.describe(channel="Tag a channel or include its ID")
    @coordinator_predicator()
    async def room_mute_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        reason: Optional[str] = "No reason provided.",
    ):
        state = StateService(interaction)
        channel_obj = None
        muted_members, pages, skipped_members, failed_members = [], [], [], []
        try:
            channel_obj = await resolve_channel(
                ctx_interaction_or_message=interaction, channel_str=channel
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        for member in channel_obj.members:
            if member.id == interaction.user.id:
                continue
            voice_mute = await VoiceMute.select(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member.id,
                target="user",
            )
            if voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_obj.id:
                    try:
                        await member.edit(mute=True)
                    except Exception as e:
                        logger.warning(
                            f"Unable to voice-mute member "
                            f"{member.display_name} ({member.id}) in channel "
                            f"{channel_obj.name} ({channel_obj.id}) in guild "
                            f"{interaction.guild.name} ({interaction.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
            voice_mute = VoiceMute(
                channel_snowflake=channel_obj.id,
                expires_in=expires_in,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member.id,
                reason=reason,
                target="user",
            )
            await voice_mute.create()
            muted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_obj.mention}",
            f"**Muted:** {len(muted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f'**Skipped:** {len(channel_obj.members) \
                - len(muted_members) \
                - len(failed_members)
            }',
        ]
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{get_random_emoji()} Room Mute Summary",
            color=discord.Color.blurple(),
        )
        pages.append(embed)

        await StateService.send_pages(obj=Vegan, pages=pages, state=state)

    # DONE
    @commands.command(name="rmute", help="Room mute (except yourself).")
    @coordinator_predicator()
    async def room_mute_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
        reason: str = commands.parameter(
            default="No reason provided.", description="Specify a reason."
        ),
    ):
        state = StateService(ctx)
        channel_obj = None
        muted_members, pages, skipped_members, failed_members = [], [], [], []
        try:
            channel_obj = await resolve_channel(
                ctx_interaction_or_message=ctx, channel_str=channel
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        for member in channel_obj.members:
            if member.id == ctx.author.id:
                continue
            voice_mute = await VoiceMute.select(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member.id,
                target="user",
            )
            if voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_obj.id:
                    try:
                        await member.edit(mute=True)
                    except Exception as e:
                        logger.warning(
                            f"Unable to voice-mute member "
                            f"{member.display_name} ({member.id}) in channel "
                            f"{channel_obj.name} ({channel_obj.id}) in guild "
                            f"{ctx.guild.name} ({ctx.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
            voice_mute = VoiceMute(
                channel_snowflake=channel_obj.id,
                expires_in=expires_in,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member.id,
                reason=reason,
                target="user",
            )
            await voice_mute.create()
            muted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_obj.mention}",
            f"**Muted:** {len(muted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f'**Skipped:** {len(channel_obj.members) \
                - len(muted_members) \
                - len(failed_members)
            }',
        ]
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{get_random_emoji()} Room Mute Summary",
            color=discord.Color.blurple(),
        )
        pages.append(embed)

        await StateService.send_pages(obj=Vegan, pages=pages, state=state)

    # DONE
    @app_commands.command(name="xrmute", description="Unmute all.")
    @app_commands.describe(channel="Tag a channel or include its ID")
    @coordinator_predicator()
    async def room_unmute_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        state = StateService(interaction)
        channel_obj = None
        failed_members, pages, skipped_members, unmuted_members = [], [], [], []
        try:
            channel_obj = await resolve_channel(
                ctx_interaction_or_message=interaction, channel_str=channel
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        for member in channel_obj.members:
            voice_mute = await VoiceMute.select(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member.id,
                target="user",
            )
            if not voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_obj.id:
                    try:
                        await member.edit(mute=False)
                    except Exception as e:
                        logger.warning(
                            f"Unable to undo voice-mute "
                            f"for member {member.display_name} ({member.id}) "
                            f"in channel {channel_obj.name} ({channel_obj.id}) "
                            f"in guild {interaction.guild.name} "
                            f"({interaction.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            await VoiceMute.delete(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                target="user",
            )
            unmuted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_obj.mention}",
            f"**Unmuted:** {len(unmuted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f'**Skipped:** {len(channel_obj.members) \
                - len(unmuted_members) \
                - len(failed_members)
            }',
        ]
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{get_random_emoji()} Room Unmute Summary",
            color=discord.Color.blurple(),
        )
        pages.append(embed)

        await StateService.send_pages(obj=VoiceMute, pages=pages, state=state)

    # DONE
    @commands.command(name="xrmute", help="Unmute all.")
    @coordinator_predicator()
    async def room_unmute_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
    ):
        state = StateService(ctx)
        channel_obj = None
        failed_members, pages, skipped_members, unmuted_members = [], [], [], []
        try:
            channel_obj = await resolve_channel(
                ctx_interaction_or_message=ctx, channel_str=channel
            )
        except Exception as e:
            return await state.end(
                warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
            )
        for member in channel_obj.members:
            voice_mute = await VoiceMute.select(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member.id,
                target="user",
            )
            if not voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_obj.id:
                    try:
                        await member.edit(mute=False)
                    except Exception as e:
                        logger.warning(
                            "Unable to voice-mute member "
                            f"{member.display_name} ({member.id}) in channel "
                            f"{channel_obj.name} ({channel_obj.id}) in guild "
                            f"{ctx.guild.name} ({ctx.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            await VoiceMute.delete(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                target="user",
            )
            unmuted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_obj.mention}",
            f"**Unmuted:** {len(unmuted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f'**Skipped:** {len(channel_obj.members) \
                - len(unmuted_members) \
                - len(failed_members)
            }',
        ]
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{get_random_emoji()} Room Unmute Summary",
            color=discord.Color.blurple(),
        )
        pages.append(embed)

        await StateService.send_pages(obj=VoiceMute, pages=pages, state=state)


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorCommands(bot))
