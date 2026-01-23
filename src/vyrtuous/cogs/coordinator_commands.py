"""coordinator_commands.py A discord.py cog containing coordinator commands for the Vyrtuous bot.

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

from datetime import datetime, timedelta, timezone

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.roles.coordinator import coordinator_predicator
from vyrtuous.db.roles.moderator import Moderator
from vyrtuous.db.roles.vegan import Vegan
from vyrtuous.db.actions.voice_mute import VoiceMute
from vyrtuous.fields.snowflake import (
    AppChannelSnowflake,
    AppMemberSnowflake,
    ChannelSnowflake,
    MemberSnowflake,
)
from vyrtuous.utils.logger import logger
from vyrtuous.utils.check import (
    has_equal_or_lower_role_wrapper,
)
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import DiscordObject
from vyrtuous.utils.emojis import get_random_emoji


class CoordinatorCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    # DONE
    @app_commands.command(name="mod", description="Grant/revoke mods.")
    @app_commands.describe(
        member="Tag a member or include their ID",
        channel="Tag a channel or include its ID.",
    )
    @coordinator_predicator()
    async def create_moderator_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake,
    ):
        action = None

        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)
        if not isinstance(channel_dict.get("object", None), discord.abc.GuildChannel):
            return await state.end(
                warning=f"Invalid channel ID ({channel})."
            )
        member_dict = await do.determine_from_target(target=member)
        if not isinstance(member_dict.get("object", None), discord.Member):
            return await state.end(
                warning=f"Invalid member ID ({channel})."
            )
        await has_equal_or_lower_role_wrapper(
            source=interaction,
            member_snowflake=member_dict.get("id", None),
            sender_snowflake=interaction.user.id,
        )
        kwargs = {}
        kwargs.update(channel_dict.get("columns", None))
        kwargs.update(member_dict.get("columns", None))

        moderator = await Moderator.select(**kwargs)
        if moderator:
            await Moderator.delete(**kwargs)
            action = "revoked"
        else:
            moderator = Moderator(**kwargs)
            await moderator.create()
            action = "granted"

        return await state.end(
            success=f"Moderator access for {member_dict.get("mention", None)} has been "
            f"{action} in {channel_dict.get("mention", None)}."
        )

    # DONE
    @commands.command(name="mod", help="Grant/revoke mods.")
    @coordinator_predicator()
    async def create_moderator_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            description="Tag a member or include their ID"
        ),
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
    ):
        action = None

        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)
        if not isinstance(channel_dict.get("object", None), discord.abc.GuildChannel):
            return await state.end(
                warning=f"Invalid channel ID ({channel})."
            )
        member_dict = await do.determine_from_target(target=member)
        if not isinstance(member_dict.get("object", None), discord.Member):
            return await state.end(
                warning=f"Invalid member ID ({channel})."
            )
        await has_equal_or_lower_role_wrapper(
            source=ctx,
            member_snowflake=member_dict.get("id", None),
            sender_snowflake=ctx.author.id,
        )
        kwargs = {}
        kwargs.update(channel_dict.get("columns", None))
        kwargs.update(member_dict.get("columns", None))

        moderator = await Moderator.select(**kwargs)
        if moderator:
            await Moderator.delete(**kwargs)
            action = "revoked"
        else:
            moderator = Moderator(**kwargs)
            await moderator.create()
            action = "granted"

        return await state.end(
            success=f"Moderator access for {member_dict.get("mention", None)} has been "
            f"{action} in {channel_dict.get("mention", None)}."
        )

    # DONE
    @app_commands.command(name="rmute", description="Room mute (except yourself).")
    @app_commands.describe(channel="Tag a channel or include its ID.")
    @coordinator_predicator()
    async def room_mute_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        reason: str = "No reason provided.",
    ):
        muted_members, pages, skipped_members, failed_members = [], [], [], []

        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)
        kwargs = channel_dict.get("columns", None)

        for member in channel_dict.get("object", None).members:
            if member.id == interaction.user.id:
                continue
            voice_mute = await VoiceMute.select(
                **kwargs,
                target="user",
            )
            if voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_dict.get("id", None):
                    try:
                        await member.edit(mute=True)
                    except Exception as e:
                        logger.warning(
                            f"Unable to voice-mute member "
                            f"{member.display_name} ({member.id}) in channel "
                            f"{channel_dict.get('name', None)} ({channel_dict.get('id', None)}) in guild "
                            f"{interaction.guild.name} ({interaction.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
            voice_mute = VoiceMute(
                expires_in=expires_in,
                member_snowflake=member.id,
                reason=reason,
                target="user",
                **kwargs,
            )
            await voice_mute.create()
            muted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_dict.get("mention", None)}",
            f"**Muted:** {len(muted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f'**Skipped:** {len(channel_dict.get('object', None).members) \
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
            description="Tag a channel or include its ID."
        ),
        *,
        reason: str = commands.parameter(
            default="No reason provided.", description="Specify a reason."
        ),
    ):
        muted_members, pages, skipped_members, failed_members = [], [], [], []
        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)
        kwargs = channel_dict.get("columns", None)

        for member in channel_dict.get("object", None).members:
            if member.id == ctx.author.id:
                continue
            voice_mute = await VoiceMute.select(
                **kwargs,
                target="user",
            )
            if voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_dict.get("id", None):
                    try:
                        await member.edit(mute=True)
                    except Exception as e:
                        logger.warning(
                            f"Unable to voice-mute member "
                            f"{member.display_name} ({member.id}) in channel "
                            f"{channel_dict.get('name', None)} ({channel_dict.get('id', None)}) in guild "
                            f"{ctx.guild.name} ({ctx.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
            voice_mute = VoiceMute(
                expires_in=expires_in,
                member_snowflake=member.id,
                reason=reason,
                target="user",
                **kwargs,
            )
            await voice_mute.create()
            muted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_dict.get("mention", None)}",
            f"**Muted:** {len(muted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f'**Skipped:** {len(channel_dict.get('object', None).members) \
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
    @app_commands.describe(channel="Tag a channel or include its ID.")
    @coordinator_predicator()
    async def room_unmute_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        unmuted_members, pages, skipped_members, failed_members = [], [], [], []

        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)
        kwargs = channel_dict.get("columns", None)

        for member in channel_dict.get("object", None).members:
            voice_mute = await VoiceMute.select(target="user", **kwargs)
            if not voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_dict.get("id", None):
                    try:
                        await member.edit(mute=False)
                    except Exception as e:
                        logger.warning(
                            f"Unable to undo voice-mute "
                            f"for member {member.display_name} ({member.id}) "
                            f"in channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) "
                            f"in guild {interaction.guild.name} "
                            f"({interaction.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            await VoiceMute.delete(target="user", **kwargs)
            unmuted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_dict.get("mention", None)}",
            f"**Unmuted:** {len(unmuted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f'**Skipped:** {len(channel_dict.get('object', None).members) \
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
            description="Tag a channel or include its ID."
        ),
    ):
        unmuted_members, pages, skipped_members, failed_members = [], [], [], []

        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)
        kwargs = channel_dict.get("columns", None)

        for member in channel_dict.get("object", None).members:
            voice_mute = await VoiceMute.select(target="user", **kwargs)
            if not voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_dict.get("id", None):
                    try:
                        await member.edit(mute=False)
                    except Exception as e:
                        logger.warning(
                            f"Unable to undo voice-mute "
                            f"for member {member.display_name} ({member.id}) "
                            f"in channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) "
                            f"in guild {ctx.guild.name} "
                            f"({ctx.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            await VoiceMute.delete(target="user", **kwargs)
            unmuted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_dict.get("mention", None)}",
            f"**Unmuted:** {len(unmuted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f'**Skipped:** {len(channel_dict.get('object', None).members) \
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
