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

from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.roles.coordinator import coordinator_predicator
from vyrtuous.db.roles.moderator import Moderator
from vyrtuous.db.actions.voice_mute import VoiceMute
from vyrtuous.fields.snowflake import (
    ChannelSnowflake,
    MemberSnowflake,
)
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import DiscordObject


class CoordinatorTextCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    @commands.command(name="mod", help="Grant/revoke mods.")
    @coordinator_predicator()
    async def toggle_moderator_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            description="Tag a member or include their ID"
        ),
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
    ):
        state = StateService(ctx=ctx)
        snowflake_kwargs = {
            "channel_snowflake": ctx.channel.id,
            "guild_snowflake": ctx.guild.id,
            "member_snowflake": ctx.author.id,
        }
        do = DiscordObject(ctx=ctx)
        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        msg = await Moderator.toggle_moderator(
            channel_dict=channel_dict,
            member_dict=member_dict,
            snowflake_kwargs=snowflake_kwargs,
        )
        return await state.end(success=msg)

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
        state = StateService(ctx=ctx)
        snowflake_kwargs = {
            "channel_snowflake": ctx.channel.id,
            "guild_snowflake": ctx.guild.id,
            "member_snowflake": ctx.author.id,
        }
        do = DiscordObject(ctx=ctx)
        channel_dict = await do.determine_from_target(target=channel)
        pages = VoiceMute.room_mute(
            channel_dict=channel_dict,
            guild_snowflake=ctx.guild.id,
            reason=reason,
            snowflake_kwargs=snowflake_kwargs,
        )
        await StateService.send_pages(plural=VoiceMute.PLURAL, pages=pages, state=state)

    @commands.command(name="xrmute", help="Unmute all.")
    @coordinator_predicator()
    async def room_unmute_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        channel_dict = await do.determine_from_target(target=channel)
        pages = VoiceMute.room_unmute(
            channel_dict=channel_dict, guild_snowflake=ctx.guild.id
        )
        await StateService.send_pages(plural=VoiceMute.PLURAL, pages=pages, state=state)


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorTextCommands(bot))
