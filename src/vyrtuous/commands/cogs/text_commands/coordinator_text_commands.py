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
from vyrtuous.commands.cogs.help_command import skip_help_discovery
from vyrtuous.commands.discord_object_service import DiscordObject
from vyrtuous.commands.fields.duration import Duration, DurationObject
from vyrtuous.commands.fields.snowflake import ChannelSnowflake, MemberSnowflake
from vyrtuous.commands.messaging.message_service import MessageService
from vyrtuous.commands.messaging.state_service import StateService
from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.roles.coord.coordinator import Coordinator
from vyrtuous.db.roles.coord.coordinator_service import coordinator_predicator
from vyrtuous.db.roles.mod.moderator_service import ModeratorService
from vyrtuous.db.rooms.stage.stage_service import StageService


class CoordinatorTextCommands(commands.Cog):

    ROLE = Coordinator

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
        default_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        do = DiscordObject(ctx=ctx)
        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(channel_dict.get("columns", None))
        await PermissionService.has_equal_or_lower_role(
            member_snowflake=int(member_dict.get("id", None)),
            updated_kwargs=updated_kwargs,
        )
        msg = await ModeratorService.toggle_moderator(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            member_dict=member_dict,
        )
        return await state.end(success=msg)

    @commands.command(name="stage", help="Start/stop stage")
    @coordinator_predicator()
    @skip_help_discovery()
    async def toggle_stage_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
        *,
        duration: Duration = commands.parameter(
            default=DurationObject("1h"),
            description="Options: (+|-)duration(m|h|d) "
            "0 - permanent / 24h - default",
        ),
    ):
        state = StateService(ctx=ctx)
        default_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        do = DiscordObject(ctx=ctx)
        channel = channel or int(ctx.channel.id)
        channel_dict = await do.determine_from_target(target=channel)
        pages = await StageService.toggle_stage(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            duration=duration,
        )
        await StateService.send_pages(title="Stage", pages=pages, state=state)


async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorTextCommands(bot))
