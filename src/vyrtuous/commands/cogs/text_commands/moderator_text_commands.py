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

from pathlib import Path

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.cogs.help_command import skip_help_discovery
from vyrtuous.commands.discord_object_service import DiscordObject
from vyrtuous.commands.fields.snowflake import (
    ChannelSnowflake,
    MemberSnowflake,
    MessageSnowflake,
)
from vyrtuous.commands.home import at_home
from vyrtuous.commands.messaging.message_service import MessageService
from vyrtuous.commands.messaging.state_service import StateService
from vyrtuous.db.alias.alias_service import AliasService
from vyrtuous.db.infractions.ban.ban_service import BanService
from vyrtuous.db.infractions.flag.flag_service import FlagService
from vyrtuous.db.infractions.tmute.text_mute_service import TextMuteService
from vyrtuous.db.infractions.vmute.voice_mute_service import VoiceMuteService
from vyrtuous.db.roles.admin.administrator_service import AdministratorService
from vyrtuous.db.roles.coord.coordinator_service import CoordinatorService
from vyrtuous.db.roles.mod.moderator import Moderator
from vyrtuous.db.roles.mod.moderator_service import (
    ModeratorService,
    moderator_predicator,
)
from vyrtuous.db.roles.vegan.vegan_service import VeganService
from vyrtuous.db.rooms.stage.stage_service import StageService
from vyrtuous.db.rooms.temp.temporary_room_service import TemporaryRoomService
from vyrtuous.utils.dir_to_classes import dir_to_classes


class ModeratorTextCommands(commands.Cog):

    ROLE = Moderator

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.message_service = MessageService(self.bot)

    @commands.command(name="admins", help="Lists admins.")
    @moderator_predicator()
    @skip_help_discovery()
    async def list_administrators_text_command(
        self,
        ctx: commands.Context,
        *,
        target: str | None = commands.parameter(
            default=None,
            description="Specify one of: `all`, " "channel ID/mention or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        target = target or int(ctx.guild.id)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=str(target))
        pages = await AdministratorService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Administrators", pages=pages, state=state)

    @commands.command(name="bans", help="List bans.")
    @moderator_predicator()
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        target: str | None = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        target = target or int(ctx.channel.id)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await BanService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Bans", pages=pages, state=state)

    @commands.command(name="cmds", help="List aliases.")
    @moderator_predicator()
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        target: str | None = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        target = target or int(ctx.channel.id)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await AliasService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Aliases", pages=pages, state=state)

    @commands.command(name="coords", help="Lists coords.")
    @moderator_predicator()
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        target: str | None = commands.parameter(
            default=None,
            description="Specify one of: `all`, channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        target = target or int(ctx.channel.id)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await CoordinatorService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Coordinators", pages=pages, state=state)

    @commands.command(name="del", help="Delete message.")
    @moderator_predicator()
    @skip_help_discovery()
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        message: MessageSnowflake = commands.parameter(description="Message snowflake"),
    ):
        state = StateService(ctx=ctx)
        try:
            msg = await ctx.channel.fetch_message(message)
        except discord.NotFound:
            return await state.end(warning=f"Message `{message}` not found.")
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())
        return await state.end(success=f"Message `{message}` deleted successfully.")

    @commands.command(name="flags", help="List flags.")
    @moderator_predicator()
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        target: str | None = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        target = target or int(ctx.channel.id)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await FlagService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Flags", pages=pages, state=state)

    @commands.command(name="ls", help="List new vegans.")
    @moderator_predicator()
    @skip_help_discovery()
    async def list_new_vegans_text_command(
        self,
        ctx: commands.Context,
        *,
        target: str | None = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        target = target or int(ctx.guild.id)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await VeganService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Vegans", pages=pages, state=state)

    @commands.command(
        name="migrate",
        help="Migrate a temporary room to a new channel by snowflake.",
    )
    @moderator_predicator()
    @skip_help_discovery()
    async def migrate_temp_room_text_command(
        self,
        ctx: commands.Context,
        old_name: str = commands.parameter(description="Provide a channel name"),
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
        msg = await TemporaryRoomService.migrate_temporary_room(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            old_name=old_name,
        )
        return await state.end(success=msg)

    @commands.command(name="mods", help="Lists mods.")
    @moderator_predicator()
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        target: str | None = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        target = target or int(ctx.channel.id)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await ModeratorService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Moderators", pages=pages, state=state)

    @commands.command(name="mutes", help="List mutes.")
    @moderator_predicator()
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        target: str | None = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        target = target or int(ctx.channel.id)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await VoiceMuteService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Voice Mutes", pages=pages, state=state)

    @commands.command(name="mstage", help="Toggle stage mute/unmute.")
    @moderator_predicator()
    @skip_help_discovery()
    async def stage_mute_text_command(
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
        channel = channel or int(ctx.channel.id)
        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(channel_dict.get("columns", None))
        updated_kwargs.update(member_dict.get("columns", None))
        msg = await StageService.toggle_stage_mute(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            member_dict=member_dict,
        )
        await state.end(success=msg)

    @commands.command(name="summary", help="List user moderation.")
    @moderator_predicator()
    async def list_moderation_summary_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            description="Specify a member ID/mention."
        ),
    ):
        pages = []
        dir_paths = []
        dir_paths.append(Path("/app/vyrtuous/db/infractions"))

        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        member_dict = await do.determine_from_target(target=member)
        for obj in dir_to_classes(dir_paths=dir_paths, parent=AliasService):
            object_pages = await obj.build_pages(
                object_dict=member_dict, is_at_home=is_at_home
            )
            if object_pages:
                pages.extend(object_pages)
        await StateService.send_pages(title="Infractions", pages=pages, state=state)

    @commands.command(name="survey", help="Survey stage members.")
    @moderator_predicator()
    @skip_help_discovery()
    async def stage_survey_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        channel = channel or int(ctx.channel.id)
        channel_dict = await do.determine_from_target(target=channel)
        pages = await StageService.survey(
            channel_dict=channel_dict, guild_snowflake=ctx.guild.id
        )
        await StateService.send_pages(title="Stage Roles", pages=pages, state=state)

    @commands.command(name="tmutes", help="List text-mutes.")
    @moderator_predicator()
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        target: str | None = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        target = target or int(ctx.channel.id)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await TextMuteService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Text Mutes", pages=pages, state=state)


async def setup(bot: DiscordBot):
    cog = ModeratorTextCommands(bot)
    await bot.add_cog(cog)
