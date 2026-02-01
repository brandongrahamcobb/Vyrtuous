"""moderator_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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
from vyrtuous.cogs.help_command import skip_help_discovery
from vyrtuous.fields.snowflake import (
    ChannelSnowflake,
    MemberSnowflake,
    MessageSnowflake,
)
from vyrtuous.service.discord_object_service import DiscordObject
from vyrtuous.service.infractions.ban_service import BanService
from vyrtuous.service.infractions.flag_service import FlagService
from vyrtuous.service.infractions.text_mute_service import TextMuteService
from vyrtuous.service.infractions.voice_mute_service import VoiceMuteService
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.mgmt.alias_service import AliasService
from vyrtuous.service.roles.administrator_service import AdministratorService
from vyrtuous.service.roles.coordinator_service import CoordinatorService
from vyrtuous.service.roles.developer_service import DeveloperService
from vyrtuous.service.roles.moderator_service import (
    ModeratorService,
    moderator_predicator,
)
from vyrtuous.service.roles.vegan_service import VeganService
from vyrtuous.service.rooms.stage_service import StageService
from vyrtuous.service.rooms.temporary_room_service import TemporaryRoomService
from vyrtuous.service.state_service import StateService
from vyrtuous.utils.dir_to_classes import dir_to_classes
from vyrtuous.utils.home import at_home


class ModeratorTextCommands(commands.Cog):
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
        target: str = commands.parameter(
            description="Specify one of: `all`, " "channel ID/mention or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=str(target))
        pages = await AdministratorService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural="Administrators", pages=pages, state=state)

    @commands.command(name="bans", description="List bans.")
    @moderator_predicator()
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await BanService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural="Bans", pages=pages, state=state)

    @commands.command(name="cmds", help="List aliases.")
    @moderator_predicator()
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await AliasService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural="Aliases", pages=pages, state=state)

    @commands.command(name="coords", help="Lists coords.")
    @moderator_predicator()
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: `all`, channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await CoordinatorService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural="Coordinators", pages=pages, state=state)

    @commands.command(name="del", help="Delete message.")
    @moderator_predicator()
    @skip_help_discovery()
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        message: MessageSnowflake = commands.parameter(description="Message snowflake"),
    ):
        state = StateService(ctx=ctx)
        for channel_obj in ctx.guild.channels:
            msg = await channel_obj.fetch_message(message)
        else:
            return await state.end(warning=f"Message `{message}` not found.")
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())
        return await state.end(success=f"Message `{message}` deleted successfully.")

    @commands.command(name="devs", help="List devs.")
    @moderator_predicator()
    @skip_help_discovery()
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        *,
        target: str = commands.parameter(
            description="'all', a specific server or user mention/ID"
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await DeveloperService.build_pages(object_dict=object_dict)
        await StateService.send_pages(plural="Developers", pages=pages, state=state)

    @commands.command(name="flags", help="List flags.")
    @moderator_predicator()
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await FlagService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural="Flags", pages=pages, state=state)

    @commands.command(name="ls", help="List new vegans.")
    @moderator_predicator()
    @skip_help_discovery()
    async def list_new_vegans_text_command(
        self,
        ctx: commands.Context,
        *,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await VeganService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural="Vegans", pages=pages, state=state)

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
        snowflake_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        do = DiscordObject(ctx=ctx)
        channel_dict = await do.determine_from_target(target=channel)
        msg = await TemporaryRoomService.migrate_temporary_room(
            channel_dict=channel_dict,
            old_name=old_name,
            snowflake_kwargs=snowflake_kwargs,
        )
        return await state.end(success=msg)

    @commands.command(name="mods", help="Lists mods.")
    @moderator_predicator()
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await ModeratorService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural="Moderators", pages=pages, state=state)

    @commands.command(name="mutes", help="List mutes.")
    @moderator_predicator()
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await VoiceMuteService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural="Voice Mutes", pages=pages, state=state)

    @commands.command(name="mstage", help="Stage mute/unmute.")
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
        snowflake_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        do = DiscordObject(ctx=ctx)
        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        msg = await StageService.toggle_stage_mute(
            channel_dict=channel_dict,
            member_dict=member_dict,
            snowflake_kwargs=snowflake_kwargs,
        )
        await state.end(success=msg)

    @commands.command(name="roleid", help="Get role by name.")
    @moderator_predicator()
    @skip_help_discovery()
    async def get_role_id_text_command(self, ctx: commands.Context, *, role_name: str):
        state = StateService(ctx=ctx)
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            return await state.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await state.end(
                warning=f"No role named `{role_name}` found in this server."
            )

    @commands.command(name="summary", description="Moderation summary.")
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
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/infractions")
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        member_dict = await do.determine_from_target(target=member)
        for obj in dir_to_classes(dir_paths=dir_paths):
            if "member" in obj.SCOPES:
                object_pages = await obj.build_pages(
                    object_dict=member_dict, is_at_home=is_at_home
                )
                if object_pages:
                    pages.extend(object_pages)
        await StateService.send_pages(plural="infractions", pages=pages, state=state)

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
        channel_dict = await do.determine_from_target(target=channel)
        pages = await StageService.survey(
            channel_dict=channel_dict, guild_snowflake=ctx.guild.id
        )
        await StateService.send_pages(plural="Stage Roles", pages=pages, state=state)

    @commands.command(name="tmutes", help="List text-mutes.")
    @moderator_predicator()
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await TextMuteService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural="Text Mutes", pages=pages, state=state)


async def setup(bot: DiscordBot):
    cog = ModeratorTextCommands(bot)
    await bot.add_cog(cog)
