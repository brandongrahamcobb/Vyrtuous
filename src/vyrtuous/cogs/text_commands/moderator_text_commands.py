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


from discord.ext import commands
import discord

from vyrtuous.cogs.help_command import skip_help_discovery
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.actions.ban import Ban
from vyrtuous.db.actions.flag import Flag
from vyrtuous.db.actions.text_mute import TextMute
from vyrtuous.db.actions.voice_mute import VoiceMute
from vyrtuous.db.roles.administrator import Administrator
from vyrtuous.db.roles.coordinator import Coordinator
from vyrtuous.db.roles.developer import Developer
from vyrtuous.db.roles.moderator import (
    Moderator,
    moderator_predicator,
)
from vyrtuous.db.roles.vegan import Vegan
from vyrtuous.db.rooms.stage import Stage
from vyrtuous.db.rooms.temporary_room import TemporaryRoom
from vyrtuous.fields.snowflake import MessageSnowflake
from vyrtuous.fields.snowflake import (
    ChannelSnowflake,
    MemberSnowflake,
)
from vyrtuous.utils.home import at_home
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import DiscordObject
from vyrtuous.utils.dir_to_classes import dir_to_classes


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
        pages = await Administrator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            plural=Administrator.PLURAL, pages=pages, state=state
        )

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
        pages = await Ban.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Ban.PLURAL, pages=pages, state=state)

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
        pages = await Alias.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Alias.PLURAL, pages=pages, state=state)

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
        pages = await Coordinator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            plural=Coordinator.PLURAL, pages=pages, state=state
        )

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
        pages = await Developer.build_pages(object_dict=object_dict)
        await StateService.send_pages(plural=Developer.PLURAL, pages=pages, state=state)

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
        pages = await Flag.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Flag.PLURAL, pages=pages, state=state)

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
        pages = await Vegan.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Vegan.PLURAL, pages=pages, state=state)

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
            "channel_snowflake": ctx.channel.id,
            "guild_snowflake": ctx.guild.id,
            "member_snowflake": ctx.author.id,
        }
        do = DiscordObject(ctx=ctx)
        channel_dict = await do.determine_from_target(target=channel)
        msg = await TemporaryRoom.migrate_temporary_room(
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
        pages = await Moderator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=Moderator.PLURAL, pages=pages, state=state)

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
        pages = await VoiceMute.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=VoiceMute.PLURAL, pages=pages, state=state)

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
            "channel_snowflake": ctx.channel.id,
            "guild_snowflake": ctx.guild.id,
            "member_snowflake": ctx.author.id,
        }
        do = DiscordObject(ctx=ctx)
        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        msg = await Stage.toggle_stage_mute(
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
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/actions")
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
        pages = await Stage.survey(
            channel_dict=channel_dict, guild_snowflake=ctx.guild.id
        )
        await StateService.send_pages(plural=Stage.PLURAL, pages=pages, state=state)

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
        pages = await TextMute.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=TextMute.PLURAL, pages=pages, state=state)


async def setup(bot: DiscordBot):
    cog = ModeratorTextCommands(bot)
    await bot.add_cog(cog)
