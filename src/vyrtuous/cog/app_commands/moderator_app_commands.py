"""!/bin/python3
moderator_app_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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
from discord import app_commands
from discord.ext import commands

from vyrtuous.administrator.administrator_service import AdministratorService
from vyrtuous.alias.alias_service import AliasService
from vyrtuous.ban.ban import Ban
from vyrtuous.ban.ban_service import BanService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.field.snowflake import (AppChannelSnowflake, AppMemberSnowflake,
                                      AppMessageSnowflake)
from vyrtuous.flag.flag_service import FlagService
from vyrtuous.modal.duration_modal import DurationModal
from vyrtuous.modal.reason_modal import ReasonModal
from vyrtuous.moderator.moderator import Moderator
from vyrtuous.moderator.moderator_service import (ModeratorService,
                                                  moderator_predicator)
from vyrtuous.stage_room.stage_service import StageService
from vyrtuous.temporary_room.temporary_room_service import TemporaryRoomService
from vyrtuous.text_mute.text_mute import TextMute
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.utils.dir_to_classes import dir_to_classes
from vyrtuous.utils.discord_object_service import DiscordObject
from vyrtuous.utils.home import at_home
from vyrtuous.utils.message_service import MessageService
from vyrtuous.utils.state_service import StateService
from vyrtuous.utils.view_context import ViewContext
from vyrtuous.vegan.vegan_service import VeganService
from vyrtuous.view.data_view import DataView
from vyrtuous.view.infraction_view import InfractionView
from vyrtuous.view.modify_infraction_view import ModifyInfractionView
from vyrtuous.voice_mute.voice_mute import VoiceMute
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService


class ModeratorAppCommands(commands.Cog):
    ROLE = Moderator

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.message_service = MessageService(self.bot)

    @app_commands.command(name="admins", description="Lists admins.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_administrators_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.guild.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await AdministratorService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Administrators", pages=pages, state=state)

    @app_commands.command(name="bans", description="List bans.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_bans_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await BanService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Bans", pages=pages, state=state)

    @app_commands.command(name="cmds", description="List aliases.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_commands_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await AliasService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Aliases", pages=pages, state=state)

    @app_commands.command(name="coords", description="Lists coords.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_coordinators_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await CoordinatorService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Coordinators", pages=pages, state=state)

    @app_commands.command(name="data", description="Create a chart.")
    @moderator_predicator()
    async def create_data_app_command(self, interaction: discord.Interaction):
        state = StateService(interaction=interaction)
        default_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        view = DataView(interaction=interaction, state=state)
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel, duration and infraction",
            view=view,
            ephemeral=True,
        )

    @app_commands.command(name="del", description="Delete message.")
    @app_commands.describe(message="Message ID")
    @moderator_predicator()
    async def delete_message_app_command(
        self, interaction: discord.Interaction, message: AppMessageSnowflake
    ):
        state = StateService(interaction=interaction)
        try:
            msg = await interaction.channel.fetch_message(message)
        except discord.NotFound:
            return await state.end(warning=f"Message `{message}` not found.")
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())
        return await state.end(success=f"Message `{message}` deleted successfully.")

    @app_commands.command(name="duration", description="Modify a duration.")
    @app_commands.describe(member="The ID or mention of the member.")
    @moderator_predicator()
    async def change_moderation_duration_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        ctx = ViewContext(interaction=interaction)
        await ctx.setup(target_member_snowflake=member_dict.get("id", None))
        view = ModifyInfractionView(
            ctx=ctx,
            modal=DurationModal,
            state=state,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and an infraction", view=view, ephemeral=True
        )

    @app_commands.command(name="flags", description="List flags.")
    @moderator_predicator()
    async def list_flags_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await FlagService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Flags", pages=pages, state=state)

    @app_commands.command(name="ls", description="List new vegans.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_new_vegans_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.guild.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await VeganService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Vegans", pages=pages, state=state)

    @app_commands.command(
        name="migrate", description="Migrate a temporary room to a new channel."
    )
    @app_commands.describe(
        old_name="Old temporary room name", channel="New channel to migrate to"
    )
    @moderator_predicator()
    async def migrate_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        old_name: str,
        channel: AppChannelSnowflake,
    ):
        state = StateService(interaction=interaction)
        default_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        await TemporaryRoomService.migrate_temporary_room(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            old_name=old_name,
        )
        return await state.end(
            success=f"Temporary room `{old_name}` migrated to {channel_dict.get('mention', None)}."
        )

    @app_commands.command(name="mods", description="Lists mods.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_moderators_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await ModeratorService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Moderators", pages=pages, state=state)

    @app_commands.command(name="mutes", description="List mutes.")
    @moderator_predicator()
    async def list_mutes_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await VoiceMuteService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Voice Mutes", pages=pages, state=state)

    @app_commands.command(name="mstage", description="Toggle stage mute/unmute.")
    @app_commands.describe(member="Tag a member or include their ID")
    @moderator_predicator()
    async def stage_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake,
    ):
        state = StateService(interaction=interaction)
        default_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        do = DiscordObject(interaction=interaction)
        channel = channel or int(interaction.channel.id)
        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        msg = await StageService.toggle_stage_mute(
            channel_dict=channel_dict,
            default_kwargs=default_kwargs,
            member_dict=member_dict,
        )
        await state.end(success=msg)

    @app_commands.command(name="reason", description="Modify a reason.")
    @app_commands.describe(member="The ID or mention of the member.")
    @moderator_predicator()
    async def change_moderation_reason_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        default_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        member_dict = await do.determine_from_target(target=member)
        ctx = ViewContext(interaction=interaction)
        await ctx.setup(target_member_snowflake=member_dict.get("id", None))
        view = ModifyInfractionView(
            ctx=ctx,
            modal=ReasonModal,
            state=state,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and an infraction", view=view, ephemeral=True
        )

    @app_commands.command(name="summary", description="Moderation summary.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
    )
    @moderator_predicator()
    async def list_moderation_summary_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
    ):
        pages = []
        dir_paths = []
        dir_paths.append(Path("/app/vyrtuous/db/infractions"))

        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        member_dict = await do.determine_from_target(target=member)
        for obj in dir_to_classes(dir_paths=dir_paths, parent=AliasService):
            object_pages = await obj.service.build_pages(
                object_dict=member_dict, is_at_home=is_at_home
            )
            pages.extend(object_pages)
        await StateService.send_pages(title="infractions", pages=pages, state=state)

    @app_commands.command(name="survey", description="Survey stage members.")
    @app_commands.describe(channel="Tag a voice/stage channel")
    @moderator_predicator()
    async def stage_survey_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        channel = channel or int(interaction.channel.id)
        channel_dict = await do.determine_from_target(target=channel)
        pages = await StageService.survey(
            channel_dict=channel_dict, guild_snowflake=interaction.guild.id
        )
        await StateService.send_pages(title="Stage Roles", pages=pages, state=state)

    @app_commands.command(name="tmutes", description="List text-mutes.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_text_mutes_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await TextMuteService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Text Mutes", pages=pages, state=state)

    @app_commands.command(name="vban", description="Create a ban.")
    @app_commands.describe(member="The ID or mention of the member.")
    @moderator_predicator()
    async def create_ban_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        ctx = ViewContext(interaction=interaction)
        ctx.record = Ban
        await ctx.setup(target_member_snowflake=member_dict.get("id", None))
        view = InfractionView(
            ctx=ctx,
            state=state,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a duration", view=view, ephemeral=True
        )

    @app_commands.command(name="vmute", description="Create a mute.")
    @app_commands.describe(member="The ID or mention of the member.")
    @moderator_predicator()
    async def create_voice_mute_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        ctx = ViewContext(interaction=interaction)
        ctx.record = VoiceMute
        await ctx.setup(target_member_snowflake=member_dict.get("id", None))
        view = InfractionView(
            ctx=ctx,
            state=state,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a duration", view=view, ephemeral=True
        )

    @app_commands.command(name="vtmute", description="Create a text-mute.")
    @app_commands.describe(member="The ID or mention of the member.")
    @moderator_predicator()
    async def create_text_mute_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        ctx = ViewContext(interaction=interaction)
        ctx.record = TextMute
        await ctx.setup(target_member_snowflake=member_dict.get("id", None))
        view = InfractionView(
            ctx=ctx,
            state=state,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a duration", view=view, ephemeral=True
        )


async def setup(bot: DiscordBot):
    cog = ModeratorAppCommands(bot)
    await bot.add_cog(cog)
