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
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.fields.snowflake import (
    AppChannelSnowflake,
    AppMemberSnowflake,
    AppMessageSnowflake,
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
from vyrtuous.utils.check import HasEqualOrLowerRole, has_equal_or_lower_role_wrapper
from vyrtuous.utils.dir_to_classes import dir_to_classes
from vyrtuous.utils.duration_modal import DurationModal
from vyrtuous.utils.home import at_home
from vyrtuous.utils.moderation_view import ModerationView
from vyrtuous.utils.reason_modal import ReasonModal


class ModeratorAppCommands(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.message_service = MessageService(self.bot)

    @app_commands.command(name="admins", description="Lists admins.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, " "or server ID."
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

    @app_commands.command(name="del", description="Delete message.")
    @app_commands.describe(message="Message ID")
    @moderator_predicator()
    async def delete_message_app_command(
        self, interaction: discord.Interaction, message: AppMessageSnowflake
    ):
        state = StateService(interaction=interaction)
        for channel_obj in interaction.guild.channels:
            msg = await channel_obj.fetch_message(message)
        else:
            return await state.end(warning=f"Message `{message}` not found.")
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())
        return await state.end(success=f"Message `{message}` deleted successfully.")

    @app_commands.command(name="devs", description="List devs.")
    @app_commands.describe(target="Specify one of: 'all', or server ID.")
    @moderator_predicator()
    async def list_developers_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or "all"
        object_dict = await do.determine_from_target(target=target)
        pages = await DeveloperService.build_pages(object_dict=object_dict)
        await StateService.send_pages(title="Developer", pages=pages, state=state)

    @app_commands.command(name="duration", description="Modify a duration.")
    @app_commands.describe(member="The ID or mention of the member.")
    @moderator_predicator()
    async def change_moderation_duration_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        try:
            await has_equal_or_lower_role_wrapper(
                source=interaction,
                member_snowflake=member_dict.get("id", None),
                sender_snowflake=interaction.user.id,
            )
        except HasEqualOrLowerRole as e:
            state = StateService(interaction=interaction)
            return await state.end(warning=str(e).capitalize())
        view = ModerationView(
            interaction=interaction,
            member_snowflake=member_dict.get("id", None),
            modal=DurationModal,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a category", view=view, ephemeral=True
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
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        await TemporaryRoomService.migrate_temporary_room(
            channel_dict=channel_dict,
            old_name=old_name,
            snowflake_kwargs=snowflake_kwargs,
        )
        return await state.end(
            success=f"Temporary room `{old_name}` migrated to {channel_dict.get("mention", None)}."
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
        snowflake_kwargs = {
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
            member_dict=member_dict,
            snowflake_kwargs=snowflake_kwargs,
        )
        await state.end(success=msg)

    @app_commands.command(name="reason", description="Modify a reason.")
    @app_commands.describe(member="The ID or mention of the member.")
    @moderator_predicator()
    async def change_moderation_reason_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        try:
            await has_equal_or_lower_role_wrapper(
                source=interaction,
                member_snowflake=member_dict.get("id", None),
                sender_snowflake=interaction.user.id,
            )
        except HasEqualOrLowerRole as e:
            state = StateService(interaction=interaction)
            return await state.end(warning=str(e).capitalize())
        view = ModerationView(
            interaction=interaction,
            member_snowflake=member_dict.get("id", None),
            modal=ReasonModal,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a category", view=view, ephemeral=True
        )

    @app_commands.command(name="roleid", description="Get role by name.")
    @app_commands.describe(role_name="The name of the role to look up")
    @moderator_predicator()
    async def get_role_id_app_command(
        self, interaction: discord.Interaction, role_name: str
    ):
        state = StateService(interaction=interaction)
        role = discord.utils.get(interaction.guild.roles, name=role_name)
        if role:
            return await state.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await state.end(
                warning=f"No role named `{role_name}` found in this server."
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
        dir_paths.append(Path(__file__).resolve().parents[2] / "service/infractions")
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        member_dict = await do.determine_from_target(target=member)
        for obj in dir_to_classes(dir_paths=dir_paths):
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


async def setup(bot: DiscordBot):
    cog = ModeratorAppCommands(bot)
    await bot.add_cog(cog)
