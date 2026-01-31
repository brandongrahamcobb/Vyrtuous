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

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.infractions.ban import Ban
from vyrtuous.db.infractions.flag import Flag
from vyrtuous.db.infractions.text_mute import TextMute
from vyrtuous.db.infractions.voice_mute import VoiceMute
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
from vyrtuous.fields.snowflake import AppMessageSnowflake
from vyrtuous.fields.snowflake import (
    AppChannelSnowflake,
    AppMemberSnowflake,
)
from vyrtuous.utils.home import at_home
from vyrtuous.utils.check import (
    has_equal_or_lower_role_wrapper,
    HasEqualOrLowerRole,
)
from vyrtuous.utils.moderation_view import ModerationView
from vyrtuous.utils.duration_modal import DurationModal
from vyrtuous.utils.reason_modal import ReasonModal
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import DiscordObject
from vyrtuous.utils.dir_to_classes import dir_to_classes


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
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Administrator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            plural=Administrator.PLURAL, pages=pages, state=state
        )

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
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Ban.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Ban.PLURAL, pages=pages, state=state)

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
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Alias.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Alias.PLURAL, pages=pages, state=state)

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
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Coordinator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            plural=Coordinator.PLURAL, pages=pages, state=state
        )

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
        object_dict = await do.determine_from_target(target=target)
        pages = await Developer.build_pages(object_dict=object_dict)
        await StateService.send_pages(plural=Developer.PLURAL, pages=pages, state=state)

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
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Flag.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Flag.PLURAL, pages=pages, state=state)

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
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Vegan.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Vegan.PLURAL, pages=pages, state=state)

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
        await TemporaryRoom.migrate_temporary_room(
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
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Moderator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=Moderator.PLURAL, pages=pages, state=state)

    @app_commands.command(name="mutes", description="List mutes.")
    @moderator_predicator()
    async def list_mutes_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await VoiceMute.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=VoiceMute.PLURAL, pages=pages, state=state)

    @app_commands.command(name="mstage", description="Stage mute/unmute.")
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
        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        msg = await Stage.toggle_stage_mute(
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
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/infractions")
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        member_dict = await do.determine_from_target(target=member)
        for obj in dir_to_classes(dir_paths=dir_paths):
            if "member" in obj.SCOPES:
                object_pages = await obj.build_pages(
                    object_dict=member_dict, is_at_home=is_at_home
                )
                pages.extend(object_pages)
        await StateService.send_pages(plural="infractions", pages=pages, state=state)

    @app_commands.command(name="survey", description="Survey stage members.")
    @app_commands.describe(channel="Tag a voice/stage channel")
    @moderator_predicator()
    async def stage_survey_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        pages = await Stage.survey(
            channel_dict=channel_dict, guild_snowflake=interaction.guild.id
        )
        await StateService.send_pages(plural=Stage.PLURAL, pages=pages, state=state)

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
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await TextMute.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=TextMute.PLURAL, pages=pages, state=state)


async def setup(bot: DiscordBot):
    cog = ModeratorAppCommands(bot)
    await bot.add_cog(cog)
