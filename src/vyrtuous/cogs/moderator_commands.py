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
from typing import Optional

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.actions.ban import Ban
from vyrtuous.db.mgmt.cap import Cap
from vyrtuous.db.actions.flag import Flag
from vyrtuous.db.actions.server_mute import ServerMute
from vyrtuous.db.actions.text_mute import TextMute
from vyrtuous.db.actions.voice_mute import VoiceMute
from vyrtuous.db.roles.administrator import Administrator, is_administrator
from vyrtuous.db.roles.coordinator import Coordinator, is_coordinator
from vyrtuous.db.roles.developer import Developer, is_developer
from vyrtuous.db.roles.guild_owner import is_guild_owner
from vyrtuous.db.roles.moderator import (
    Moderator,
    is_moderator,
    moderator_predicator,
)
from vyrtuous.db.roles.sysadmin import is_sysadmin
from vyrtuous.db.roles.vegan import Vegan
from vyrtuous.db.rooms.stage import Stage
from vyrtuous.db.rooms.temporary_room import TemporaryRoom
from vyrtuous.fields.snowflake import AppMessageSnowflake, MessageSnowflake
from vyrtuous.fields.snowflake import (
    AppChannelSnowflake,
    AppMemberSnowflake,
    ChannelSnowflake,
    MemberSnowflake,
)
from vyrtuous.utils.home import at_home
from vyrtuous.utils.check import (
    check,
    has_equal_or_lower_role_wrapper,
    HasEqualOrLowerRole,
)
from vyrtuous.utils.logger import logger
from vyrtuous.utils.moderation_view import ModerationView
from vyrtuous.utils.duration_modal import DurationModal
from vyrtuous.utils.reason_modal import ReasonModal
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import DiscordObject
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.dir_to_classes import dir_to_classes


class ModeratorCommands(commands.Cog):
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
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Administrator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            plural=Administrator.PLURAL, pages=pages, state=state
        )

    # DONE
    @commands.command(name="admins", help="Lists admins.")
    @moderator_predicator()
    async def list_administrators_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            description="Specify one of: `all`, " "channel ID/mention or server ID.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=str(target))
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
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Ban.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Ban.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="bans", description="List bans.")
    @moderator_predicator()
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention or server ID.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await Ban.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Ban.PLURAL, pages=pages, state=state)

    # DONE
    @app_commands.command(name="cmds", description="List aliases.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_commands_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Alias.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Alias.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="cmds", help="List aliases.")
    @moderator_predicator()
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await Alias.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Alias.PLURAL, pages=pages, state=state)

    @app_commands.command(name="coords", description="Lists coords.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_coordinators_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Coordinator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            plural=Coordinator.PLURAL, pages=pages, state=state
        )

    # DONE
    @commands.command(name="coords", help="Lists coords.")
    @moderator_predicator()
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: `all`, channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await Coordinator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            plural=Coordinator.PLURAL, pages=pages, state=state
        )

    # DONE
    @app_commands.command(name="del", description="Delete message.")
    @app_commands.describe(message="Message ID")
    @moderator_predicator()
    async def delete_message_app_command(
        self, interaction: discord.Interaction, message: AppMessageSnowflake
    ):
        state = StateService(source=interaction)

        for channel_obj in interaction.guild.channels:
            msg = await channel_obj.fetch_message(message)
        else:
            return await state.end(warning=f"Message `{message}` not found.")
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())

        return await state.end(success=f"Message `{message}` deleted successfully.")

    # DONE
    @commands.command(name="del", help="Delete message.")
    @moderator_predicator()
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        message: MessageSnowflake = commands.parameter(description="Message snowflake"),
    ):
        state = StateService(source=ctx)

        for channel_obj in ctx.guild.channels:
            msg = await channel_obj.fetch_message(message)
        else:
            return await state.end(warning=f"Message `{message}` not found.")
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())

        return await state.end(success=f"Message `{message}` deleted successfully.")

    # DONE
    @app_commands.command(name="devs", description="List devs.")
    @app_commands.describe(target="Specify one of: 'all', or server ID.")
    @moderator_predicator()
    async def list_developers_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Developer.build_pages(object_dict=object_dict)
        await StateService.send_pages(plural=Developer.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="devs", help="List devs.")
    @moderator_predicator()
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            description="'all', a specific server or user mention/ID"
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
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
            state = StateService(source=interaction)
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

    # DONE
    @app_commands.command(name="flags", description="List flags.")
    @moderator_predicator()
    async def list_flags_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Flag.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Flag.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="flags", help="List flags.")
    @moderator_predicator()
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await Flag.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Flag.PLURAL, pages=pages, state=state)

    # DONE
    @app_commands.command(name="ls", description="List new vegans.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_new_vegans_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Vegan.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Vegan.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="ls", help="List new vegans.")
    @moderator_predicator()
    async def list_new_vegans_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await Vegan.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Vegan.PLURAL, pages=pages, state=state)

    # DONE
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
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        old_room = await TemporaryRoom.select(
            guild_snowflake=interaction.guild.id, room_name=old_name
        )
        if old_room:
            channel_dict = await do.determine_from_target(target=channel)
            is_owner = old_room.member_snowflake == interaction.user.id
            await check(source=interaction, lowest_role="Administrator") or is_owner
            set_kwargs = {"channel_snowflake": channel_dict.get("id", None)}
            temp_where_kwargs = {
                "channel_snowflake": old_room.channel_snowflake,
                "guild_snowflake": interaction.guild.id,
                "room_name": channel_dict.get("name", None),
            }
            where_kwargs = {
                "channel_snowflake": old_room.channel_snowflake,
                "guild_snowflake": interaction.guild.id,
            }
            kwargs = {
                "set_kwargs": set_kwargs,
                "where_kwargs": where_kwargs,
            }
            await TemporaryRoom.update(
                set_kwargs=set_kwargs,
                where_kwargs=temp_where_kwargs,
            )
            await Alias.update(**kwargs)
            await Ban.update(**kwargs)
            await Cap.update(**kwargs)
            await Coordinator.update(**kwargs)
            await Flag.update(**kwargs)
            await Moderator.update(**kwargs)
            await Stage.update(**kwargs)
            await TextMute.update(**kwargs)
            await Vegan.update(**kwargs)
            await VoiceMute.update(**kwargs)
            return await state.end(
                success=f"Temporary room `{old_name}` migrated to {channel_dict.get("mention", None)}."
            )
        return await state.end(
            warning=f"No temporary rooms found called `{old_name}` in {interaction.guild.name}."
        )

    # DONE
    @commands.command(
        name="migrate",
        help="Migrate a temporary room to a new channel by snowflake.",
        hidden=True,
    )
    @moderator_predicator()
    async def migrate_temp_room_text_command(
        self,
        ctx: commands.Context,
        old_name: str = commands.parameter(description="Provide a channel name"),
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        old_room = await TemporaryRoom.select(
            guild_snowflake=ctx.guild.id, room_name=old_name
        )
        if old_room:
            channel_dict = await do.determine_from_target(target=channel)
            is_owner = old_room.member_snowflake == ctx.author.id
            await check(source=ctx, lowest_role="Administrator") or is_owner
            set_kwargs = {"channel_snowflake": channel_dict.get("id", None)}
            temp_where_kwargs = {
                "channel_snowflake": old_room.channel_snowflake,
                "guild_snowflake": ctx.guild.id,
                "room_name": channel_dict.get("name", None),
            }
            where_kwargs = {
                "channel_snowflake": old_room.channel_snowflake,
                "guild_snowflake": ctx.guild.id,
            }
            kwargs = {
                "set_kwargs": set_kwargs,
                "where_kwargs": where_kwargs,
            }
            await TemporaryRoom.update(
                set_kwargs=set_kwargs,
                where_kwargs=temp_where_kwargs,
            )
            await Alias.update(**kwargs)
            await Ban.update(**kwargs)
            await Cap.update(**kwargs)
            await Coordinator.update(**kwargs)
            await Flag.update(**kwargs)
            await Moderator.update(**kwargs)
            await Stage.update(**kwargs)
            await TextMute.update(**kwargs)
            await Vegan.update(**kwargs)
            return await state.end(
                success=f"Temporary room `{old_name}` migrated to {channel_dict.get("mention", None)}."
            )

        return await state.end(
            warning=f"No temporary rooms found called `{old_name}` in {ctx.guild.name}."
        )

    # DONE
    @app_commands.command(name="mods", description="Lists mods.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_moderators_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Moderator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=Moderator.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="mods", help="Lists mods.")
    @moderator_predicator()
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await Moderator.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=Moderator.PLURAL, pages=pages, state=state)

    # DONE
    @app_commands.command(name="mutes", description="List mutes.")
    @moderator_predicator()
    async def list_mutes_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await VoiceMute.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=VoiceMute.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="mutes", help="List mutes.")
    @moderator_predicator()
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await VoiceMute.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=VoiceMute.PLURAL, pages=pages, state=state)

    # DONE
    @app_commands.command(name="mstage", description="Stage mute/unmute.")
    @app_commands.describe(member="Tag a member or include their ID")
    async def stage_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake,
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        await has_equal_or_lower_role_wrapper(
            source=interaction,
            member_snowflake=member_dict.get("id", None),
            sender_snowflake=interaction.user.id,
        )
        kwargs = channel_dict.get("columns", None)

        stage = await Stage.select(**kwargs)
        if not stage:
            return await state.end(
                warning=f"No active stage found in {channel_dict.get("mention", None)}."
            )
        try:
            await member_dict.get("object", None).edit(
                mute=not member_dict.get("object", None).voice.mute
            )
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @commands.command(name="mstage", help="Stage mute/unmute.")
    @moderator_predicator()
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
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        await has_equal_or_lower_role_wrapper(
            source=ctx,
            member_snowflake=member_dict.get("id", None),
            sender_snowflake=ctx.author.id,
        )
        kwargs = channel_dict.get("columns", None)

        stage = await Stage.select(**kwargs)
        if not stage:
            return await state.end(
                warning=f"No active stage found in {channel_dict.get("mention", None)}."
            )
        try:
            await member_dict.get("object", None).edit(
                mute=not member_dict.get("object", None).voice.mute
            )
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())

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
            state = StateService(source=interaction)
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

    # DONE
    @app_commands.command(name="roleid", description="Get role by name.")
    @app_commands.describe(role_name="The name of the role to look up")
    @moderator_predicator()
    async def get_role_id_app_command(
        self, interaction: discord.Interaction, role_name: str
    ):
        state = StateService(source=interaction)
        role = discord.utils.get(interaction.guild.roles, name=role_name)
        if role:
            return await state.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await state.end(
                warning=f"No role named `{role_name}` found in this server."
            )

    # DONE
    @commands.command(name="roleid", help="Get role by name.")
    @moderator_predicator()
    async def get_role_id_text_command(self, ctx: commands.Context, *, role_name: str):
        state = StateService(source=ctx)
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            return await state.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await state.end(
                warning=f"No role named `{role_name}` found in this server."
            )

    @app_commands.command(name="smutes", description="List mutes.")
    @moderator_predicator()
    async def list_server_mutes_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await ServerMute.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            plural=ServerMute.PLURAL, pages=pages, state=state
        )

    # DONE
    @commands.command(name="smutes", help="List mutes.")
    @moderator_predicator()
    async def list_server_mutes_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await ServerMute.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            plural=ServerMute.PLURAL, pages=pages, state=state
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
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/actions")
        state = StateService(source=interaction)
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
        state = StateService(source=ctx)
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

    @app_commands.command(name="survey", description="Get all.")
    @app_commands.describe(channel="Tag a voice/stage channel")
    @moderator_predicator()
    async def stage_survey_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        chunk_size, pages = 7, []
        sysadmins = developers = guild_owners = administrators = coordinators = (
            moderators
        ) = []

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)

        for member in channel_dict.get("object", None).members:
            try:
                if is_sysadmin(member.id):
                    sysadmins.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_developer(member.id):
                    developers.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_guild_owner(interaction.guild.id, member.id):
                    guild_owners.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_administrator(interaction.guild.id, member.id):
                    administrators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_coordinator(
                    channel_dict.get("id", None), interaction.guild.id, member.id
                ):
                    coordinators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_moderator(
                    channel_dict.get("id", None), interaction.guild.id, member.id
                ):
                    moderators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
        sysadmins_chunks = [
            sysadmins[i : i + chunk_size] for i in range(0, len(sysadmins), chunk_size)
        ]
        guild_owners_chunks = [
            guild_owners[i : i + chunk_size]
            for i in range(0, len(guild_owners), chunk_size)
        ]
        developers_chunks = [
            developers[i : i + chunk_size]
            for i in range(0, len(developers), chunk_size)
        ]
        administrators_chunks = [
            administrators[i : i + chunk_size]
            for i in range(0, len(administrators), chunk_size)
        ]
        coordinators_chunks = [
            coordinators[i : i + chunk_size]
            for i in range(0, len(coordinators), chunk_size)
        ]
        moderators_chunks = [
            moderators[i : i + chunk_size]
            for i in range(0, len(moderators), chunk_size)
        ]
        roles_chunks = [
            ("Sysadmins", sysadmins, sysadmins_chunks),
            ("Developers", developers, developers_chunks),
            ("Guild Owners", guild_owners, guild_owners_chunks),
            ("Administrators", administrators, administrators_chunks),
            ("Coordinators", coordinators, coordinators_chunks),
            ("Moderators", moderators, moderators_chunks),
        ]
        max_pages = max(len(c[2]) for c in roles_chunks)
        for page in range(max_pages):
            embed = discord.Embed(
                title=f"{get_random_emoji()} Survey results for {channel_dict.get('name', None)}",
                description=f"Total surveyed: {len(channel_dict.get('object', None).members)}",
                color=discord.Color.blurple(),
            )
            for role_name, role_list, chunks in roles_chunks:
                chunk = chunks[page] if page < len(chunks) else []
                embed.add_field(
                    name=f"{role_name} ({len(chunk)}/{len(role_list)})",
                    value=", ".join(u.mention for u in chunk) if chunk else "*None*",
                    inline=False,
                )
            pages.append(embed)

        if pages:
            return await state.end(success=pages)
        else:
            return await state.end(warning="No permissions found.")

    # DONE
    @commands.command(name="survey", help="Get all.")
    @moderator_predicator()
    async def stage_survey_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
    ):
        chunk_size, pages = 7, []
        sysadmins = developers = guild_owners = administrators = coordinators = (
            moderators
        ) = []

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)
        for member in channel_dict.get("object", None).members:
            try:
                if is_sysadmin(member.id):
                    sysadmins.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_developer(member.id):
                    developers.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_guild_owner(ctx.guild.id, member.id):
                    guild_owners.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_administrator(ctx.guild.id, member.id):
                    administrators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_coordinator(
                    channel_dict.get("id", None), ctx.guild.id, member.id
                ):
                    coordinators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_moderator(
                    channel_dict.get("id", None), ctx.guild.id, member.id
                ):
                    moderators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
        sysadmins_chunks = [
            sysadmins[i : i + chunk_size] for i in range(0, len(sysadmins), chunk_size)
        ]
        guild_owners_chunks = [
            guild_owners[i : i + chunk_size]
            for i in range(0, len(guild_owners), chunk_size)
        ]
        developers_chunks = [
            developers[i : i + chunk_size]
            for i in range(0, len(developers), chunk_size)
        ]
        administrators_chunks = [
            administrators[i : i + chunk_size]
            for i in range(0, len(administrators), chunk_size)
        ]
        coordinators_chunks = [
            coordinators[i : i + chunk_size]
            for i in range(0, len(coordinators), chunk_size)
        ]
        moderators_chunks = [
            moderators[i : i + chunk_size]
            for i in range(0, len(moderators), chunk_size)
        ]
        roles_chunks = [
            ("Sysadmins", sysadmins, sysadmins_chunks),
            ("Developers", developers, developers_chunks),
            ("Guild Owners", guild_owners, guild_owners_chunks),
            ("Administrators", administrators, administrators_chunks),
            ("Coordinators", coordinators, coordinators_chunks),
            ("Moderators", moderators, moderators_chunks),
        ]
        max_pages = max(len(c[2]) for c in roles_chunks)
        for page in range(max_pages):
            embed = discord.Embed(
                title=f"{get_random_emoji()} Survey results for {channel_dict.get('name', None)}",
                description=f"Total surveyed: {len(channel_dict.get('object', None).members)}",
                color=discord.Color.blurple(),
            )
            for role_name, role_list, chunks in roles_chunks:
                chunk = chunks[page] if page < len(chunks) else []
                embed.add_field(
                    name=f"{role_name} ({len(chunk)}/{len(role_list)})",
                    value=", ".join(u.mention for u in chunk) if chunk else "*None*",
                    inline=False,
                )
            pages.append(embed)

        if pages:
            return await state.end(success=pages)
        else:
            return await state.end(warning="No permissions found.")

    # DONE
    @app_commands.command(name="tmutes", description="List text-mutes.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_text_mutes_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await TextMute.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=TextMute.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="tmutes", help="List text-mutes.")
    @moderator_predicator()
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        object_dict = await do.determine_from_target(target=target)
        pages = await TextMute.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=TextMute.PLURAL, pages=pages, state=state)


async def setup(bot: DiscordBot):
    cog = ModeratorCommands(bot)
    await bot.add_cog(cog)
