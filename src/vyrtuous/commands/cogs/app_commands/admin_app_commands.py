"""!/bin/python3

admin_app_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.clear_service import ClearService
from vyrtuous.commands.discord_object_service import DiscordObject
from vyrtuous.commands.fields.category import AppCategory
from vyrtuous.commands.fields.snowflake import (
    AppChannelSnowflake,
    AppMemberSnowflake,
    AppRoleSnowflake,
)
from vyrtuous.commands.home import at_home
from vyrtuous.commands.messaging.cancel_confirm import VerifyView
from vyrtuous.commands.messaging.message_service import MessageService
from vyrtuous.commands.messaging.state_service import StateService
from vyrtuous.commands.permission_check_service import PermissionCheckService
from vyrtuous.db.base.alias.alias_service import AliasService
from vyrtuous.db.infractions.smute.server_mute_service import ServerMuteService
from vyrtuous.db.infractions.vmute.voice_mute_service import VoiceMuteService
from vyrtuous.db.mgmt.cap.cap_service import CapService
from vyrtuous.db.mgmt.stream.stream_service import StreamService
from vyrtuous.db.roles.admin.administrator_service import (
    AdministratorRoleService,
    administrator_predicator,
)
from vyrtuous.db.roles.coord.coordinator_service import CoordinatorService
from vyrtuous.db.roles.permissions.check import check
from vyrtuous.db.rooms.stage.stage_service import StageService
from vyrtuous.db.rooms.temp.temporary_room_service import TemporaryRoomService
from vyrtuous.db.rooms.video.video_room_service import VideoRoomService
from vyrtuous.inc.helpers import PATH_LOG
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger


class AdminAppCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    @app_commands.command(name="alias", description="Alias creation.")
    @administrator_predicator()
    @app_commands.describe(
        alias_name="Alias/Pseudonym",
        category="Specify a category for a `ban`, `flag`, `role`, `tmute`, `vegan` or `vmute` action.",
        channel="Tag a channel or include the ID.",
        role="Tag a channel, role or include the ID.",
    )
    async def create_alias_app_command(
        self,
        interaction: discord.Interaction,
        category: AppCategory,
        alias_name: str,
        channel: AppChannelSnowflake,
        *,
        role: AppRoleSnowflake,
    ):
        role_dict = None
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        if category == "role":
            role_dict = await do.determine_from_target(target=role)
        msg = await AliasService.create_alias(
            alias_name=alias_name,
            category=category,
            channel_dict=channel_dict,
            role_dict=role_dict,
        )
        return await state.end(success=msg)

    @app_commands.command(name="aroles", description="Administrator roles.")
    @app_commands.describe(target="Specify one of: 'all', or server ID.")
    @administrator_predicator()
    async def list_administrator_roles_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.guild.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await AdministratorRoleService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            title="Administrator Roles", pages=pages, state=state
        )

    @app_commands.command(name="cap", description="Cap alias duration for mods.")
    @administrator_predicator()
    @app_commands.describe(
        channel="Tag a channel or include its ID.",
        category="One of: `mute`, `ban`, `tmute`",
        hours="(+|-)duration(m|h|d), 0=permanent, default=24h",
    )
    async def cap_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        category: AppCategory,
        hours: int,
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        msg = await CapService.toggle_cap(
            category=category, channel_dict=channel_dict, hours=hours
        )
        return await state.end(success=msg)

    @app_commands.command(name="caps", description="List caps.")
    @administrator_predicator()
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    async def list_caps_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await CapService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Caps", pages=pages, state=state)

    @app_commands.command(name="clear", description="Reset database.")
    @app_commands.describe(
        target="Specify 'all', tag a channel/guild/member or include its ID",
        category="Specify one of: `admin`, `alias`, `arole`, `all`, `ban`, `coord`, "
        "flag`, `mod`, `troom`, `tmute`, `stage`, `stream`, `vegan`, `vmute` or `vroom`.",
    )
    @administrator_predicator()
    async def clear_channel_access_app_command(
        self,
        interaction: discord.Interaction,
        target: str,
        category: AppCategory,
    ):
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        do = DiscordObject(interaction=interaction)
        object_dict = await do.determine_from_target(target=target)
        where_kwargs = object_dict.get("columns", None)
        view = VerifyView(
            category=str(category),
            mention=object_dict.get("mention", "All"),
            author_snowflake=interaction.user.id,
            **where_kwargs,
        )
        embed = view.build_embed()
        await interaction.response.send_message(embed=embed, view=view)
        await view.wait()
        state = StateService(interaction=interaction)
        msg = await ClearService.clear(
            category=category,
            object_dict=object_dict,
            snowflake_kwargs=snowflake_kwargs,
            where_kwargs=where_kwargs,
            target=target,
            view=view,
        )
        return await state.end(success=msg)

    @app_commands.command(name="coord", description="Grant/revoke coords.")
    @app_commands.describe(
        member="Tag a member or include their ID",
        channel="Tag a channel or include its ID.",
    )
    @administrator_predicator()
    async def toggle_coordinator_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake,
    ):
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)
        msg = await CoordinatorService.toggle_coordinator(
            channel_dict=channel_dict,
            member_dict=member_dict,
            snowflake_kwargs=snowflake_kwargs,
        )
        return await state.end(success=msg)

    @app_commands.command(
        name="debug", description="Shows the last `n` number of logging."
    )
    @administrator_predicator()
    async def debug_app_command(self, interaction: discord.Interaction, lines: int = 3):
        state = StateService(interaction=interaction)
        if lines <= 0:
            return await state.end(warning="Lines must be greater than 0")
        try:
            with open(PATH_LOG, "r") as f:
                content = f.readlines()[-lines:]
                content = [line.split(" - ", 3)[-1] for line in content]
        except FileNotFoundError:
            return await state.end(warning="Log file not found")
        output = "".join(content)
        if len(output) > 1900:
            output = output[-1900:]
        return await state.end(success=f"```log\n{output}\n```")

    @app_commands.command(name="pc", description="View permissions.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @administrator_predicator()
    async def list_permissions_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        is_at_home = at_home(source=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        object_dict = await do.determine_from_target(target=target)
        if target and str(target).lower() == "all":
            await check(snowflake_kwargs=snowflake_kwargs, lowest_role="Developer")
        if target and str(target).lower() == "all":
            await check(snowflake_kwargs=snowflake_kwargs, lowest_role="Developer")
        if target and str(target).lower() == "all":
            channel_objs = [
                channel_obj
                for guild in self.bot.guilds
                for channel_obj in guild.channels
            ]
        elif hasattr(object_dict.get("object", None), "channels"):
            channel_objs = object_dict.get("object", None).channels
        else:
            channel_objs = [object_dict.get("object", None)]
        pages = await PermissionCheckService.build_pages(
            channel_objs=channel_objs,
            is_at_home=is_at_home,
            snowflake_kwargs=snowflake_kwargs,
        )
        if pages:
            return await state.end(warning=pages)
        else:
            return await state.end(
                success=f"{self.bot.user.display_name} has all permissions for `{target}`."
            )

    @app_commands.command(name="rmute", description="Room mute (except yourself).")
    @app_commands.describe(channel="Tag a channel or include its ID.")
    @administrator_predicator()
    async def room_mute_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        reason: str = "No reason provided.",
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
        pages = await VoiceMuteService.room_mute(
            channel_dict=channel_dict,
            guild_snowflake=interaction.guild.id,
            reason=reason,
            snowflake_kwargs=snowflake_kwargs,
        )
        await StateService.send_pages(title="Voice Mutes", pages=pages, state=state)

    @app_commands.command(name="xrmute", description="Unmute all.")
    @app_commands.describe(channel="Tag a channel or include its ID.")
    @administrator_predicator()
    async def room_unmute_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        channel = channel or int(interaction.channel.id)
        channel_dict = await do.determine_from_target(target=channel)
        pages = await VoiceMuteService.room_unmute(
            channel_dict=channel_dict, guild_snowflake=interaction.guild.id
        )
        await StateService.send_pages(title="Voice Unmutes", pages=pages, state=state)

    @app_commands.command(name="rmv", description="VC move.")
    @app_commands.describe(
        source_channel="Tag the source channel or include its ID.",
        target_channel="Tag the target channel or include its ID.",
    )
    @administrator_predicator()
    async def room_move_all_app_command(
        self,
        interaction: discord.Interaction,
        source_channel: AppChannelSnowflake,
        target_channel: AppChannelSnowflake,
    ):
        failed, moved = [], []
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        source_channel_dict = await do.determine_from_target(target=source_channel)
        target_channel_dict = await do.determine_from_target(target=target_channel)
        for member in source_channel_dict.get("object", None).members:
            try:
                await member.move_to(target_channel_dict.get("object", None))
                moved.append(member)
            except discord.Forbidden as e:
                failed.append(member)
                logger.warning(
                    f"Unable to move member "
                    f"{member.display_name} ({member.id}) from channel "
                    f"{source_channel_dict.get('name', None)} ({source_channel}) to channel "
                    f"{target_channel_dict.get('name', None)} ({target_channel}) in guild "
                    f"{interaction.guild.name} ({interaction.guild.id}). "
                    f"{str(e).capitalize()}"
                )
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"Moved {source_channel_dict.get("mention", None)} to "
            f"{target_channel_dict.get("mention", None)}",
            color=discord.Color.green(),
        )
        if moved:
            embed.add_field(
                name=f"Successfully Moved (`{len(moved)}`)",
                value="\n".join(member.mention for member in moved),
                inline=False,
            )
        else:
            embed.add_field(name="Successfully Moved", value="None", inline=False)
        if failed:
            embed.add_field(
                name=f"Failed to Move ({len(failed)})",
                value="\n".join(member.mention for member in failed),
                inline=False,
            )
        embed.set_footer(
            text=f"Moved from {source_channel_dict.get("name", None)} "
            f"to {target_channel_dict.get("name", None)}"
        )
        return await state.end(success=embed)

    @app_commands.command(name="smute", description="Server mute/server unmute.")
    @app_commands.describe(
        member="Tag a member or include their ID",
        reason="Optional reason (required for 7 days or more)",
    )
    @administrator_predicator()
    async def toggle_server_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        reason: str = "No reason provided",
    ):
        state = StateService(interaction=interaction)
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        msg = await ServerMuteService.toggle_server_mute(
            member_dict=member_dict, reason=reason, snowflake_kwargs=snowflake_kwargs
        )
        return await state.end(success=msg)

    @app_commands.command(name="smutes", description="List mutes.")
    @administrator_predicator()
    async def list_server_mutes_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.guild.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await ServerMuteService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Server Mutes", pages=pages, state=state)

    @app_commands.command(name="stages", description="List stages.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @administrator_predicator()
    async def list_stages_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await StageService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Stages", pages=pages, state=state)

    @app_commands.command(
        name="temp", description="Toggle a temporary room and assign an owner."
    )
    @app_commands.describe(
        channel="Tag a channel or include its ID.",
    )
    @administrator_predicator()
    async def toggle_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        msg = await TemporaryRoomService.toggle_temporary_room(
            channel_dict=channel_dict
        )
        return await state.end(success=msg)

    @app_commands.command(
        name="temps", description="List temporary rooms with matching command aliases."
    )
    @administrator_predicator()
    async def list_temp_rooms_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await TemporaryRoomService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Temporary Rooms", pages=pages, state=state)

    @app_commands.command(name="stream", description="Setup streaming.")
    @app_commands.describe(
        channel="Tag a channel or include its ID where the messages will be sent.",
        action="create | modify | delete.",
        entry_type="all | channel | general.",
        snowflakes="Optional list of channel/member IDs to push events from.",
    )
    @administrator_predicator()
    async def modify_streaming_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        action: str,
        entry_type: str,
        snowflakes: str,
    ):
        channel_mentions, failed_snowflakes, resolved_channels = [], [], []
        state = StateService(interaction=interaction)
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        if snowflakes:
            for snowflake in snowflakes:
                try:
                    snowflake_dict = await do.determine_from_target(target=channel)
                    resolved_channels.append(snowflake_dict.get("id", None))
                    channel_mentions.append(snowflake_dict.get("mention", None))
                except Exception:
                    failed_snowflakes.append(snowflake)
                    continue
        pages = await StreamService.modify_stream(
            action=action,
            channel_dict=channel_dict,
            channel_mentions=channel_mentions,
            entry_type=entry_type,
            failed_snowflakes=failed_snowflakes,
            resolved_channels=resolved_channels,
            snowflake_kwargs=snowflake_kwargs,
        )
        await StateService.send_pages(
            title="Streaming Routes", pages=pages, state=state
        )

    @app_commands.command(name="streams", description="List streaming routes.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @administrator_predicator()
    async def list_streaming_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await StreamService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(
            title="Streaming Routes", pages=pages, state=state
        )

    @app_commands.command(name="vr", description="Start/stop video-only room.")
    @app_commands.describe(channel="Tag a channel or include the ID")
    @administrator_predicator()
    async def toggle_video_room_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        msg = await VideoRoomService.toggle_video_room(channel_dict=channel_dict)
        return await state.end(success=msg)

    @app_commands.command(name="vrs", description="List video rooms.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, or server ID."
    )
    @administrator_predicator()
    async def list_video_rooms_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction=interaction)
        do = DiscordObject(interaction=interaction)
        target = target or int(interaction.channel.id)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await VideoRoomService.build_pages(
            object_dict=object_dict, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Video Rooms", pages=pages, state=state)

    @app_commands.command(name="xalias", description="Delete alias.")
    @administrator_predicator()
    @app_commands.describe(alias_name="Include an alias name")
    async def delete_alias_app_command(
        self, interaction: discord.Interaction, alias_name: str
    ):
        state = StateService(interaction=interaction)
        snowflake_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        msg = await AliasService.delete_alias(
            alias_name=alias_name, snowflake_kwargs=snowflake_kwargs
        )
        return await state.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(AdminAppCommands(bot))
