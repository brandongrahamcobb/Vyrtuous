"""admin_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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

from datetime import datetime, timezone
from typing import Optional

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.actions.alias import Alias
from vyrtuous.database.actions.server_mute import ServerMute
from vyrtuous.database.actions.voice_mute import VoiceMute
from vyrtuous.database.logs.history import History
from vyrtuous.database.roles.administrator import AdministratorRole
from vyrtuous.database.roles.coordinator import Coordinator
from vyrtuous.database.roles.moderator import Moderator
from vyrtuous.database.roles.permission_role import PermissionRole
from vyrtuous.database.rooms.stage import Stage
from vyrtuous.database.rooms.room import (
    member_relevant_objects_dict,
    room_relevant_objects_dict,
)
from vyrtuous.database.rooms.temporary_room import TemporaryRoom
from vyrtuous.database.rooms.video_room import VideoRoom
from vyrtuous.database.settings.cap import Cap
from vyrtuous.properties.duration import AppDuration, Duration, DurationObject
from vyrtuous.properties.moderation_type import AppModerationType, ModerationType
from vyrtuous.properties.snowflake import (
    AppChannelSnowflake,
    AppMemberSnowflake,
    AppRoleSnowflake,
    ChannelSnowflake,
    MemberSnowflake,
    RoleSnowflake,
)
from vyrtuous.service.check_service import (
    administrator_predicator,
    at_home,
    check,
    NotGuildOwner,
)
from vyrtuous.service.logging_service import logger
from vyrtuous.service.messaging.message_service import MessageService
from vyrtuous.service.messaging.state_service import StateService
from vyrtuous.service.resolution.discord_object_service import (
    DiscordObject,
    DiscordObjectNotFound,
)

from vyrtuous.service.scope_service import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    generate_skipped_roles,
    clean_guild_dictionary,
    flush_page,
)
from vyrtuous.utils.cancel_confirm import VerifyView
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.permission import TARGET_PERMISSIONS


class AdminCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot, self.bot.db_pool)

    # DONE
    @app_commands.command(name="alias", description="Alias creation.")
    @administrator_predicator()
    @app_commands.describe(
        alias_name="Alias/Pseudonym",
        moderation_type="One of: vegan, carnist, vmute, unvmute,"
        "ban, unban, flag, unflag, tmute, untmute, role, unrole",
        channel="Tag a channel or include its ID",
        role="Role ID (only for role/unrole)",
    )
    async def create_alias_app_command(
        self,
        interaction: discord.Interaction,
        moderation_type: AppModerationType,
        alias_name: Optional[str],
        channel: AppChannelSnowflake,
        role: AppRoleSnowflake = None,
    ):
        state = StateService(source=interaction)

        do = DiscordObject(source=interaction)

        channel_dict = await do.determine_from_target(target=channel)
        kwargs = channel_dict["columns"]
        try:
            role_dict = await do.determine_from_target(target=role)
        except DiscordObjectNotFound as e:
            logger.warning(str(e).capitalize())
        else:
            kwargs.update(role_dict["columns"])

        alias = await Alias.select(
            alias_name=alias_name, guild_snowflake=interaction.guild.id
        )
        if alias:
            return await state.end(
                warning=f"Alias `{alias.alias_name}` already exists in {interaction.guild.name}."
            )

        alias = Alias(alias_name=alias_name, alias_type=moderation_type, **kwargs)
        await alias.create()

        if role_dict:
            msg = (
                f"Alias `{alias_name}` of type `{moderation_type}` "
                f"created successfully for channel {channel_dict['mention']} "
                f"and role {role_dict['mention']}."
            )
        else:
            msg = (
                f"Alias `{alias_name}`  of type `{moderation_type}` "
                f"created successfully for channel {channel_dict['mention']}."
            )
        return await state.end(success=msg)

    # DONE
    @commands.command(
        name="alias",
        help="Set an alias for a vegan, carnist, mute, unmute, ban, "
        "unban, flag, unflag, tmute, untmute, role, or unrole action.",
    )
    @administrator_predicator()
    async def create_alias_text_command(
        self,
        ctx: commands.Context,
        moderation_type: ModerationType = commands.parameter(
            default=None,
            description="One of: `vegan`, `carnist`, `vmute`, "
            "`unvmute`, `ban`, `unban`, `flag`, "
            "`unflag`, `tmute`, `untmute`, `role`, `unrole`",
        ),
        alias_name: str = commands.parameter(
            default=None, description="Alias/Pseudonym"
        ),
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
        *,
        role: Optional[RoleSnowflake] = commands.parameter(
            default=None, description="Role ID (only for role/unrole)"
        ),
    ):
        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)
        channel_dict = do.determine_from_target(target=channel)
        kwargs = channel_dict["columns"]
        try:
            role_dict = await do.determine_from_target(target=role)
        except DiscordObjectNotFound as e:
            logger.warning(str(e).capitalize())
        else:
            kwargs.update(role_dict["columns"])

        alias = await Alias.select(alias_name=alias_name, guild_snowflake=ctx.guild.id)
        if alias:
            return await state.end(
                warning=f"Alias `{alias.alias_name}` already exists in {ctx.guild.name}."
            )

        alias = Alias(alias_name=alias_name, alias_type=moderation_type, **kwargs)
        await alias.create()

        if role_dict:
            msg = (
                f"Alias `{alias_name}` of type `{moderation_type}` "
                f"created successfully for channel {channel_dict['mention']} "
                f"and role {role_dict['mention']}."
            )
        else:
            msg = (
                f"Alias `{alias_name}`  of type `{moderation_type}` "
                f"created successfully for channel {channel_dict['mention']}."
            )
        return await state.end(success=msg)

    @app_commands.command(name="aroles", description="Administrator roles.")
    @app_commands.describe(target="Specify one of: 'all', server ID or empty.")
    @administrator_predicator()
    async def list_administrator_roles_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {AdministratorRole.PLURAL}"

        state = StateService(source=interaction)

        is_at_home = at_home(source=interaction)

        do = DiscordObject(source=interaction)
        object_dict = do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        administrator_roles = await AdministratorRole.select(**kwargs)

        for administrator_role in administrator_roles:
            guild_dictionary.setdefault(
                administrator_role.guild_snowflake, {"roles": {}}
            )
            guild_dictionary[administrator_role.guild_snowflake]["roles"].setdefault(
                administrator_role.role_snowflake, {}
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_roles = generate_skipped_roles(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_roles=skipped_roles,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for role_snowflake, entry in guild_data.get("roles", {}).items():
                role = guild.get_role(role_snowflake)
                if field_count >= chunk_size:
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                embed.add_field(name=role.name, value=role.mention, inline=False)
                field_count += 1
            pages.append(embed)

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_roles:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_roles,
                    title="Skipped Roles in Server",
                )

        await StateService.send_pages(obj=AdministratorRole, pages=pages, state=state)

    @commands.command(name="aroles", help="Administrator roles.")
    @administrator_predicator()
    async def list_administrator_roles_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: 'all', "
            "channel ID/mention, server ID or empty.",
        ),
    ):
        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {AdministratorRole.PLURAL}"

        state = StateService(source=ctx)

        is_at_home = at_home(source=ctx)

        do = DiscordObject(source=ctx)
        object_dict = do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        administrator_roles = await AdministratorRole.select(**kwargs)

        for administrator_role in administrator_roles:
            guild_dictionary.setdefault(
                administrator_role.guild_snowflake, {"roles": {}}
            )
            guild_dictionary[administrator_role.guild_snowflake]["roles"].setdefault(
                administrator_role.role_snowflake, {}
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_roles = generate_skipped_roles(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_roles=skipped_roles,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for role_snowflake, role_dict in guild_data.get("roles", {}).items():
                role = guild.get_role(role_snowflake)
                embed, field_count = flush_page(embed, pages, title, guild.name)
                embed.add_field(name=role.name, value=role.mention, inline=False)
                field_count += 1
            pages.append(embed)

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_roles:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_roles,
                    title="Skipped Roles in Server",
                )

        await StateService.send_pages(obj=AdministratorRole, pages=pages, state=state)

    @app_commands.command(name="cap", description="Cap alias duration for mods.")
    @administrator_predicator()
    @app_commands.describe(
        channel="Tag a channel or include its ID",
        moderation_type="One of: `mute`, `ban`, `tmute`",
        hours="(+|-)duration(m|h|d), 0=permanent, default=24h",
    )
    async def cap_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        moderation_type: AppModerationType,
        hours: int,
    ):
        state = StateService(source=interaction)

        seconds = int(hours) * 3600

        do = DiscordObject(source=interaction)
        channel_dict = do.determine_from_target(target=channel)

        kwargs = channel_dict["columns"]
        kwargs.update("moderation_type", moderation_type)

        cap = await Cap.select(**kwargs)
        if cap and seconds:
            await Cap.update(
                set_kwargs={"duration_seconds": seconds}, where_kwargs=kwargs
            )
            msg = f"Cap `{moderation_type}` modified for {channel_dict['mention']}."
        elif cap:
            await Cap.delete(**kwargs)
            return await state.end(
                success=f"Cap of type {moderation_type} "
                f"and channel {channel_dict['mention']} deleted successfully."
            )
        else:
            kwargs.update("duration_seconds", seconds)
            cap = Cap(**kwargs)
            await cap.create()
            msg = (
                f"Cap `{moderation_type}` created for "
                f"{channel_dict['mention']} successfully."
            )
        return await state.end(success=msg)

    # DONE
    @commands.command(name="cap", help="Cap alias duration for mods.")
    @administrator_predicator()
    async def cap_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
        moderation_type: ModerationType = commands.parameter(
            default=None, description="One of: `mute`, `ban`, `tmute`"
        ),
        *,
        hours: int = commands.parameter(default=24, description="# of hours"),
    ):
        state = StateService(source=ctx)

        seconds = int(hours) * 3600

        do = DiscordObject(source=ctx)
        channel_dict = do.determine_from_target(target=channel)

        kwargs = channel_dict["columns"]
        kwargs.update("moderation_type", moderation_type)

        cap = await Cap.select(**kwargs)
        if cap and seconds:
            await Cap.update(
                set_kwargs={"duration_seconds": seconds}, where_kwargs=kwargs
            )
            msg = f"Cap `{moderation_type}` modified for {channel_dict['mention']}."
        elif cap:
            await Cap.delete(**kwargs)
            return await state.end(
                success=f"Cap of type {moderation_type} "
                f"and channel {channel_dict['mention']} deleted successfully."
            )
        else:
            kwargs.update("duration_seconds", seconds)
            cap = Cap(**kwargs)
            await cap.create()
            msg = (
                f"Cap `{moderation_type}` created for "
                f"{channel_dict['mention']} successfully."
            )
        return await state.end(success=msg)

    # DONE
    @app_commands.command(
        name="chown", description="Change the owner of a temporary room."
    )
    @app_commands.describe(
        member="Tag a user or provide their ID",
        channel="Tag a channel or provide it's ID",
    )
    @administrator_predicator()
    async def change_temp_room_owner_app_command(
        self, interaction, channel: AppChannelSnowflake, member: AppMemberSnowflake
    ):
        state = StateService(source=interaction)
        do = DiscordObject(source=interaction)

        channel_dict = do.determine_from_target(target=channel)
        member_dict = do.determine_from_target(target=member)

        kwargs = {}
        kwargs.update(channel_dict["columns"])
        kwargs.update(member_dict["columns"])

        await TemporaryRoom.update_owner(**kwargs)

        return await state.end(
            success=f"Temporary room {channel_dict['object'].mention} ownership "
            f"transferred to {member_dict['object'].mention}."
        )

    # DONE
    @commands.command(
        name="chown", help="Change the owner of a temporary room.", hidden=True
    )
    @administrator_predicator()
    async def change_temp_room_owner_text_command(
        self,
        ctx,
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or provide it's ID"
        ),
        member: MemberSnowflake = commands.parameter(
            default=None, description="Tag a user or provide their ID"
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(source=ctx)

        channel_dict = do.determine_from_target(target=channel)
        member_dict = do.determine_from_target(target=member)

        kwargs = {}
        kwargs.update(channel_dict["columns"])
        kwargs.update(member_dict["columns"])

        await TemporaryRoom.update_owner(**kwargs)

        return await state.end(
            success=f"Temporary room {channel_dict['object'].mention} ownership "
            f"transferred to {member_dict['object'].mention}."
        )

    # DONE
    @app_commands.command(name="clear", description="Reset channel/member.")
    @app_commands.describe(
        target="Tag a channel/member or include the ID",
        action_type="Specify one of: `alias`, `all`, `ban`, "
        "`coord`, `flag`, `mod`, `temp`, `tmute`, `track`, `vegan`, `vmute` or `vr`.",
    )
    @administrator_predicator()
    async def clear_channel_access_app_command(
        self, interaction: discord.Interaction, target: str, action_type: str
    ):
        do = DiscordObject(source=interaction)
        object_dict = do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        if object_dict["object_type"] == discord.Member:

            view = VerifyView(
                action_type=action_type, author_snowflake=interaction.user.id, **kwargs
            )
            embed = view.build_embed(action_type=view.action_type, target=view.target)
            await interaction.response.send_message(embed=embed, view=view)
            await view.wait()

            state = StateService(source=interaction)

            if view.result:
                if action_type == "all":
                    for m_obj in member_relevant_objects_dict.values():
                        m_obj.delete(**kwargs)
                    msg = (
                        f"Deleted all associated moderation actions and "
                        f"roles for {object_dict['object'].mention}."
                    )
                for singular, m_obj in member_relevant_objects_dict.items():
                    if action_type.lower() == singular:
                        m_obj.delete(**kwargs)
                        msg = f"Deleted all associated {m_obj.PLURAL} on {object_dict['object'].mention}."

        elif object_dict["object_type"] == discord.VoiceChannel:
            try:
                await check(source=interaction, lowest_role="Guild Owner")
            except NotGuildOwner() as e:
                logger.warning(str(e).capitalize())

            view = VerifyView(
                action_type=action_type, author_snowflake=interaction.user.id, **kwargs
            )
            embed = view.build_embed(action_type=view.action_type, target=view.target)
            await interaction.response.send_message(embed=embed, view=view)
            await view.wait()

            state = StateService(source=interaction)

            if view.result:
                if action_type == "all":
                    for r_obj in room_relevant_objects_dict.values():
                        r_obj.delete(**kwargs)
                    msg = (
                        f"Deleted all associated moderation actions and "
                        f"roles in {object_dict['object'].mention}."
                    )
                for singular, r_obj in room_relevant_objects_dict.items():
                    if action_type.lower() == singular:
                        r_obj.delete(**kwargs)
                        msg = f"Deleted all associated {r_obj.PLURAL} in {object_dict['object'].mention}."

        return await state.end(success=msg)

    # DONE
    @commands.command(name="clear", help="Reset channel/member.")
    @administrator_predicator()
    async def clear_channel_access_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            default=None, description="Tag a channel, a member or include the ID"
        ),
        *,
        action_type: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: `alias`, `all`, `ban`, `coord`, "
            "flag`, `mod`, `temp`, `tmute`, `track`, `vegan`, `vmute` or `vr`.",
        ),
    ):
        do = DiscordObject(source=ctx)
        object_dict = do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        if object_dict["object_type"] == discord.Member:

            view = VerifyView(
                action_type=action_type, author_snowflake=ctx.author.id, **kwargs
            )
            embed = view.build_embed(action_type=view.action_type, target=view.target)
            await ctx.reply(embed=embed, view=view)
            await view.wait()

            state = StateService(source=ctx)

            if view.result:
                if action_type == "all":
                    for m_obj in member_relevant_objects_dict.values():
                        m_obj.delete(**kwargs)
                    msg = (
                        f"Deleted all associated moderation actions and "
                        f"roles for {object_dict['mention']}."
                    )
                for singular, m_obj in member_relevant_objects_dict.items():
                    if action_type.lower() == singular:
                        m_obj.delete(**kwargs)
                        msg = f"Deleted all associated {m_obj.PLURAL} on {object_dict['mention']}."

        elif object_dict["object_type"] == discord.VoiceChannel:
            try:
                await check(source=ctx, lowest_role="Guild Owner")
            except NotGuildOwner() as e:
                logger.warning(str(e).capitalize())

            view = VerifyView(
                action_type=action_type, author_snowflake=ctx.author.id, **kwargs
            )
            embed = view.build_embed(action_type=view.action_type, target=view.target)
            await ctx.reply(embed=embed, view=view)
            await view.wait()

            state = StateService(source=ctx)

            if view.result:
                if action_type == "all":
                    for r_obj in room_relevant_objects_dict.values():
                        r_obj.delete(**kwargs)
                    msg = (
                        f"Deleted all associated moderation actions and "
                        f"roles in {object_dict['mention']}."
                    )
                for singular, r_obj in room_relevant_objects_dict.items():
                    if action_type.lower() == singular:
                        r_obj.delete(**kwargs)
                        msg = f"Deleted all associated {r_obj.PLURAL} in {object_dict['mention']}."

        return await state.end(success=msg)

    # DONE
    @app_commands.command(name="coord", description="Grant/revoke coords.")
    @app_commands.describe(
        member="Tag a member or include their ID",
        channel="Tag a channel or include its ID",
    )
    @administrator_predicator()
    async def create_coordinator_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake,
    ):
        state = StateService(source=interaction)

        do = DiscordObject(source=interaction)
        channel_dict = do.determine_from_target(target=channel)
        member_dict = do.determine_from_target(target=member)
        kwargs = {}
        kwargs.update(channel_dict["columns"])
        kwargs.update(member_dict["columns"])

        coordinator = await Coordinator.select(**kwargs)
        if coordinator:
            await Coordinator.delete(**kwargs)
            action = "revoked"
        else:
            coordinator = Coordinator(**kwargs)
            await coordinator.create()
            action = "granted"

        return await state.end(
            success=f"Coordinator access has been {action} for {member_dict['mention']} "
            f"in {channel_dict['mention']}."
        )

    # DONE
    @commands.command(name="coord", help="Grant/revoke coords.")
    @administrator_predicator()
    async def create_coordinator_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            default=None, description="Tag a member or include their ID"
        ),
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
    ):
        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)
        channel_dict = do.determine_from_target(target=channel)
        member_dict = do.determine_from_target(target=member)
        kwargs = {}
        kwargs.update(channel_dict["columns"])
        kwargs.update(member_dict["columns"])

        coordinator = await Coordinator.select(**kwargs)
        if coordinator:
            await Coordinator.delete(**kwargs)
            action = "revoked"
        else:
            coordinator = Coordinator(**kwargs)
            await coordinator.create()
            action = "granted"

        return await state.end(
            success=f"Coordinator access has been {action} for {member_dict['mention']} "
            f"in {channel_dict['mention']}."
        )

    # @app_commands.command(name='pc', description='View permissions.')
    # @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, server ID or empty.")

    @commands.command(name="pc", help="View permissions.")
    @administrator_predicator()
    async def list_permissions_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            default=None,
            description="Specify one of: `all`, channel "
            "ID/mention, server ID or empty.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {self.bot.user.display_name} Permissions"

        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)
        object_dict = do.determine_from_target(target=target)

        if target and target.lower() == "all":
            await check(source=ctx, lowest_role="Developer")
            channel_objs = [
                channel_obj
                for guild in self.bot.guilds
                for channel_obj in guild.channels
            ]
        elif object_dict["object_type"] == discord.Guild:
            channel_objs = object_dict["object"].channels
        elif object_dict["object_type"] in (
            discord.StageChannel,
            discord.TextChannel,
            discord.VoiceChannel,
        ):
            channel_objs = [object_dict["object"]]

        for channel in channel_objs:
            permissions = channel.permissions_for(self.bot.me)
            missing = []
            for permission in TARGET_PERMISSIONS:
                if not getattr(permissions, permission):
                    missing.append(permission)
            if not missing:
                continue
            guild_dictionary.setdefault(channel.guild.id, {})
            guild_dictionary[channel.guild.id].setdefault("channels", {})
            guild_dictionary[channel.guild.id]["channels"].setdefault(channel.id, {})
            guild_dictionary[channel.guild.id]["channels"][channel.id].update(
                {"Missing Permissions": missing}
            )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                for section_name, permissions in channel_data.items():
                    lines.append(section_name)
                    for permission in permissions:
                        lines.append(f"  ↳ {permission}")
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        await StateService.send_pages(obj=PermissionRole, pages=pages, state=state)

    # DONE
    @app_commands.command(name="rmv", description="VC move.")
    @app_commands.describe(
        source_channel="Tag the source channel or include its ID",
        target_channel="Tag the target channel or include its ID",
    )
    @administrator_predicator()
    async def room_move_all_app_command(
        self,
        interaction: discord.Interaction,
        source_channel: AppChannelSnowflake,
        target_channel: AppChannelSnowflake,
    ):
        failed, moved = [], []
        state = StateService(source=interaction)

        do = DiscordObject(source=interaction)

        source_channel_dict = do.determine_from_target(target=source_channel)
        target_channel_dict = do.determine_from_target(target=target_channel)

        for member in source_channel_dict["object"].members:
            try:
                await member.move_to(target_channel_dict["object"])
                moved.append(member)
            except discord.Forbidden as e:
                failed.append(member)
                logger.warning(
                    f"Unable to move member "
                    f"{member.display_name} ({member.id}) from channel "
                    f"{source_channel_dict['name']} ({source_channel}) to channel "
                    f"{target_channel_dict['name']} ({target_channel}) in guild "
                    f"{interaction.guild.name} ({interaction.guild.id}). "
                    f"{str(e).capitalize()}"
                )
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"Moved {source_channel_dict['mention']} to "
            f"{target_channel_dict['mention']}",
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
            text=f"Moved from {source_channel_dict['mention']} "
            f"to {target_channel_dict['mention']}"
        )

        return await state.end(success=embed)

    # DONE
    @commands.command(name="rmv", help="VC move.")
    @administrator_predicator()
    async def room_move_all_text_command(
        self,
        ctx: commands.Context,
        source_channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
        target_channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
    ):
        failed, moved = [], []
        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)

        source_channel_dict = do.determine_from_target(target=source_channel)
        target_channel_dict = do.determine_from_target(target=target_channel)

        for member in source_channel_dict["object"].members:
            try:
                await member.move_to(target_channel_dict["object"])
                moved.append(member)
            except discord.Forbidden as e:
                failed.append(member)
                logger.warning(
                    f"Unable to move member "
                    f"{member.display_name} ({member.id}) from channel "
                    f"{source_channel_dict['name']} ({source_channel}) to channel "
                    f"{target_channel_dict['name']} ({target_channel}) in guild "
                    f"{ctx.guild.name} ({ctx.guild.id}). "
                    f"{str(e).capitalize()}"
                )
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"Moved {source_channel_dict['mention']} to "
            f"{target_channel_dict['mention']}",
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
            text=f"Moved from {source_channel_dict['mention']} "
            f"to {target_channel_dict['mention']}"
        )

        return await state.end(success=embed)

    # DONE
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
        reason: Optional[str] = "No reason provided",
    ):
        state = StateService(source=interaction)

        do = DiscordObject(source=interaction)

        member_dict = do.determine_from_target(target=member)
        kwargs = member_dict["columns"]

        server_mute = await ServerMute.select(**kwargs)
        if not server_mute:
            server_mute = ServerMute(
                **kwargs,
                reason=reason,
            )
            await server_mute.create()
            action = "muted"
            should_be_muted = True
        else:
            await ServerMute.delete(**kwargs)
            action = "unmuted"
            should_be_muted = False

        if member_dict["object"].voice and member_dict["object"].voice.channel:
            try:
                await member_dict["object"].edit(mute=should_be_muted)
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())

        return await state.end(
            success=f"Successfully server {action} {member_dict['mention']}."
        )

    # DONE
    @commands.command(name="smute", help="Server mute/server unmute.")
    @administrator_predicator()
    async def toggle_server_mute_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            default=None, description="Tag a member or include their ID"
        ),
        *,
        reason: Optional[str] = commands.parameter(
            default="No reason provided",
            description="Optional reason (required for 7 days or more)",
        ),
    ):
        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)

        member_dict = do.determine_from_target(target=member)
        kwargs = member_dict["columns"]

        server_mute = await ServerMute.select(**kwargs)
        if not server_mute:
            server_mute = ServerMute(
                **kwargs,
                reason=reason,
            )
            await server_mute.create()
            action = "muted"
            should_be_muted = True
        else:
            await ServerMute.delete(**kwargs)
            action = "unmuted"
            should_be_muted = False

        if member_dict["object"].voice and member_dict["object"].voice.channel:
            try:
                await member_dict["object"].edit(mute=should_be_muted)
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())

        return await state.end(
            success=f"Successfully server {action} {member_dict['mention']}."
        )

    # DONE
    @app_commands.command(name="stage", description="Start/stop stage.")
    @app_commands.describe(
        channel="Tag a voice/stage channel",
        duration="Duration of the stage (e.g., 1h, 30m)",
    )
    @administrator_predicator()
    async def stage_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        duration: AppDuration = DurationObject("1h"),
    ):
        failed = pages = skipped = succeeded = []
        is_modification = False

        state = StateService(source=interaction)
        do = DiscordObject(source=interaction)

        channel_dict = do.determine_from_target(target=channel)

        if not isinstance(duration, DurationObject):
            duration = DurationObject(duration)
        if duration.is_modification:
            is_modification = True

        kwargs = channel_dict["columns"]

        stage = await Stage.select(**kwargs)
        if is_modification and stage:
            delta = duration.expires_in - datetime.now(timezone.utc)
            if delta.total_seconds() < 0:
                return await state.end(
                    warning="You are not authorized to decrease the duration "
                    "below the current time."
                )
            if stage:
                set_kwargs = {"expires_in": duration.expired_in}
                await Stage.update(set_kwargs=set_kwargs, where_kwargs=kwargs)
                description_lines = [
                    f"**Channel:** {channel_dict['mention']}",
                    f"**Expires:** {duration}",
                ]
                embed = discord.Embed(
                    description="\n".join(description_lines),
                    title=f"{get_random_emoji()} Stage Modified",
                    color=discord.Color.blurple(),
                )
        elif stage:
            title = f"{get_random_emoji()} Stage Ended in {channel_dict['mention']}"
            await Stage.delete(**kwargs)
            for member in channel_dict["object"].members:
                await VoiceMute.delete(
                    **kwargs,
                    member_snowflake=member.id,
                    target="room",
                )
                voice_mute = await VoiceMute.select(
                    **kwargs,
                    member_snowflake=member.id,
                    target="user",
                )
                if not voice_mute and member.voice and member.voice.mute:
                    try:
                        await member.edit(
                            mute=False,
                            reason="Stage ended — no user-specific mute found",
                        )
                        succeeded.append(member)
                    except discord.Forbidden as e:
                        logger.warning(
                            f"Unable to undo voice-mute "
                            f"for member {member.display_name} ({member.id}) in "
                            f"channel {channel_dict['name']} ({channel_dict['id']}) in "
                            f"guild {interaction.guild.name} ({interaction.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict['mention']}",
                f"**Unmuted:** {len(succeeded)} users",
            ]
            if failed:
                description_lines.append(f"**Failed:** {len(failed)}")
            embed = discord.Embed(
                description="\n".join(description_lines),
                title=title,
                color=discord.Color.blurple(),
            )
            pages.append(embed)
        else:
            stage = Stage(
                **kwargs,
                expires_in=duration.expires_in,
            )
            await stage.create()
            for member in channel_dict["object"].members:
                if (
                    await check(
                        source=interaction,
                        member_snowflake=member.id,
                    )
                    or member.id == interaction.user.id
                ):
                    skipped.append(member)
                    continue
                voice_mute = await VoiceMute(
                    **kwargs,
                    expires_in=duration.expires_in,
                    member_snowflake=member.id,
                    target="room",
                    reason="Stage mute",
                )
                await voice_mute.create()
                try:
                    if member.voice and member.voice.channel.id == channel_dict["id"]:
                        await member.edit(mute=True)
                    succeeded.append(member)
                except Exception as e:
                    logger.warning(
                        f"Unable to voice-mute "
                        f"member {member.display_name} ({member.id}) "
                        f"in channel {channel_dict['name']} ({channel_dict['id']}) "
                        f"in guild {interaction.guild.name} ({interaction.guild.id}). "
                        f"{str(e).capitalize()}"
                    )
                    failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict['mention']}",
                f"**Expires:** {duration}",
                f"**Muted:** {len(succeeded)} users",
                f"**Skipped:** {len(skipped)}",
            ]
            if failed:
                description_lines.append(f"**Failed:** {len(failed)}")
            embed = discord.Embed(
                description="\n".join(description_lines),
                title=f"{get_random_emoji()} Stage Created in {channel_dict['name']}",
                color=discord.Color.blurple(),
            )
            pages.append(embed)

        await StateService.send_pages(obj=Stage, pages=pages, state=state)

    # DONE
    @commands.command(name="stage", help="Start/stop stage")
    @administrator_predicator()
    async def stage_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
        *,
        duration: Duration = commands.parameter(
            default=DurationObject("1h"),
            description="Options: (+|-)duration(m|h|d) "
            "0 - permanent / 24h - default",
        ),
    ):
        failed = pages = skipped = succeeded = []
        is_modification = False

        state = StateService(source=ctx)
        do = DiscordObject(source=ctx)

        channel_dict = do.determine_from_target(target=channel)
        kwargs = channel_dict["columns"]

        if not isinstance(duration, DurationObject):
            duration = DurationObject(duration)
        if duration.is_modification:
            is_modification = True

        stage = await Stage.select(**kwargs)
        if is_modification and stage:
            delta = duration.expires_in - datetime.now(timezone.utc)
            if delta.total_seconds() < 0:
                return await state.end(
                    warning="You are not authorized to decrease the duration "
                    "below the current time."
                )
            if stage:
                set_kwargs = {"expires_in": duration.expired_in}
                await Stage.update(set_kwargs=set_kwargs, where_kwargs=kwargs)
                description_lines = [
                    f"**Channel:** {channel_dict['mention']}",
                    f"**Expires:** {duration}",
                ]
                embed = discord.Embed(
                    description="\n".join(description_lines),
                    title=f"{get_random_emoji()} Stage Modified",
                    color=discord.Color.blurple(),
                )
        elif stage:
            title = f"{get_random_emoji()} Stage Ended in {channel_dict['mention']}"
            await Stage.delete(**kwargs)
            for member in channel_dict["object"].members:
                await VoiceMute.delete(
                    **kwargs,
                    member_snowflake=member.id,
                    target="room",
                )
                voice_mute = await VoiceMute.select(
                    **kwargs,
                    member_snowflake=member.id,
                    target="user",
                )
                if not voice_mute and member.voice and member.voice.mute:
                    try:
                        await member.edit(
                            mute=False,
                            reason="Stage ended — no user-specific mute found",
                        )
                        succeeded.append(member)
                    except discord.Forbidden as e:
                        logger.warning(
                            f"Unable to undo voice-mute "
                            f"for member {member.display_name} ({member.id}) in "
                            f"channel {channel_dict['name']} ({channel_dict['id']}) in "
                            f"guild {ctx.guild.name} ({ctx.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict['mention']}",
                f"**Unmuted:** {len(succeeded)} users",
            ]
            if failed:
                description_lines.append(f"**Failed:** {len(failed)}")
            embed = discord.Embed(
                description="\n".join(description_lines),
                title=title,
                color=discord.Color.blurple(),
            )
            pages.append(embed)
        else:
            stage = Stage(
                **kwargs,
                expires_in=duration.expires_in,
            )
            await stage.create()
            for member in channel_dict["object"].members:
                if (
                    await check(
                        source=ctx,
                        lowest_role="Coordinator",
                        member_snowflake=member.id,
                    )
                    or member.id == ctx.author.id
                ):
                    skipped.append(member)
                    continue
                voice_mute = await VoiceMute(
                    **kwargs,
                    expires_in=duration.expires_in,
                    member_snowflake=member.id,
                    target="room",
                    reason="Stage mute",
                )
                await voice_mute.create()
                try:
                    if member.voice and member.voice.channel.id == channel_dict["id"]:
                        await member.edit(mute=True)
                    succeeded.append(member)
                except Exception as e:
                    logger.warning(
                        f"Unable to voice-mute "
                        f"member {member.display_name} ({member.id}) "
                        f"in channel {channel_dict['name']} ({channel_dict['id']}) "
                        f"in guild {ctx.guild.name} ({ctx.guild.id}). "
                        f"{str(e).capitalize()}"
                    )
                    failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict['mention']}",
                f"**Expires:** {duration}",
                f"**Muted:** {len(succeeded)} users",
                f"**Skipped:** {len(skipped)}",
            ]
            if failed:
                description_lines.append(f"**Failed:** {len(failed)}")
            embed = discord.Embed(
                description="\n".join(description_lines),
                title=f"{get_random_emoji()} Stage Created in {channel_dict['name']}",
                color=discord.Color.blurple(),
            )
            pages.append(embed)

        await StateService.send_pages(obj=Stage, pages=pages, state=state)

    # DONE
    @app_commands.command(
        name="temp", description="Toggle a temporary room and assign an owner."
    )
    @app_commands.describe(
        channel="Tag a channel or include its ID",
        owner="Tag a member or include their ID",
    )
    @administrator_predicator()
    async def toggle_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        owner: AppMemberSnowflake,
    ):
        action = None

        state = StateService(source=interaction)
        do = DiscordObject(source=interaction)

        channel_dict = do.determine_from_target(target=channel)
        member_dict = do.determine_from_target(target=owner)

        kwargs = {}
        kwargs.update(channel_dict["columns"])

        temporary_room = await TemporaryRoom.select(**kwargs)
        kwargs.update(member_dict["columns"])
        if temporary_room:
            if temporary_room.member_snowflake:
                await Moderator.delete(**kwargs)
            await TemporaryRoom.delete(**kwargs)
            action = "removed"
        else:
            moderator = Moderator(**kwargs)
            await moderator.create()
            temporary_room = TemporaryRoom(
                **kwargs,
                room_name=channel_dict["name"],
            )
            await temporary_room.create()
            action = "created"

        return await state.end(
            success=f"Temporary room {action} in {channel_dict['mention']}."
        )

    # DONE
    @commands.command(
        name="temp", help="Toggle a temporary room and assign an owner.", hidden=True
    )
    @administrator_predicator()
    async def toggle_temp_room_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
        owner: MemberSnowflake = commands.parameter(
            default=None, description="Tag a member or include their Discord ID"
        ),
    ):
        action = None

        state = StateService(source=ctx)
        do = DiscordObject(source=ctx)

        channel_dict = do.determine_from_target(target=channel)
        member_dict = do.determine_from_target(target=owner)

        kwargs = {}
        kwargs.update(channel_dict["columns"])

        temporary_room = await TemporaryRoom.select(**kwargs)
        kwargs.update(member_dict["columns"])
        if temporary_room:
            if temporary_room.member_snowflake:
                await Moderator.delete(**kwargs)
            await TemporaryRoom.delete(**kwargs)
            action = "removed"
        else:
            moderator = Moderator(**kwargs)
            await moderator.create()
            temporary_room = TemporaryRoom(
                **kwargs,
                room_name=channel_dict["name"],
            )
            await temporary_room.create()
            action = "created"

        return await state.end(
            success=f"Temporary room {action} in {channel_dict['mention']}."
        )

    # DONE
    @app_commands.command(
        name="temps", description="List temporary rooms with matching command aliases."
    )
    @administrator_predicator()
    async def list_temp_rooms_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {TemporaryRoom.PLURAL}"
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)

        do = DiscordObject(source=interaction)

        object_dict = do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        aliases = Alias.select(**kwargs)
        temporary_rooms = TemporaryRoom.select(**kwargs)

        for temporary_room in temporary_rooms:
            guild_dictionary.setdefault(
                temporary_room.guild_snowflake, {"channels": {}}
            )
            guild_dictionary[temporary_room.guild_snowflake]["channels"].setdefault(
                temporary_room.channel_snowflake, {}
            )
            if aliases:
                for alias in aliases:
                    if (
                        alias.guild_snowflake == temporary_room.guild_snowflake
                        and alias.channel_snowflake == temporary_room.channel_snowflake
                    ):
                        guild_dictionary[temporary_room.guild_snowflake]["channels"][
                            temporary_room.channel_snowflake
                        ].setdefault(alias.alias_type, [])
                        guild_dictionary[temporary_room.guild_snowflake]["channels"][
                            temporary_room.channel_snowflake
                        ][alias.alias_type].append(alias.alias_name)

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                lines.append(f"Channel: {channel.mention}")
                field_count += 1
                for alias_type, alias_names in channel_data.items():
                    lines.append(f"{alias_type}")
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= chunk_size:
                            embed.add_field(
                                name="Information",
                                value="\n\n".join(lines),
                                inline=False,
                            )
                            embed, field_count = flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if is_at_home:
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )

        await StateService.send_pages(obj=TemporaryRoom, pages=pages, state=state)

    # DONE
    @commands.command(
        name="temps",
        help="List temporary rooms with matching command aliases.",
        hidden=True,
    )
    @administrator_predicator()
    async def list_temp_rooms_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: `all`, channel ID/mention, "
            "server ID or empty.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {TemporaryRoom.PLURAL}"
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)

        object_dict = do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        aliases = Alias.select(**kwargs)
        temporary_rooms = TemporaryRoom.select(**kwargs)

        for temporary_room in temporary_rooms:
            guild_dictionary.setdefault(
                temporary_room.guild_snowflake, {"channels": {}}
            )
            guild_dictionary[temporary_room.guild_snowflake]["channels"].setdefault(
                temporary_room.channel_snowflake, {}
            )
            if aliases:
                for alias in aliases:
                    if (
                        alias.guild_snowflake == temporary_room.guild_snowflake
                        and alias.channel_snowflake == temporary_room.channel_snowflake
                    ):
                        guild_dictionary[temporary_room.guild_snowflake]["channels"][
                            temporary_room.channel_snowflake
                        ].setdefault(alias.alias_type, [])
                        guild_dictionary[temporary_room.guild_snowflake]["channels"][
                            temporary_room.channel_snowflake
                        ][alias.alias_type].append(alias.alias_name)

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                lines.append(f"Channel: {channel.mention}")
                field_count += 1
                for alias_type, alias_names in channel_data.items():
                    lines.append(f"{alias_type}")
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= chunk_size:
                            embed.add_field(
                                name="Information",
                                value="\n\n".join(lines),
                                inline=False,
                            )
                            embed, field_count = flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if is_at_home:
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )

        await StateService.send_pages(obj=TemporaryRoom, pages=pages, state=state)

    # DONE
    @app_commands.command(name="track", description="Setup tracking.")
    @app_commands.describe(
        channel="Tag a channel or include its ID where the messages will be sent.",
        scope="create | modify | delete.",
        entry_type="all | channel | general.",
        snowflakes="Optional list of channel/member IDs to be tracked.",
    )
    @administrator_predicator()
    async def modify_tracking_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        scope: Optional[str] = None,
        entry_type: Optional[str] = None,
        snowflakes: Optional[str] = None,
    ):
        channel_mentions = failed_snowflakes = []

        state = StateService(source=interaction)

        do = DiscordObject(source=interaction)

        channel_dict = do.determine_from_target(target=channel)
        kwargs = channel_dict["columns"]

        if scope is None and entry_type is None:
            history = await History.select(**kwargs)
            if not history:
                return await state.end(
                    warning=f"No tracking exists for {channel_dict['mention']}."
                )
            enabled = not history[0].enabled
            action = "enabled" if enabled else "disabled"
            await History.update_by_channel_enabled_and_guild(
                channel_snowflake=channel_dict["id"],
                enabled=enabled,
                guild_snowflake=interaction.guild.id,
            )
            return await state.end(
                success=f"Tracking has been {action} in {channel_dict['mention']}."
            )

        if scope and entry_type:
            match scope.lower():
                case "create" | "modify":
                    resolved_channels = []
                    for snowflake in snowflakes:
                        try:
                            channel_dict = do.determine_from_target(target=channel)
                            resolved_channels.append(channel_dict["id"])
                            channel_mentions.append(channel_dict["mention"])
                        except Exception:
                            failed_snowflakes.append(snowflake)
                            continue
                    if scope.lower() == "create":
                        history_entry = History(
                            **kwargs,
                            enabled=True,
                            entry_type=entry_type,
                            snowflakes=resolved_channels,
                        )
                        await history_entry.create()
                        action = "created"
                    else:
                        await History.update_by_channel_guild_and_type(
                            **kwargs,
                            entry_type=entry_type,
                            snowflakes=resolved_channels,
                        )
                        action = "modified"
                case "delete":
                    await History.delete(**kwargs)
                    action = "deleted"
                case _:
                    return await state.end(
                        warning="\U000026a0\U0000fe0f Scope must be one of `create`, `delete` or `modify`."
                    )

        embed = discord.Embed(
            title=f"{get_random_emoji()} Tracking {action.capitalize()} for {channel_dict['mention']}",
            color=0x00FF00,
        )
        if channel_mentions:
            embed.add_field(
                name="Processed Channels",
                value=", ".join(channel_mentions),
                inline=False,
            )
        if failed_snowflakes:
            embed.add_field(
                name="Failed IDs",
                value=", ".join(str(s) for s in failed_snowflakes),
                inline=False,
            )

        return await state.end(success=embed)

    # DONE
    @commands.command(name="track", help="Setup tracking.")
    @administrator_predicator()
    async def modify_tracking_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID."
        ),
        scope: Optional[str] = commands.parameter(
            default=None, description="create | modify | delete."
        ),
        entry_type: Optional[str] = commands.parameter(
            default=None, description="all | channel."
        ),
        *snowflakes: Optional[int],
    ):
        channel_mentions = failed_snowflakes = []

        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)

        channel_dict = do.determine_from_target(target=channel)
        kwargs = channel_dict["columns"]

        if scope is None and entry_type is None:
            history = await History.select(**kwargs)
            if not history:
                return await state.end(
                    warning=f"No tracking exists for {channel_dict['mention']}."
                )
            enabled = not history[0].enabled
            action = "enabled" if enabled else "disabled"
            await History.update_by_channel_enabled_and_guild(
                channel_snowflake=channel_dict["id"],
                enabled=enabled,
                guild_snowflake=ctx.guild.id,
            )
            return await state.end(
                success=f"Tracking has been {action} in {channel_dict['mention']}."
            )

        if scope and entry_type:
            match scope.lower():
                case "create" | "modify":
                    resolved_channels = []
                    for snowflake in snowflakes:
                        try:
                            channel_dict = do.determine_from_target(target=channel)
                            resolved_channels.append(channel_dict["id"])
                            channel_mentions.append(channel_dict["mention"])
                        except Exception:
                            failed_snowflakes.append(snowflake)
                            continue
                    if scope.lower() == "create":
                        history_entry = History(
                            **kwargs,
                            enabled=True,
                            entry_type=entry_type,
                            snowflakes=resolved_channels,
                        )
                        await history_entry.create()
                        action = "created"
                    else:
                        await History.update_by_channel_guild_and_type(
                            **kwargs,
                            entry_type=entry_type,
                            snowflakes=resolved_channels,
                        )
                        action = "modified"
                case "delete":
                    await History.delete(**kwargs)
                    action = "deleted"
                case _:
                    return await state.end(
                        warning="Scope must be one of `create`, `delete` or `modify`."
                    )

        embed = discord.Embed(
            title=f"{get_random_emoji()} Tracking {action.capitalize()} for {channel_dict['mention']}",
            color=0x00FF00,
        )
        if channel_mentions:
            embed.add_field(
                name="Processed Channels",
                value=", ".join(channel_mentions),
                inline=False,
            )
        if failed_snowflakes:
            embed.add_field(
                name="Failed IDs",
                value=", ".join(str(s) for s in failed_snowflakes),
                inline=False,
            )

        return await state.end(success=embed)

    # DONE
    @app_commands.command(name="tracks", description="List history channels.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, server ID or empty."
    )
    @administrator_predicator()
    async def list_tracking_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=interaction)
        title = f"{get_random_emoji} {History.PLURAL}"

        state = StateService(source=interaction)

        do = DiscordObject(source=interaction)

        object_dict = do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        histories = await History.select(**kwargs)

        for history in histories:
            guild_dictionary.setdefault(history.guild_snowflake, {"channels", {}})
            guild_dictionary[history.guild_snowflake]["channels"][
                history.channel_snowflake
            ] = {
                "enabled": history.enabled,
                "entry_type": history.entry_type,
                "snowflakes": history.snowflakes,
            }

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, entries in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                for entry in entries:
                    status = "\u2705" if entry["enabled"] else "\u26d4"
                    channel = guild.get_channel(entry["channel_snowflake"])
                    lines.append(
                        f"{status}**Channel:** {channel.mention}\n**Type:** {entry['entry_type']}"
                    )
                    if object_dict["type"] in (
                        discord.StageChannel,
                        discord.TextChannel,
                        discord.VoiceChannel,
                    ):
                        lines.append(f"**Snowflakes:** {entry['snowflakes']}")
                    field_count += 1
                    if field_count >= chunk_size:
                        embed.add_field(
                            name="Information",
                            value="\n\n".join(lines),
                            inline=False,
                        )
                        embed, field_count = flush_page(embed, pages, title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if is_at_home:
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )

        await StateService.send_pages(obj=History, pages=pages, state=state)

    # DONE
    @commands.command(name="tracks", help="List history channels.")
    @administrator_predicator()
    async def list_tracking_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: `all`, channel ID/mention, server ID or empty.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=ctx)
        title = f"{get_random_emoji} {History.PLURAL}"

        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)

        object_dict = do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        histories = await History.select(**kwargs)

        for history in histories:
            guild_dictionary.setdefault(history.guild_snowflake, {"channels", {}})
            guild_dictionary[history.guild_snowflake]["channels"][
                history.channel_snowflake
            ] = {
                "enabled": history.enabled,
                "entry_type": history.entry_type,
                "snowflakes": history.snowflakes,
            }

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, entries in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                for entry in entries:
                    status = "\u2705" if entry["enabled"] else "\u26d4"
                    channel = guild.get_channel(entry["channel_snowflake"])
                    lines.append(
                        f"{status}**Channel:** {channel.mention}\n**Type:** {entry['entry_type']}"
                    )
                    if object_dict["type"] in (
                        discord.StageChannel,
                        discord.TextChannel,
                        discord.VoiceChannel,
                    ):
                        lines.append(f"**Snowflakes:** {entry['snowflakes']}")
                    field_count += 1
                    if field_count >= chunk_size:
                        embed.add_field(
                            name="Information",
                            value="\n\n".join(lines),
                            inline=False,
                        )
                        embed, field_count = flush_page(embed, pages, title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if is_at_home:
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )

        await StateService.send_pages(obj=History, pages=pages, state=state)

    # DONE
    @app_commands.command(name="vr", description="Start/stop video-only room.")
    @app_commands.describe(channel="Tag a channel or include the ID")
    @administrator_predicator()
    async def toggle_video_room_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        action = None

        state = StateService(source=interaction)

        do = DiscordObject(source=interaction)

        channel_dict = do.determine_from_target(target=channel)
        kwargs = channel_dict["columns"]

        video_room = await VideoRoom.select(**kwargs)

        if video_room:
            action = "removed"
            VideoRoom.video_rooms = [
                vr
                for vr in VideoRoom.video_rooms
                if vr.channel_snowflake != video_room.channel_snowflake
            ]
            await VideoRoom.delete(**kwargs)
        else:
            video_room = VideoRoom(**kwargs)
            await video_room.create()
            VideoRoom.video_rooms.append(video_room)
            action = "created"

        return await state.end(
            success=f"Video-only room {action} in {channel_dict['mention']}."
        )

    @commands.command(name="vr", help="Start/stop video-only room.")
    @administrator_predicator()
    async def toggle_video_room_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include the ID"
        ),
    ):
        action = None

        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)

        channel_dict = do.determine_from_target(target=channel)
        kwargs = channel_dict["columns"]

        video_room = await VideoRoom.select(**kwargs)

        if video_room:
            action = "removed"
            VideoRoom.video_rooms = [
                vr
                for vr in VideoRoom.video_rooms
                if vr.channel_snowflake != video_room.channel_snowflake
            ]
            await VideoRoom.delete(**kwargs)
        else:
            video_room = VideoRoom(**kwargs)
            await video_room.create()
            VideoRoom.video_rooms.append(video_room)
            action = "created"

        return await state.end(
            success=f"Video-only room {action} in {channel_dict['mention']}."
        )

    # DONE
    @app_commands.command(name="vrs", description="List video rooms.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, server ID or empty."
    )
    @administrator_predicator()
    async def list_video_rooms_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {VideoRoom.PLURAL}"
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)

        do = DiscordObject(source=interaction)

        object_dict = do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        aliases = Alias.select(**kwargs)
        video_rooms = VideoRoom.select(**kwargs)

        for video_room in video_rooms:
            guild_dictionary.setdefault(video_room.guild_snowflake, {"channels": {}})
            guild_dictionary[video_room.guild_snowflake]["channels"].setdefault(
                video_room.channel_snowflake, {}
            )
            if aliases:
                for alias in aliases:
                    if (
                        alias.guild_snowflake == video_room.guild_snowflake
                        and alias.channel_snowflake == video_room.channel_snowflake
                    ):
                        guild_dictionary[video_room.guild_snowflake]["channels"][
                            video_room.channel_snowflake
                        ].setdefault(alias.alias_type, [])
                        guild_dictionary[video_room.guild_snowflake]["channels"][
                            video_room.channel_snowflake
                        ][alias.alias_type].append(alias.alias_name)

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                lines.append(f"Channel: {channel.mention}")
                field_count += 1
                for alias_type, alias_names in channel_data.items():
                    lines.append(f"{alias_type}")
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= chunk_size:
                            embed.add_field(
                                name="Information",
                                value="\n\n".join(lines),
                                inline=False,
                            )
                            embed, field_count = flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if is_at_home:
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )

        await StateService.send_pages(obj=VideoRoom, pages=pages, state=state)

    # DONE
    @commands.command(name="vrs", help="List video rooms.")
    @administrator_predicator()
    async def list_video_rooms_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None, description="Include `all`, channel or server ID."
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {VideoRoom.PLURAL}"
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)

        object_dict = do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        aliases = Alias.select(**kwargs)
        video_rooms = VideoRoom.select(**kwargs)

        for video_room in video_rooms:
            guild_dictionary.setdefault(video_room.guild_snowflake, {"channels": {}})
            guild_dictionary[video_room.guild_snowflake]["channels"].setdefault(
                video_room.channel_snowflake, {}
            )
            if aliases:
                for alias in aliases:
                    if (
                        alias.guild_snowflake == video_room.guild_snowflake
                        and alias.channel_snowflake == video_room.channel_snowflake
                    ):
                        guild_dictionary[video_room.guild_snowflake]["channels"][
                            video_room.channel_snowflake
                        ].setdefault(alias.alias_type, [])
                        guild_dictionary[video_room.guild_snowflake]["channels"][
                            video_room.channel_snowflake
                        ][alias.alias_type].append(alias.alias_name)

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                lines.append(f"Channel: {channel.mention}")
                field_count += 1
                for alias_type, alias_names in channel_data.items():
                    lines.append(f"{alias_type}")
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= chunk_size:
                            embed.add_field(
                                name="Information",
                                value="\n\n".join(lines),
                                inline=False,
                            )
                            embed, field_count = flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if is_at_home:
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )

        await StateService.send_pages(obj=VideoRoom, pages=pages, state=state)

    # DONE
    @app_commands.command(name="xalias", description="Delete alias.")
    @administrator_predicator()
    @app_commands.describe(alias_name="Include an alias name")
    async def delete_alias_app_command(
        self, interaction: discord.Interaction, alias_name: str
    ):
        state = StateService(source=interaction)

        do = DiscordObject(source=interaction)

        kwargs = {"alias_name": alias_name, "guild_snowflake": interaction.guild.id}

        alias = await Alias.select(**kwargs)
        if not alias:
            return await state.end(warning=f"No aliases found for `{alias_name}`.")
        await Alias.delete(**kwargs)

        channel_dict = do.determine_from_target(target=alias.channel_snowflake)

        if alias.role_snowflake:
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.alias_type}` for channel {channel_dict['mention']} "
                f" and role {alias.role_mention} deleted successfully."
            )
        else:
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.alias_type}` for channel {channel_dict['mention']} "
                f"deleted successfully."
            )

        return await state.end(success=msg)

    # DONE
    @commands.command(name="xalias", help="Delete alias.")
    @administrator_predicator()
    async def delete_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: str = commands.parameter(
            default=None, description="Include an alias name"
        ),
    ):
        state = StateService(source=ctx)

        do = DiscordObject(source=ctx)

        kwargs = {"alias_name": alias_name, "guild_snowflake": ctx.guild.id}

        alias = await Alias.select(**kwargs)
        if not alias:
            return await state.end(warning=f"No aliases found for `{alias_name}`.")
        await Alias.delete(**kwargs)

        channel_dict = do.determine_from_target(target=alias.channel_snowflake)

        if alias.role_snowflake:
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.alias_type}` for channel {channel_dict['mention']} "
                f" and role {alias.role_mention} deleted successfully."
            )
        else:
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.alias_type}` for channel {channel_dict['mention']} "
                f"deleted successfully."
            )

        return await state.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(AdminCommands(bot))
