"""admin_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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

from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.actions.ban import Ban
from vyrtuous.db.actions.text_mute import TextMute
from vyrtuous.db.actions.server_mute import ServerMute
from vyrtuous.db.actions.voice_mute import VoiceMute
from vyrtuous.db.roles.administrator import (
    AdministratorRole,
    administrator_predicator,
)
from vyrtuous.db.roles.coordinator import Coordinator
from vyrtuous.db.roles.moderator import Moderator
from vyrtuous.db.roles.guild_owner import NotGuildOwner
from vyrtuous.db.roles.sysadmin import is_sysadmin_wrapper
from vyrtuous.db.rooms.stage import Stage
from vyrtuous.db.rooms.temporary_room import TemporaryRoom
from vyrtuous.db.rooms.video_room import VideoRoom
from vyrtuous.db.mgmt.cap import Cap
from vyrtuous.db.mgmt.stream import Streaming
from vyrtuous.inc.helpers import PATH_LOG
from vyrtuous.fields.duration import AppDuration, Duration, DurationObject
from vyrtuous.fields.category import AppCategory, Category
from vyrtuous.fields.snowflake import (
    AppChannelSnowflake,
    AppRoleSnowflake,
    AppMemberSnowflake,
    ChannelSnowflake,
    MemberSnowflake,
    RoleSnowflake,
)
from vyrtuous.utils.home import at_home
from vyrtuous.utils.check import (
    check,
    has_equal_or_lower_role_wrapper,
    HasEqualOrLowerRole,
)
from vyrtuous.utils.dir_to_classes import dir_to_classes
from vyrtuous.utils.logger import logger
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import (
    DiscordObject,
    DiscordObjectNotFound,
)
from vyrtuous.utils.guild_dictionary import (
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
from vyrtuous.inc.helpers import TARGET_PERMISSIONS


class AdminCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    # DONE
    @app_commands.command(name="alias", description="Alias creation.")
    @administrator_predicator()
    @app_commands.describe(
        alias_name="Alias/Pseudonym",
        category="Specify an category for a `ban`, `flag`, `hide`, `role`, `tmute`, `vegan` or `vmute` action.",
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
        role: AppRoleSnowflake = None,
    ):
        msg = ""
        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)
        kwargs = {}
        try:
            channel_dict = await do.determine_from_target(target=str(channel))
        except (DiscordObjectNotFound, TypeError) as e:
            logger.warning(str(e).capitalize())
        else:
            kwargs.update(channel_dict.get("columns", None))
            msg = (
                f"Alias `{alias_name}`  of type `{category}` "
                f"created successfully for channel {channel_dict.get("mention", None)}."
            )
        if category in ("ban", "hide", "tmute", "role"):
            if not role and str(category) in ("hide", "tmute"):
                for role in interaction.guild.roles:
                    if role.name == alias_name:
                        return await state.end(
                            warning=f"{role.name} already exists. You must specify it to override."
                        )
                role_obj = await interaction.guild.create_role(name=alias_name)
                role = str(role_obj.id)
            elif not role and category == "ban":
                return await state.end(
                    warning="You must provide a preexisting text-mute role to associate with the ban."
                )
            if role:
                try:
                    role_dict = await do.determine_from_target(target=role)
                except (DiscordObjectNotFound, TypeError) as e:
                    logger.warning(str(e).capitalize())
                    return await state.end(
                        warning="Provided role ({role}) was not found."
                    )
                else:
                    if category == "ban":
                        text_mute = await TextMute.select(role_snowflake=role_dict["id"])
                        if not text_mute:
                            return await state.end(
                                warning=f"{role_dict.get("mention", None)} is not a text-mute role. You must provide a preexisting text-mute role to associate with the ban."
                            )
                    kwargs.update(role_dict.get("columns", None))
                    msg = (
                        f"Alias `{alias_name}` of type `{category}` "
                        f"created successfully "
                        f"for role {role_dict.get("mention", None)}."
                    )
            else:
                return await state.end(
                    warning=f"The alias type ({category}) requires a role."
                )

        alias = await Alias.select(
            alias_name=alias_name, guild_snowflake=interaction.guild.id, singular=True
        )

        if alias:
            return await state.end(
                warning=f"Alias `{alias.alias_name}` already exists in {interaction.guild.name}."
            )

        alias = Alias(alias_name=alias_name, category=str(category), **kwargs)
        await alias.create()
        return await state.end(success=msg)

    # DONE
    @commands.command(
        name="alias",
        help="Alias creation.",
    )
    @administrator_predicator()
    async def create_alias_text_command(
        self,
        ctx: commands.Context,
        category: Category = commands.parameter(
            description="Specify an category for a `ban`, `flag`, `hide`, `role`, `tmute`, `vegan` or `vmute` action."
        ),
        alias_name: str = commands.parameter(description="Alias/Pseudonym"),
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include the ID"
        ),
        *,
        role: RoleSnowflake = commands.parameter(
            default=None, description="Tag a role or include the ID."
        ),
    ):
        msg = ""
        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)
        kwargs = {}
        try:
            channel_dict = await do.determine_from_target(target=channel)
        except (DiscordObjectNotFound, TypeError) as e:
            logger.warning(str(e).capitalize())
        else:
            kwargs.update(channel_dict.get("columns", None))
            msg = (
                f"Alias `{alias_name}`  of type `{category}` "
                f"created successfully for channel {channel_dict.get("mention", None)}."
            )
        if category in ("ban", "hide", "tmute", "role"):
            if not role and category in ("hide", "tmute"):
                for role in ctx.guild.roles:
                    if role.name == alias_name:
                        return await state.end(
                            warning=f"{role.name} already exists. You must specify it to override."
                        )
                role_obj = await ctx.guild.create_role(name=alias_name)
                role = str(role_obj.id)
            elif not role and category == "ban":
                return await state.end(
                    warning="You must provide a preexisting text-mute role to associate with the ban."
                )
            if role:
                try:
                    role_dict = await do.determine_from_target(target=role)
                except (DiscordObjectNotFound, TypeError) as e:
                    logger.warning(str(e).capitalize())
                    return await state.end(
                        warning="Provided role ({role}) was not found."
                    )
                else:
                    if category == "ban":
                        text_mute = await TextMute.select(role_snowflake=role_dict["id"])
                        if not text_mute:
                            return await state.end(
                                warning=f"{role_dict.get("mention", None)} is not a text-mute role. You must provide a preexisting text-mute role to associate with the ban."
                            )
                    kwargs.update(role_dict.get("columns", None))
                    msg = (
                        f"Alias `{alias_name}` of type `{category}` "
                        f"created successfully "
                        f"for role {role_dict.get("mention", None)}."
                    )
            else:
                return await state.end(
                    warning=f"The alias type ({category}) requires a role."
                )

        alias = await Alias.select(
            alias_name=alias_name, guild_snowflake=ctx.guild.id, singular=True
        )

        if alias:
            return await state.end(
                warning=f"Alias `{alias.alias_name}` already exists in {ctx.guild.name}."
            )

        alias = Alias(alias_name=alias_name, category=str(category), **kwargs)
        await alias.create()
        return await state.end(success=msg)

    @app_commands.command(name="aroles", description="Administrator roles.")
    @app_commands.describe(target="Specify one of: 'all', or server ID.")
    @administrator_predicator()
    async def list_administrator_roles_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {AdministratorRole.PLURAL}"

        state = StateService(source=interaction)

        is_at_home = at_home(source=interaction)

        do = DiscordObject(interaction=interaction)
        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict.get("columns", None)

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

        await StateService.send_pages(
            plural=AdministratorRole.PLURAL, pages=pages, state=state
        )

    @commands.command(name="aroles", help="Administrator roles.")
    @administrator_predicator()
    async def list_administrator_roles_text_command(
        self,
        ctx: commands.Context,
        *,
        target: str = commands.parameter(
            description="Specify one of: 'all', " "channel ID/mention, or server ID.",
        ),
    ):
        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {AdministratorRole.PLURAL}"

        state = StateService(source=ctx)

        is_at_home = at_home(source=ctx)

        do = DiscordObject(ctx=ctx)
        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict.get("columns", None)

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

        await StateService.send_pages(
            plural=AdministratorRole.PLURAL, pages=pages, state=state
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
        state = StateService(source=interaction)

        seconds = int(hours) * 3600

        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)

        kwargs = channel_dict.get("columns", None)
        kwargs.update({"category": category})

        cap = await Cap.select(**kwargs)
        if cap and seconds:
            await Cap.update(
                set_kwargs={"duration_seconds": seconds}, where_kwargs=kwargs
            )
            msg = f"Cap `{category}` modified for {channel_dict.get("mention", None)}."
        elif cap:
            await Cap.delete(**kwargs)
            return await state.end(
                success=f"Cap of type {category} "
                f"and channel {channel_dict.get("mention", None)} deleted successfully."
            )
        else:
            kwargs.update({"duration_seconds": seconds})
            cap = Cap(**kwargs)
            await cap.create()
            msg = (
                f"Cap `{category}` created for "
                f"{channel_dict.get("mention", None)} successfully."
            )

        return await state.end(success=msg)

    # DONE
    @commands.command(name="cap", help="Cap alias duration for mods.")
    @administrator_predicator()
    async def cap_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
        category: Category = commands.parameter(
            description="One of: `mute`, `ban`, `tmute`"
        ),
        *,
        hours: int = commands.parameter(default=24, description="# of hours"),
    ):
        state = StateService(source=ctx)

        seconds = int(hours) * 3600

        do = DiscordObject(ctx=ctx)
        channel_dict = await do.determine_from_target(target=channel)

        kwargs = channel_dict.get("columns", None)
        kwargs.update({"category": category})

        cap = await Cap.select(**kwargs)
        if cap and seconds:
            await Cap.update(
                set_kwargs={"duration_seconds": seconds}, where_kwargs=kwargs
            )
            msg = f"Cap `{category}` modified for {channel_dict.get("mention", None)}."
        elif cap:
            await Cap.delete(**kwargs)
            return await state.end(
                success=f"Cap of type {category} "
                f"and channel {channel_dict.get("mention", None)} deleted successfully."
            )
        else:
            kwargs.update({"duration_seconds": seconds})
            cap = Cap(**kwargs)
            await cap.create()
            msg = (
                f"Cap `{category}` created for "
                f"{channel_dict.get("mention", None)} successfully."
            )
        return await state.end(success=msg)

    # DONE
    @app_commands.command(name="caps", description="List caps.")
    @administrator_predicator()
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    async def list_caps_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        object_dict = await do.determine_from_target(target=target)
        pages = await Cap.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Cap.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="caps", help="List caps.")
    @administrator_predicator()
    async def list_caps_text_command(
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
        pages = await Cap.build_pages(object_dict=object_dict, is_at_home=is_at_home)
        await StateService.send_pages(plural=Cap.PLURAL, pages=pages, state=state)

    # DONE
    @app_commands.command(
        name="chown", description="Change the owner of a temporary room."
    )
    @app_commands.describe(
        member="Tag a user or provide their ID",
        channel="Tag a channel or provide its ID.",
    )
    @administrator_predicator()
    async def change_temp_room_owner_app_command(
        self, interaction, channel: AppChannelSnowflake, member: AppMemberSnowflake
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)

        where_kwargs = channel_dict.get("columns", None)
        set_kwargs = member_dict.get("columns", None)
        await TemporaryRoom.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)

        return await state.end(
            success=f"Temporary room {channel_dict.get('object', None).mention} ownership "
            f"transferred to {member_dict.get('object', None).mention}."
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
            description="Tag a channel or provide its ID."
        ),
        member: MemberSnowflake = commands.parameter(
            description="Tag a user or provide their ID"
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)
        member_dict = await do.determine_from_target(target=member)

        where_kwargs = channel_dict.get("columns", None)
        set_kwargs = member_dict.get("columns", None)
        await TemporaryRoom.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)

        return await state.end(
            success=f"Temporary room {channel_dict.get('object', None).mention} ownership "
            f"transferred to {member_dict.get('object', None).mention}."
        )

    # DONE
    @app_commands.command(name="clear", description="Reset database.")
    @app_commands.describe(
        target="Specify 'all', tag a channel/guild/member or include its ID",
        category="Specify one of: `admin`, `alias`, `arole`, `all`, `ban`, `coord`, "
        "flag`, `mod`, `temp`, `tmute`, `stage`, `stream`, `vegan`, `vmute` or `vr`.",
    )
    @administrator_predicator()
    async def clear_channel_access_app_command(
        self,
        interaction: discord.Interaction,
        target: str,
        category: AppCategory = "all",
    ):
        do = DiscordObject(interaction=interaction)
        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict.get("columns", None)

        dir_paths = []
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/actions")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/mgmt")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/roles")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/rooms")
        view = VerifyView(
            category=str(category),
            mention=object_dict.get("mention", "All"),
            author_snowflake=interaction.user.id,
            **kwargs,
        )
        embed = view.build_embed()
        await interaction.response.send_message(embed=embed, view=view)
        await view.wait()
        state = StateService(source=interaction)
        if isinstance(object_dict.get("object", None), discord.Member):
            try:
                await has_equal_or_lower_role_wrapper(
                    source=interaction,
                    member_snowflake=object_dict.get("id"),
                    sender_snowflake=interaction.user.id,
                )
            except HasEqualOrLowerRole() as e:
                logger.warning(str(e).capitalize())
                return state.end(warning=str(e).capitalize())

            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if "member" in obj.SCOPES:
                        if str(category) == "all":
                            await obj.delete(**kwargs)
                            msg = f"Deleted all associated database information for {object_dict.get('mention', None)}."
                        elif str(category).lower() == obj.CATEGORY:
                            await obj.delete(**kwargs)
                            msg = f"Deleted all associated {obj.PLURAL.lower()} for {object_dict.get('mention', None)}."
                        if isinstance(obj, Ban):
                            bans = Ban.select(
                                guild_snowflake=interaction.guild.id,
                                member_snowflake=object_dict.get("id", None),
                            )
                            for ban in bans:
                                await obj.revoke_role(
                                    guild_snowflake=interaction.guild.id,
                                    member_snowflake=ban.member_snowflake,
                                    role_snowflake=ban.role_snowflake,
                                )
                        elif isinstance(obj, TextMute):
                            text_mutes = TextMute.select(
                                guild_snowflake=interaction.guild.id,
                                member_snowflake=object_dict.get("id", None),
                            )
                            for text_mute in text_mutes:
                                await obj.revoke_role(
                                    guild_snowflake=interaction.guild.id,
                                    member_snowflake=text_mute.member_snowflake,
                                    role_snowflake=text_mute.role_snowflake,
                                )
        elif isinstance(object_dict.get("object", None), discord.abc.GuildChannel):
            try:
                await check(source=interaction, lowest_role="Guild Owner")
            except NotGuildOwner() as e:
                logger.warning(str(e).capitalize())
                return state.end(warning=str(e).capitalize())

            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if "channel" in obj.SCOPES:
                        if category == "all":
                            await obj.delete(**kwargs)
                            msg = f"Deleted all database information for {object_dict.get('mention')}."
                        elif str(category).lower() == obj.CATEGORY:
                            await obj.delete(**kwargs)
                            msg = f"Deleted all associated {obj.PLURAL.lower()} in {object_dict.get('mention', None)}."
                        if isinstance(obj, Ban):
                            bans = Ban.select(
                                channel_snowflake=object_dict.get("id", None),
                                guild_snowflake=interaction.guild.id,
                            )
                            for ban in bans:
                                await obj.revoke_role(
                                    guild_snowflake=interaction.guild.id,
                                    member_snowflake=ban.member_snowflake,
                                    role_snowflake=ban.role_snowflake,
                                )
                        elif isinstance(obj, TextMute):
                            text_mutes = TextMute.select(
                                channel_snowflake=object_dict.get("id", None),
                                guild_snowflake=interaction.guild.id,
                            )
                            for text_mute in text_mutes:
                                await obj.revoke_role(
                                    guild_snowflake=interaction.guild.id,
                                    member_snowflake=text_mute.member_snowflake,
                                    role_snowflake=text_mute.role_snowflake,
                                )
        elif isinstance(object_dict.get("object", None), discord.Guild):
            try:
                await check(source=interaction, lowest_role="Guild Owner")
            except NotGuildOwner() as e:
                logger.warning(str(e).capitalize())
                return state.end(warning=str(e).capitalize())

            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if any(
                        scope in obj.SCOPES for scope in ("guild", "channel", "member")
                    ):
                        if str(category).lower() == "all":
                            await obj.delete(**kwargs)
                            msg = f"Deleted all database information for {object_dict.get('name')}."
                        elif str(category).lower() == obj.CATEGORY:
                            await obj.delete(**kwargs)
                            msg = f"Deleted all associated {obj.PLURAL.lower()} in {object_dict.get('name', None)}."
                        if isinstance(obj, Ban):
                            bans = Ban.select(
                                guild_snowflake=interaction.guild.id,
                            )
                            for ban in bans:
                                await obj.revoke_role(
                                    guild_snowflake=interaction.guild.id,
                                    member_snowflake=ban.member_snowflake,
                                    role_snowflake=ban.role_snowflake,
                                )
                        elif isinstance(obj, TextMute):
                            text_mutes = TextMute.select(
                                guild_snowflake=interaction.guild.id,
                            )
                            for text_mute in text_mutes:
                                await obj.revoke_role(
                                    guild_snowflake=interaction.guild.id,
                                    member_snowflake=text_mute.member_snowflake,
                                    role_snowflake=text_mute.role_snowflake,
                                )
                        elif isinstance(obj, AdministratorRole):
                            administrator_roles = AdministratorRole.select(
                                guild_snowflake=interaction.guild.id,
                            )
                            for administrator_role in administrator_roles:
                                await obj.revoke_role(
                                    guild_snowflake=interaction.guild.id,
                                    role_snowflake=administrator_role.role_snowflake,
                                )
        elif target == "all" and await is_sysadmin_wrapper(source=interaction):
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    await obj.delete(**kwargs)
                    msg = "Deleted all database entries."
                    if isinstance(obj, Ban):
                        bans = Ban.select(
                            guild_snowflake=interaction.guild.id,
                        )
                        for ban in bans:
                            await obj.revoke_role(
                                guild_snowflake=interaction.guild.id,
                                member_snowflake=ban.member_snowflake,
                                role_snowflake=ban.role_snowflake,
                            )
                    elif isinstance(obj, TextMute):
                        text_mutes = TextMute.select(
                            guild_snowflake=interaction.guild.id,
                        )
                        for text_mute in text_mutes:
                            await obj.revoke_role(
                                guild_snowflake=interaction.guild.id,
                                member_snowflake=text_mute.member_snowflake,
                                role_snowflake=text_mute.role_snowflake,
                            )
                    elif isinstance(obj, AdministratorRole):
                        administrator_roles = AdministratorRole.select(
                            guild_snowflake=interaction.guild.id,
                        )
                        for administrator_role in administrator_roles:
                            await obj.revoke_role(
                                guild_snowflake=interaction.guild.id,
                                role_snowflake=administrator_role.role_snowflake,
                            )
        else:
            state = StateService(source=interaction)
            msg = f"Invalid target ({target})."

        return await state.end(success=msg)

    # DONE
    @commands.command(name="clear", help="Reset database.")
    @administrator_predicator()
    async def clear_channel_access_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify 'all', tag a channel/guild/member or include its ID"
        ),
        *,
        category: Category = commands.parameter(
            default="all",
            description="Specify one of: `alias`, `arole`, `all`, `ban`, `coord`, "
            "flag`, `mod`, `temp`, `tmute`, `stage`, `stream`, `vegan`, `vmute` or `vr`.",
        ),
    ):
        do = DiscordObject(ctx=ctx)
        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict.get("columns", None)

        dir_paths = []
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/actions")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/mgmt")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/roles")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/rooms")
        view = VerifyView(
            category=str(category),
            mention=object_dict.get("mention", "All"),
            author_snowflake=ctx.author.id,
            **kwargs,
        )
        embed = view.build_embed()
        await ctx.reply(embed=embed, view=view)
        await view.wait()
        state = StateService(source=ctx)
        if isinstance(object_dict.get("object", None), discord.Member):
            try:
                await has_equal_or_lower_role_wrapper(
                    source=ctx,
                    member_snowflake=object_dict.get("id"),
                    sender_snowflake=ctx.author.id,
                )
            except HasEqualOrLowerRole() as e:
                logger.warning(str(e).capitalize())

            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if "member" in obj.SCOPES:
                        if str(category) == "all":
                            await obj.delete(**kwargs)
                            msg = f"Deleted all associated database information for {object_dict.get('mention', None)}."
                        elif str(category).lower() == obj.CATEGORY:
                            await obj.delete(**kwargs)
                            msg = f"Deleted all associated {obj.PLURAL.lower()} for {object_dict.get('mention', None)}."
                        if isinstance(obj, Ban):
                            bans = Ban.select(
                                guild_snowflake=obj.guild.id,
                                member_snowflake=object_dict.get("id", None),
                            )
                            for ban in bans:
                                await obj.revoke_role(
                                    guild_snowflake=obj.guild.id,
                                    member_snowflake=ban.member_snowflake,
                                    role_snowflake=ban.role_snowflake,
                                )
                        elif isinstance(obj, TextMute):
                            text_mutes = TextMute.select(
                                guild_snowflake=obj.guild.id,
                                member_snowflake=object_dict.get("id", None),
                            )
                            for text_mute in text_mutes:
                                await obj.revoke_role(
                                    guild_snowflake=obj.guild.id,
                                    member_snowflake=text_mute.member_snowflake,
                                    role_snowflake=text_mute.role_snowflake,
                                )
        elif isinstance(object_dict.get("object", None), discord.abc.GuildChannel):
            try:
                await check(source=ctx, lowest_role="Guild Owner")
            except NotGuildOwner() as e:
                logger.warning(str(e).capitalize())
                return state.end(warning=str(e).capitalize())

            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if "channel" in obj.SCOPES:
                        if category == "all":
                            await obj.delete(**kwargs)
                            msg = f"Deleted all database information for {object_dict.get('mention')}."
                        elif str(category).lower() == obj.CATEGORY:
                            await obj.delete(**kwargs)
                            msg = f"Deleted all associated {obj.PLURAL.lower()} in {object_dict.get('mention', None)}."
                        if isinstance(obj, Ban):
                            bans = Ban.select(
                                channel_snowflake=object_dict.get("id", None),
                                guild_snowflake=obj.guild.id,
                            )
                            for ban in bans:
                                await obj.revoke_role(
                                    guild_snowflake=obj.guild.id,
                                    member_snowflake=ban.member_snowflake,
                                    role_snowflake=ban.role_snowflake,
                                )
                        elif isinstance(obj, TextMute):
                            text_mutes = TextMute.select(
                                channel_snowflake=object_dict.get("id", None),
                                guild_snowflake=obj.guild.id,
                            )
                            for text_mute in text_mutes:
                                await obj.revoke_role(
                                    guild_snowflake=obj.guild.id,
                                    member_snowflake=text_mute.member_snowflake,
                                    role_snowflake=text_mute.role_snowflake,
                                )
        elif isinstance(object_dict.get("object", None), discord.Guild):
            try:
                await check(source=ctx, lowest_role="Guild Owner")
            except NotGuildOwner() as e:
                logger.warning(str(e).capitalize())
                return state.end(warning=str(e).capitalize())

            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if any(
                        scope in obj.SCOPES for scope in ("guild", "channel", "member")
                    ):
                        if str(category).lower() == "all":
                            await obj.delete(**kwargs)
                            msg = f"Deleted all database information for {object_dict.get('name')}."
                        elif str(category).lower() == obj.CATEGORY:
                            await obj.delete(**kwargs)
                            msg = f"Deleted all associated {obj.PLURAL.lower()} in {object_dict.get('name', None)}."
                        if isinstance(obj, Ban):
                            bans = Ban.select(
                                guild_snowflake=obj.guild.id,
                            )
                            for ban in bans:
                                await obj.revoke_role(
                                    guild_snowflake=obj.guild.id,
                                    member_snowflake=ban.member_snowflake,
                                    role_snowflake=ban.role_snowflake,
                                )
                        elif isinstance(obj, TextMute):
                            text_mutes = TextMute.select(
                                guild_snowflake=obj.guild.id,
                            )
                            for text_mute in text_mutes:
                                await obj.revoke_role(
                                    guild_snowflake=obj.guild.id,
                                    member_snowflake=text_mute.member_snowflake,
                                    role_snowflake=text_mute.role_snowflake,
                                )
                        elif isinstance(obj, AdministratorRole):
                            administrator_roles = AdministratorRole.select(
                                guild_snowflake=obj.guild.id,
                            )
                            for administrator_role in administrator_roles:
                                await obj.revoke_role(
                                    guild_snowflake=obj.guild.id,
                                    role_snowflake=administrator_role.role_snowflake,
                                )
        elif target == "all" and await is_sysadmin_wrapper(source=ctx):
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    await obj.delete(**kwargs)
                msg = "Deleted all database entries."
                if isinstance(obj, Ban):
                    bans = Ban.select(
                        guild_snowflake=obj.guild.id,
                    )
                    for ban in bans:
                        await obj.revoke_role(
                            guild_snowflake=obj.guild.id,
                            member_snowflake=ban.member_snowflake,
                            role_snowflake=ban.role_snowflake,
                        )
                elif isinstance(obj, TextMute):
                    text_mutes = TextMute.select(
                        guild_snowflake=obj.guild.id,
                    )
                    for text_mute in text_mutes:
                        await obj.revoke_role(
                            guild_snowflake=obj.guild.id,
                            member_snowflake=text_mute.member_snowflake,
                            role_snowflake=text_mute.role_snowflake,
                        )
                elif isinstance(obj, AdministratorRole):
                    administrator_roles = AdministratorRole.select(
                        guild_snowflake=obj.guild.id,
                    )
                    for administrator_role in administrator_roles:
                        await obj.revoke_role(
                            guild_snowflake=obj.guild.id,
                            role_snowflake=administrator_role.role_snowflake,
                        )
        else:
            state = StateService(source=ctx)
            msg = f"Invalid target ({target})."

        return await state.end(success=msg)

    # DONE
    @app_commands.command(name="coord", description="Grant/revoke coords.")
    @app_commands.describe(
        member="Tag a member or include their ID",
        channel="Tag a channel or include its ID.",
    )
    @administrator_predicator()
    async def create_coordinator_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake,
    ):
        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)
        channel_dict = await do.determine_from_target(target=channel)
        if not isinstance(channel_dict.get("object", None), discord.abc.GuildChannel):
            return await state.end(warning=f"Invalid channel ID ({channel}).")
        member_dict = await do.determine_from_target(target=member)
        if not isinstance(member_dict.get("object", None), discord.Member):
            return await state.end(warning=f"Invalid member ID ({channel}).")
        await has_equal_or_lower_role_wrapper(
            source=interaction,
            member_snowflake=member_dict.get("id", None),
            sender_snowflake=interaction.user.id,
        )
        kwargs = {}
        kwargs.update(channel_dict.get("columns", None))
        kwargs.update(member_dict.get("columns", None))

        coordinator = await Coordinator.select(**kwargs)
        if coordinator:
            await Coordinator.delete(**kwargs)
            action = "revoked"
        else:
            coordinator = Coordinator(**kwargs)
            await coordinator.create()
            action = "granted"

        return await state.end(
            success=f"Coordinator access has been {action} for {member_dict.get("mention", None)} "
            f"in {channel_dict.get("mention", None)}."
        )

    # DONE
    @commands.command(name="coord", help="Grant/revoke coords.")
    @administrator_predicator()
    async def create_coordinator_text_command(
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
        if not isinstance(channel_dict.get("object", None), discord.abc.GuildChannel):
            return await state.end(warning=f"Invalid channel ID ({channel}).")
        member_dict = await do.determine_from_target(target=member)
        if not isinstance(member_dict.get("object", None), discord.Member):
            return await state.end(warning=f"Invalid member ID ({channel}).")

        await has_equal_or_lower_role_wrapper(
            source=ctx,
            member_snowflake=member_dict.get("id", None),
            sender_snowflake=ctx.author.id,
        )
        kwargs = {}
        kwargs.update(channel_dict.get("columns", None))
        kwargs.update(member_dict.get("columns", None))
        coordinator = await Coordinator.select(**kwargs)
        if coordinator:
            await Coordinator.delete(**kwargs)
            action = "revoked"
        else:
            coordinator = Coordinator(**kwargs)
            await coordinator.create()
            action = "granted"

        return await state.end(
            success=f"Coordinator access has been {action} for {member_dict.get("mention", None)} "
            f"in {channel_dict.get("mention", None)}."
        )

    @app_commands.command(
        name="debug", description="Shows the last `n` number of logging."
    )
    @administrator_predicator()
    async def debug_app_command(self, interaction: discord.Interaction, lines: int = 3):
        state = StateService(source=interaction)
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

    @commands.command(name="debug", help="Shows the last `n` number of logging.")
    @administrator_predicator()
    async def debug_text_command(
        self,
        ctx,
        *,
        lines: Optional[int] = commands.parameter(
            default=3, description="Specify the number of lines"
        ),
    ):
        state = StateService(source=ctx)
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
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {self.bot.user.display_name} Missing Permissions"

        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)
        object_dict = await do.determine_from_target(target=target)

        if target and target.lower() == "all":
            await check(source=interaction, lowest_role="Developer")
            channel_objs = [
                channel_obj
                for guild in self.bot.guilds
                for channel_obj in guild.channels
            ]
        elif hasattr(object_dict.get("object", None), "channels"):
            channel_objs = object_dict.get("object", None).channels
        else:
            channel_objs = [object_dict.get("object", None)]
        for channel in channel_objs:
            permissions = channel.permissions_for(interaction.guild.me)
            missing = []
            for permission in TARGET_PERMISSIONS:
                if not getattr(permissions, permission):
                    missing.append(permission)
            if not missing:
                continue
            guild_dictionary.setdefault(channel.guild.id, {"channels": {}})
            guild_dictionary[channel.guild.id]["channels"].setdefault(channel.id, {})
            guild_dictionary[channel.guild.id]["channels"][channel.id].update(
                {"permissions": missing}
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
                for section_name, permissions in channel_data.items():
                    for permission in permissions:
                        lines.append(f"  ↳ {permission}")
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if pages:
            return await state.end(warning=pages)
        else:
            return await state.end(
                success=f"{self.bot.user.display_name} has all permissions for `{target}`."
            )

    @commands.command(name="pc", help="View permissions.")
    @administrator_predicator()
    async def list_permissions_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: `all`, channel " "ID/mention or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {self.bot.user.display_name} Missing Permissions"

        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)
        object_dict = await do.determine_from_target(target=target)

        if target and target.lower() == "all":
            await check(source=ctx, lowest_role="Developer")
            channel_objs = [
                channel_obj
                for guild in self.bot.guilds
                for channel_obj in guild.channels
            ]
        elif hasattr(object_dict.get("object", None), "channels"):
            channel_objs = object_dict.get("object", None).channels
        else:
            channel_objs = [object_dict.get("object", None)]
        for channel in channel_objs:
            permissions = channel.permissions_for(ctx.guild.me)
            missing = []
            for permission in TARGET_PERMISSIONS:
                if not getattr(permissions, permission):
                    missing.append(permission)
            if not missing:
                continue
            guild_dictionary.setdefault(channel.guild.id, {"channels": {}})
            guild_dictionary[channel.guild.id]["channels"].setdefault(channel.id, {})
            guild_dictionary[channel.guild.id]["channels"][channel.id].update(
                {"permissions": missing}
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
                for section_name, permissions in channel_data.items():
                    for permission in permissions:
                        lines.append(f"  ↳ {permission}")
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if pages:
            return await state.end(warning=pages)
        else:
            return await state.end(
                success=f"{self.bot.user.display_name} has all permissions for `{target}`."
            )

    # DONE
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
        state = StateService(source=interaction)

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

    # DONE
    @commands.command(name="rmv", help="VC move.")
    @administrator_predicator()
    async def room_move_all_text_command(
        self,
        ctx: commands.Context,
        source_channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
        target_channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
    ):
        failed, moved = [], []
        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

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
                    f"{ctx.guild.name} ({ctx.guild.id}). "
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
        reason: str = "No reason provided",
    ):
        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)

        member_dict = await do.determine_from_target(target=member)
        await has_equal_or_lower_role_wrapper(
            source=interaction,
            member_snowflake=member_dict.get("id", None),
            sender_snowflake=interaction.user.id,
        )
        kwargs = member_dict.get("columns", None)

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

        if (
            member_dict.get("object", None).voice
            and member_dict.get("object", None).voice.channel
        ):
            try:
                await member_dict.get("object", None).edit(mute=should_be_muted)
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())

        return await state.end(
            success=f"Successfully server {action} {member_dict.get("mention", None)}."
        )

    # DONE
    @commands.command(name="smute", help="Server mute/server unmute.")
    @administrator_predicator()
    async def toggle_server_mute_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            description="Tag a member or include their ID"
        ),
        *,
        reason: Optional[str] = commands.parameter(
            default="No reason provided",
            description="Optional reason (required for 7 days or more)",
        ),
    ):
        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        member_dict = await do.determine_from_target(target=member)
        await has_equal_or_lower_role_wrapper(
            source=ctx,
            member_snowflake=member_dict.get("id", None),
            sender_snowflake=ctx.author.id,
        )
        kwargs = member_dict.get("columns", None)

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

        if (
            member_dict.get("object", None).voice
            and member_dict.get("object", None).voice.channel
        ):
            try:
                await member_dict.get("object", None).edit(mute=should_be_muted)
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())

        return await state.end(
            success=f"Successfully server {action} {member_dict.get("mention", None)}."
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
        failed, pages, skipped, succeeded = [], [], [], []
        is_modification = False

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)

        if not isinstance(duration, DurationObject):
            duration = DurationObject(duration)
        if duration.is_modification:
            is_modification = True

        kwargs = channel_dict.get("columns", None)

        stage = await Stage.select(**kwargs)
        if is_modification and stage:
            delta = duration.expires_in - datetime.now(timezone.utc)
            if delta.total_seconds() < 0:
                return await state.end(
                    warning="Member is not authorized to decrease the duration "
                    "below the current time."
                )
            if stage:
                set_kwargs = {"expires_in": duration.expired_in}
                await Stage.update(set_kwargs=set_kwargs, where_kwargs=kwargs)
                description_lines = [
                    f"**Channel:** {channel_dict.get("mention", None)}",
                    f"**Expires:** {duration}",
                ]
                embed = discord.Embed(
                    description="\n".join(description_lines),
                    title=f"{get_random_emoji()} Stage Modified",
                    color=discord.Color.blurple(),
                )
        elif stage:
            title = f"{get_random_emoji()} Stage Ended in {channel_dict.get("mention", None)}"
            await Stage.delete(**kwargs)
            for member in channel_dict.get("object", None).members:
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
                            f"channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) in "
                            f"guild {interaction.guild.name} ({interaction.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict.get("mention", None)}",
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
            for member in channel_dict.get("object", None).members:
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
                    if member.voice and member.voice.channel.id == channel_dict.get(
                        "id", None
                    ):
                        await member.edit(mute=True)
                    succeeded.append(member)
                except Exception as e:
                    logger.warning(
                        f"Unable to voice-mute "
                        f"member {member.display_name} ({member.id}) "
                        f"in channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) "
                        f"in guild {interaction.guild.name} ({interaction.guild.id}). "
                        f"{str(e).capitalize()}"
                    )
                    failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict.get("mention", None)}",
                f"**Expires:** {duration}",
                f"**Muted:** {len(succeeded)} users",
                f"**Skipped:** {len(skipped)}",
            ]
            if failed:
                description_lines.append(f"**Failed:** {len(failed)}")
            embed = discord.Embed(
                description="\n".join(description_lines),
                title=f"{get_random_emoji()} Stage Created in {channel_dict.get('name', None)}",
                color=discord.Color.blurple(),
            )
            pages.append(embed)

        await StateService.send_pages(plural=Stage.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="stage", help="Start/stop stage")
    @administrator_predicator()
    async def stage_text_command(
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
        failed, pages, skipped, succeeded = [], [], [], []
        is_modification = False

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)
        kwargs = channel_dict.get("columns", None)

        if not isinstance(duration, DurationObject):
            duration = DurationObject(duration)
        if duration.is_modification:
            is_modification = True
        stage = await Stage.select(**kwargs, singular=True)
        if is_modification and stage:
            delta = duration.expires_in - datetime.now(timezone.utc)
            if delta.total_seconds() < 0:
                return await state.end(
                    warning="Member is not authorized to decrease the duration "
                    "below the current time."
                )
            if stage:
                set_kwargs = {"expires_in": duration.expired_in}
                await Stage.update(set_kwargs=set_kwargs, where_kwargs=kwargs)
                description_lines = [
                    f"**Channel:** {channel_dict.get("mention", None)}",
                    f"**Expires:** {duration}",
                ]
                embed = discord.Embed(
                    description="\n".join(description_lines),
                    title=f"{get_random_emoji()} Stage Modified",
                    color=discord.Color.blurple(),
                )
        elif stage:
            title = f"{get_random_emoji()} Stage Ended in {channel_dict.get("mention", None)}"
            await Stage.delete(**kwargs)
            for member in channel_dict.get("object", None).members:
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
                            f"channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) in "
                            f"guild {ctx.guild.name} ({ctx.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict.get("mention", None)}",
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
            for member in channel_dict.get("object", None).members:
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
                    if member.voice and member.voice.channel.id == channel_dict.get(
                        "id", None
                    ):
                        await member.edit(mute=True)
                    succeeded.append(member)
                except Exception as e:
                    logger.warning(
                        f"Unable to voice-mute "
                        f"member {member.display_name} ({member.id}) "
                        f"in channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) "
                        f"in guild {ctx.guild.name} ({ctx.guild.id}). "
                        f"{str(e).capitalize()}"
                    )
                    failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict.get("mention", None)}",
                f"**Expires:** {duration}",
                f"**Muted:** {len(succeeded)} users",
                f"**Skipped:** {len(skipped)}",
            ]
            if failed:
                description_lines.append(f"**Failed:** {len(failed)}")
            embed = discord.Embed(
                description="\n".join(description_lines),
                title=f"{get_random_emoji()} Stage Created in {channel_dict.get('name', None)}",
                color=discord.Color.blurple(),
            )
            pages.append(embed)

        await StateService.send_pages(plural=Stage.PLURAL, pages=pages, state=state)

    # DONE
    @app_commands.command(
        name="temp", description="Toggle a temporary room and assign an owner."
    )
    @app_commands.describe(
        channel="Tag a channel or include its ID.",
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
        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)

        kwargs = {}
        kwargs.update(channel_dict.get("columns", None))

        temporary_room = await TemporaryRoom.select(**kwargs, singular=True)
        if temporary_room:
            if temporary_room.member_snowflake:
                kwargs.update({"member_snowflake": temporary_room.member_snowflake})
                await Moderator.delete(**kwargs)
            await TemporaryRoom.delete(**kwargs)
            action = "removed"
        else:
            member_dict = await do.determine_from_target(target=owner)
            kwargs.update(member_dict.get("columns", None))
            moderator = Moderator(**kwargs)
            await moderator.create()
            temporary_room = TemporaryRoom(
                **kwargs,
                room_name=channel_dict.get("name", None),
            )
            await temporary_room.create()
            action = "created"

        return await state.end(
            success=f"Temporary room {action} in {channel_dict.get("mention", None)}."
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
            description="Tag a channel or include its ID."
        ),
        *,
        owner: MemberSnowflake = commands.parameter(
            default=None, description="Tag a member or include their Discord ID"
        ),
    ):
        action = None

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)

        kwargs = {}
        kwargs.update(channel_dict.get("columns", None))

        temporary_room = await TemporaryRoom.select(**kwargs, singular=True)
        if temporary_room:
            if temporary_room.member_snowflake:
                kwargs.update({"member_snowflake": temporary_room.member_snowflake})
                await Moderator.delete(**kwargs)
            await TemporaryRoom.delete(**kwargs)
            action = "removed"
        else:
            member_dict = await do.determine_from_target(target=owner)
            kwargs.update(member_dict.get("columns", None))
            moderator = Moderator(**kwargs)
            await moderator.create()
            temporary_room = TemporaryRoom(
                **kwargs,
                room_name=channel_dict.get("name", None),
            )
            await temporary_room.create()
            action = "created"

        return await state.end(
            success=f"Temporary room {action} in {channel_dict.get("mention", None)}."
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

        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict.get("columns", None)

        aliases = await Alias.select(**kwargs)
        temporary_rooms = await TemporaryRoom.select(**kwargs)

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
                        ].setdefault(alias.category, [])
                        guild_dictionary[temporary_room.guild_snowflake]["channels"][
                            temporary_room.channel_snowflake
                        ][alias.category].append(alias.alias_name)

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
                for category, alias_names in channel_data.items():
                    lines.append(f"{category}")
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= chunk_size:
                            embed.add_field(
                                name="Information",
                                value="\n".join(lines),
                                inline=False,
                            )
                            embed, field_count = flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
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

        await StateService.send_pages(
            plural=TemporaryRoom.PLURAL, pages=pages, state=state
        )

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
        target: str = commands.parameter(
            description="Specify one of: `all`, channel ID/mention, " "or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {TemporaryRoom.PLURAL}"
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict.get("columns", None)

        aliases = await Alias.select(**kwargs)
        temporary_rooms = await TemporaryRoom.select(**kwargs)

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
                        ].setdefault(alias.category, [])
                        guild_dictionary[temporary_room.guild_snowflake]["channels"][
                            temporary_room.channel_snowflake
                        ][alias.category].append(alias.alias_name)

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
                for category, alias_names in channel_data.items():
                    lines.append(f"{category}")
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= chunk_size:
                            embed.add_field(
                                name="Information",
                                value="\n".join(lines),
                                inline=False,
                            )
                            embed, field_count = flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
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

        await StateService.send_pages(
            plural=TemporaryRoom.PLURAL, pages=pages, state=state
        )

    # DONE
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
        action: str = None,
        entry_type: str = None,
        snowflakes: str = None,
    ):
        channel_mentions, failed_snowflakes = [], []

        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)
        kwargs = channel_dict.get("columns", None)

        if action is None and entry_type is None:
            stream = await Streaming.select(**kwargs)
            if not stream:
                return await state.end(
                    warning=f"No streaming exists for {channel_dict.get("mention", None)}."
                )
            enabled = not stream[0].enabled
            action = "enabled" if enabled else "disabled"
            where_kwargs = {
                "channel_snowflake": kwargs["channel_snowflake"],
                "entry_type": entry_type,
                "guild_snowflake": interaction.guild.id,
            }
            set_kwargs = {"enabled": enabled}
            await Streaming.update(
                channel_snowflake=channel_dict.get("id", None),
                enabled=enabled,
                guild_snowflake=interaction.guild.id,
            )
            return await state.end(
                success=f"Streaming has been {action} in {channel_dict.get("mention", None)}."
            )

        if action and entry_type:
            match action.lower():
                case "create" | "modify":
                    resolved_channels = []
                    if snowflakes:
                        for snowflake in snowflakes:
                            try:
                                channel_dict = await do.determine_from_target(
                                    target=channel
                                )
                                resolved_channels.append(channel_dict.get("id", None))
                                channel_mentions.append(
                                    channel_dict.get("mention", None)
                                )
                            except Exception:
                                failed_snowflakes.append(snowflake)
                                continue
                    if action.lower() == "create":
                        stream = Streaming(
                            **kwargs,
                            enabled=True,
                            entry_type=entry_type,
                            snowflakes=resolved_channels,
                        )
                        await stream.create()
                        action = "created"
                    else:
                        where_kwargs = {
                            "channel_snowflake": kwargs["channel_snowflake"],
                            "entry_type": entry_type,
                            "guild_snowflake": interaction.guild.id,
                        }
                        set_kwargs = {"snowflakes": resolved_channels}
                        await Streaming.update(
                            set_kwargs=set_kwargs, where_kwargs=where_kwargs
                        )
                        action = "modified"
                case "delete":
                    await Streaming.delete(**kwargs)
                    action = "deleted"
                case _:
                    return await state.end(
                        warning="Scope must be one of `create`, `delete` or `modify`."
                    )
        embed = discord.Embed(
            title=f"{get_random_emoji()} Tracking {action.capitalize()} for {channel_dict.get("mention", None)}",
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
    @commands.command(name="stream", help="Setup streaming.")
    @administrator_predicator()
    async def modify_streaming_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include its ID."
        ),
        action: str = commands.parameter(description="create | modify | delete."),
        entry_type: str = commands.parameter(description="all | channel."),
        *snowflakes: Optional[int],
    ):
        channel_mentions, failed_snowflakes = [], []

        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)
        kwargs = channel_dict.get("columns", None)

        if action is None and entry_type is None:
            stream = await Streaming.select(**kwargs)
            if not stream:
                return await state.end(
                    warning=f"No streaming exists for {channel_dict.get("mention", None)}."
                )
            enabled = not stream[0].enabled
            action = "enabled" if enabled else "disabled"
            where_kwargs = {
                "channel_snowflake": kwargs["channel_snowflake"],
                "entry_type": entry_type,
                "guild_snowflake": ctx.guild.id,
            }
            set_kwargs = {"enabled": enabled}
            await Streaming.update(
                channel_snowflake=channel_dict.get("id", None),
                enabled=enabled,
                guild_snowflake=ctx.guild.id,
            )
            return await state.end(
                success=f"Streaming has been {action} in {channel_dict.get("mention", None)}."
            )

        if action and entry_type:
            match action.lower():
                case "create" | "modify":
                    resolved_channels = []
                    if snowflakes:
                        for snowflake in snowflakes:
                            try:
                                channel_dict = await do.determine_from_target(
                                    target=channel
                                )
                                resolved_channels.append(channel_dict.get("id", None))
                                channel_mentions.append(
                                    channel_dict.get("mention", None)
                                )
                            except Exception:
                                failed_snowflakes.append(snowflake)
                                continue
                    if action.lower() == "create":
                        stream = Streaming(
                            **kwargs,
                            enabled=True,
                            entry_type=entry_type,
                            snowflakes=resolved_channels,
                        )

                        await stream.create()
                        action = "created"
                    else:
                        where_kwargs = {
                            "channel_snowflake": kwargs["channel_snowflake"],
                            "entry_type": entry_type,
                            "guild_snowflake": ctx.guild.id,
                        }
                        set_kwargs = {"snowflakes": resolved_channels}
                        await Streaming.update(
                            set_kwargs=set_kwargs, where_kwargs=where_kwargs
                        )
                        action = "modified"
                case "delete":
                    await Streaming.delete(**kwargs)
                    action = "deleted"
                case _:
                    return await state.end(
                        warning="Scope must be one of `create`, `delete` or `modify`."
                    )
        embed = discord.Embed(
            title=f"{get_random_emoji()} Tracking {action.capitalize()} for {channel_dict.get("mention", None)}",
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
    @app_commands.command(name="streams", description="List streaming routes.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @administrator_predicator()
    async def list_streaming_app_command(
        self, interaction: discord.Interaction, target: str = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=interaction)
        title = f"{get_random_emoji()} {Streaming.PLURAL}"

        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict.get("columns", None)

        streaming = await Streaming.select(**kwargs)

        for stream in streaming:
            guild_dictionary.setdefault(stream.guild_snowflake, {"channels": {}})
            guild_dictionary[stream.guild_snowflake]["channels"][
                stream.channel_snowflake
            ] = {
                "enabled": stream.enabled,
                "entry_type": stream.entry_type,
                "snowflakes": stream.snowflakes,
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
            for channel_snowflake, entry in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                status = "\u2705" if entry["enabled"] else "\u26d4"
                lines.append(
                    f"{status}**Channel:** {channel.mention}\n**Type:** {entry['entry_type']}"
                )
                if isinstance(
                    object_dict.get("object", None), discord.abc.GuildChannel
                ):
                    lines.append(f"**Snowflakes:** {entry['snowflakes']}")
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
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

        await StateService.send_pages(plural=Streaming.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="streams", help="List streaming routes.")
    @administrator_predicator()
    async def list_streaming_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            description="Specify one of: `all`, channel ID/mention, or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=ctx)
        title = f"{get_random_emoji()} {Streaming.PLURAL}"

        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict.get("columns", None)

        streaming = await Streaming.select(**kwargs)

        for stream in streaming:
            guild_dictionary.setdefault(stream.guild_snowflake, {"channels": {}})
            guild_dictionary[stream.guild_snowflake]["channels"][
                stream.channel_snowflake
            ] = {
                "enabled": stream.enabled,
                "entry_type": stream.entry_type,
                "snowflakes": stream.snowflakes,
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
            for channel_snowflake, entry in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                status = "\u2705" if entry["enabled"] else "\u26d4"
                lines.append(
                    f"{status}**Channel:** {channel.mention}\n**Type:** {entry['entry_type']}"
                )
                if isinstance(
                    object_dict.get("object", None), discord.abc.GuildChannel
                ):
                    lines.append(f"**Snowflakes:** {entry['snowflakes']}")
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
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
        await StateService.send_pages(plural=Streaming.PLURAL, pages=pages, state=state)

    # DONE
    @app_commands.command(name="vr", description="Start/stop video-only room.")
    @app_commands.describe(channel="Tag a channel or include the ID")
    @administrator_predicator()
    async def toggle_video_room_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        action = None

        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)
        kwargs = channel_dict.get("columns", None)

        video_room = await VideoRoom.select(**kwargs, singular=True)

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
            success=f"Video-only room {action} in {channel_dict.get("mention", None)}."
        )

    @commands.command(name="vr", help="Start/stop video-only room.")
    @administrator_predicator()
    async def toggle_video_room_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            description="Tag a channel or include the ID"
        ),
    ):
        action = None

        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)
        kwargs = channel_dict.get("columns", None)

        video_room = await VideoRoom.select(**kwargs, singular=True)

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
            success=f"Video-only room {action} in {channel_dict.get("mention", None)}."
        )

    # DONE
    @app_commands.command(name="vrs", description="List video rooms.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, or server ID."
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

        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict.get("columns", None)

        aliases = await Alias.select(**kwargs)
        video_rooms = await VideoRoom.select(**kwargs)

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
                        ].setdefault(alias.category, [])
                        guild_dictionary[video_room.guild_snowflake]["channels"][
                            video_room.channel_snowflake
                        ][alias.category].append(alias.alias_name)

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
                for category, alias_names in channel_data.items():
                    lines.append(f"{category}")
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= chunk_size:
                            embed.add_field(
                                name="Information",
                                value="\n".join(lines),
                                inline=False,
                            )
                            embed, field_count = flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
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

        await StateService.send_pages(plural=VideoRoom.PLURAL, pages=pages, state=state)

    # DONE
    @commands.command(name="vrs", help="List video rooms.")
    @administrator_predicator()
    async def list_video_rooms_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            description="Include `all`, channel or server ID."
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {VideoRoom.PLURAL}"
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict.get("columns", None)

        aliases = await Alias.select(**kwargs)
        video_rooms = await VideoRoom.select(**kwargs)

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
                        ].setdefault(alias.category, [])
                        guild_dictionary[video_room.guild_snowflake]["channels"][
                            video_room.channel_snowflake
                        ][alias.category].append(alias.alias_name)

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
                for category, alias_names in channel_data.items():
                    lines.append(f"{category}")
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= chunk_size:
                            embed.add_field(
                                name="Information",
                                value="\n".join(lines),
                                inline=False,
                            )
                            embed, field_count = flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
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

        await StateService.send_pages(plural=VideoRoom.PLURAL, pages=pages, state=state)

    # DONE
    @app_commands.command(name="xalias", description="Delete alias.")
    @administrator_predicator()
    @app_commands.describe(alias_name="Include an alias name")
    async def delete_alias_app_command(
        self, interaction: discord.Interaction, alias_name: str
    ):
        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)

        kwargs = {"alias_name": alias_name, "guild_snowflake": interaction.guild.id}

        alias = await Alias.select(singular=True, **kwargs)
        if not alias:
            return await state.end(warning=f"No aliases found for `{alias_name}`.")
        await Alias.delete(**kwargs)

        channel_dict = await do.determine_from_target(
            target=str(alias.channel_snowflake)
        )
        if getattr(alias, "role_snowflake"):
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel_dict.get("mention", None)} "
                f" and role {alias.role_mention} deleted successfully."
            )
        else:
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel_dict.get("mention", None)} "
                f"deleted successfully."
            )

        return await state.end(success=msg)

    # DONE
    @commands.command(name="xalias", help="Delete alias.")
    @administrator_predicator()
    async def delete_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: str = commands.parameter(description="Include an alias name"),
    ):
        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        kwargs = {"alias_name": alias_name, "guild_snowflake": ctx.guild.id}

        alias = await Alias.select(singular=True, **kwargs)
        if not alias:
            return await state.end(warning=f"No aliases found for `{alias_name}`.")
        await Alias.delete(**kwargs)

        channel_dict = await do.determine_from_target(
            target=str(alias.channel_snowflake)
        )
        if getattr(alias, "role_snowflake"):
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel_dict.get("mention", None)} "
                f" and role {alias.role_mention} deleted successfully."
            )
        else:
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel_dict.get("mention", None)} "
                f"deleted successfully."
            )

        return await state.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(AdminCommands(bot))
