"""moderator_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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

from typing import Optional

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.actions.alias import Alias
from vyrtuous.database.actions.ban import Ban
from vyrtuous.database.settings.cap import Cap
from vyrtuous.database.actions.flag import Flag
from vyrtuous.database.actions.text_mute import TextMute
from vyrtuous.database.actions.voice_mute import VoiceMute
from vyrtuous.database.roles.administrator import Administrator
from vyrtuous.database.roles.coordinator import Coordinator
from vyrtuous.database.roles.developer import Developer
from vyrtuous.database.roles.moderator import Moderator
from vyrtuous.database.roles.permission_role import PermissionRole
from vyrtuous.database.roles.vegan import Vegan
from vyrtuous.database.rooms.stage import Stage
from vyrtuous.database.rooms.temporary_room import TemporaryRoom
from vyrtuous.properties.duration import DurationObject
from vyrtuous.properties.snowflake import AppMessageSnowflake, MessageSnowflake
from vyrtuous.properties.snowflake import (
    AppChannelSnowflake,
    AppMemberSnowflake,
    ChannelSnowflake,
    MemberSnowflake,
)
from vyrtuous.service.check_service import (
    at_home,
    check,
    has_equal_or_lower_role,
    member_is_administrator,
    member_is_coordinator,
    member_is_developer,
    member_is_guild_owner,
    member_is_moderator,
    member_is_system_owner,
    moderator_predicator,
)
from vyrtuous.service.logging_service import logger
from vyrtuous.service.messaging.message_service import MessageService
from vyrtuous.service.scope_service import member_relevant_objects_dict
from vyrtuous.service.messaging.state_service import StateService
from vyrtuous.service.resolution.discord_object_service import DiscordObject
from vyrtuous.service.scope_service import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    generate_skipped_members,
    clean_guild_dictionary,
    flush_page,
)
from vyrtuous.utils.emojis import get_random_emoji


class ModeratorCommands(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.message_service = MessageService(self.bot, self.bot.db_pool)

    @app_commands.command(name="admins", description="Lists admins.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, " "or server ID."
    )
    @moderator_predicator()
    async def list_administrators_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Administrator.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        administrators = await Administrator.select(**kwargs)

        for administrator in administrators:
            guild_dictionary.setdefault(administrator.guild_snowflake, {"members": {}})
            guild_dictionary[administrator.guild_snowflake]["members"].setdefault(
                administrator.member_snowflake, {"administrators": {}}
            )

            for role_snowflake in administrator.role_snowflakes:
                guild_dictionary[administrator.guild_snowflake]["members"][
                    administrator.member_snowflake
                ]["administrators"].update({role_snowflake: True})

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, administrators_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                primary_dictionary = administrators_dictionary.get("administrators", {})
                role_mentions = [
                    guild.get_role(role_snowflake).mention
                    for role_snowflake in primary_dictionary
                    if guild.get_role(role_snowflake)
                ]
                if not thumbnail and object_dict.get("type", None) == discord.Member:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                else:
                    member_line = f"**User:** {member.mention}"
                    if role_mentions:
                        member_line += "**Roles:** " + "\n".join(role_mentions)
                lines.append(member_line)
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n\n".join(lines), inline=False
                )
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
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Administrator, pages=pages, state=state)

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
        print("test")
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Administrator.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        administrators = await Administrator.select(**kwargs)
        for administrator in administrators:
            guild_dictionary.setdefault(administrator.guild_snowflake, {"members": {}})
            guild_dictionary[administrator.guild_snowflake]["members"].setdefault(
                administrator.member_snowflake, {"administrators": {}}
            )

            for role_snowflake in administrator.role_snowflakes:
                guild_dictionary[administrator.guild_snowflake]["members"][
                    administrator.member_snowflake
                ]["administrators"].update({role_snowflake: True})

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, administrators_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                primary_dictionary = administrators_dictionary.get("administrators", {})
                role_mentions = [
                    guild.get_role(role_snowflake).mention
                    for role_snowflake in primary_dictionary
                    if guild.get_role(role_snowflake)
                ]
                if not thumbnail and object_dict.get("type", None) == discord.Member:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                else:
                    member_line = f"**User:** {member.mention}"
                    if role_mentions:
                        member_line += "**Roles:** " + "\n".join(role_mentions)
                        lines.append(member_line)
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n\n".join(lines), inline=False
                )
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
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )
        await StateService.send_pages(obj=Administrator, pages=pages, state=state)

    @app_commands.command(name="bans", description="List bans.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_bans_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Ban.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        bans = await Ban.select(**kwargs)

        for ban in bans:
            guild_dictionary.setdefault(ban.guild_snowflake, {"members": {}})
            guild_dictionary[ban.guild_snowflake]["members"].setdefault(
                ban.member_snowflake, {"bans": {}}
            )
            guild_dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                "bans"
            ].setdefault(ban.channel_snowflake, {})
            guild_dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                "bans"
            ][ban.channel_snowflake] = {
                "reason": ban.reason,
                "expires_in": DurationObject.from_expires_in(ban.expires_in),
            }

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, ban_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in ban_dictionary.get(
                    "bans"
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
                    lines.append(f"**Expires in:** {channel_dictionary['expires_in']}")
                    if object_dict.get("type", None) == discord.Member:
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Ban, pages=pages, state=state)

    # DONE
    @commands.command(name="bans", description="List bans.")
    @moderator_predicator()
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            description="Specify one of: 'all', channel ID/mention or server ID.",
        ),
    ):
        print("test")
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        title = f"{get_random_emoji()} {Ban.PLURAL} ({target})"
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        bans = await Ban.select(**kwargs)

        for ban in bans:
            guild_dictionary.setdefault(ban.guild_snowflake, {"members": {}})
            guild_dictionary[ban.guild_snowflake]["members"].setdefault(
                ban.member_snowflake, {"bans": {}}
            )
            guild_dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                "bans"
            ].setdefault(ban.channel_snowflake, {})
            guild_dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                "bans"
            ][ban.channel_snowflake] = {
                "reason": ban.reason,
                "expires_in": DurationObject.from_expires_in(ban.expires_in),
            }

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, ban_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in ban_dictionary.get(
                    "bans"
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
                    lines.append(f"**Expires in:** {channel_dictionary['expires_in']}")
                    if object_dict.get("type", None) == discord.Member:
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Ban, pages=pages, state=state)

    # DONE
    @app_commands.command(name="caps", description="List caps.")
    @moderator_predicator()
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    async def list_caps_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=interaction)
        title = f"{get_random_emoji()} {Cap.PLURAL}"

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        caps = await Cap.select(**kwargs)
        for cap in caps:
            guild_dictionary.setdefault(cap.guild_snowflake, {"channels": {}})
            guild_dictionary[cap.guild_snowflake]["channels"].setdefault(
                cap.channel_snowflake, {"caps": {}}
            )
            guild_dictionary[cap.guild_snowflake]["channels"][cap.channel_snowflake][
                "caps"
            ][cap.moderation_type] = cap.duration_seconds

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
            for channel_snowflake, cap_dictionary in guild_data.get("channels").items():
                channel = guild.get_channel(channel_snowflake)
                for moderation_type, duration_seconds in cap_dictionary.get(
                    "caps", {}
                ).items():
                    lines.append(
                        f"  ↳ {moderation_type} ({DurationObject.from_seconds(duration_seconds)})"
                    )
                    field_count += 1
                    if field_count >= chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed, field_count = flush_page(embed, pages, title, guild.name)
                        lines = []
                if lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
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

        await StateService.send_pages(obj=Cap, pages=pages, state=state)

    # DONE
    @commands.command(name="caps", help="List caps.")
    @moderator_predicator()
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            description="Specify one of: 'all', channel ID/mention or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=ctx)
        title = f"{get_random_emoji()} {Cap.PLURAL}"

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        caps = await Cap.select(**kwargs)

        for cap in caps:
            guild_dictionary.setdefault(cap.guild_snowflake, {"channels": {}})
            guild_dictionary[cap.guild_snowflake]["channels"].setdefault(
                cap.channel_snowflake, {"caps": {}}
            )
            guild_dictionary[cap.guild_snowflake]["channels"][cap.channel_snowflake][
                "caps"
            ][cap.moderation_type] = cap.duration_seconds

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
            for channel_snowflake, cap_dictionary in guild_data.get("channels").items():
                channel = guild.get_channel(channel_snowflake)
                for moderation_type, duration_seconds in cap_dictionary.get(
                    "caps", {}
                ).items():
                    lines.append(
                        f"  ↳ {moderation_type} ({DurationObject.from_seconds(duration_seconds)})"
                    )
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
                if lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                lines = []
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

        await StateService.send_pages(obj=Cap, pages=pages, state=state)

    # DONE
    @app_commands.command(name="cmds", description="List aliases.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_commands_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=interaction)
        title = f"{get_random_emoji()} Command Aliases"

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        aliases = await Alias.select(**kwargs)

        for alias in aliases:
            guild_dictionary.setdefault(alias.guild_snowflake, {"channels": {}})
            guild_dictionary[alias.guild_snowflake]["channels"].setdefault(
                alias.channel_snowflake, {"aliases": {}}
            )
            guild_dictionary[alias.guild_snowflake]["channels"][
                alias.channel_snowflake
            ]["aliases"].setdefault(alias.alias_type, []).append(alias.alias_name)

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
            for channel_snowflake, dictionary in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                channel_lines = []
                for alias_type, alias_names in dictionary["aliases"].items():
                    channel_lines.append(f"{alias_type}")
                    for name in alias_names:
                        channel_lines.append(f"  ↳ {name}")
                if len(lines) >= chunk_size:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
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
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )

        await StateService.send_pages(obj=Alias, pages=pages, state=state)

    # DONE
    @commands.command(name="cmds", help="List aliases.")
    @moderator_predicator()
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=ctx)
        title = f"{get_random_emoji()} Command Aliases"

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        aliases = await Alias.select(**kwargs)

        for alias in aliases:
            guild_dictionary.setdefault(alias.guild_snowflake, {"channels": {}})
            guild_dictionary[alias.guild_snowflake]["channels"].setdefault(
                alias.channel_snowflake, {"aliases": {}}
            )
            guild_dictionary[alias.guild_snowflake]["channels"][
                alias.channel_snowflake
            ]["aliases"].setdefault(alias.alias_type, []).append(alias.alias_name)

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
            for channel_snowflake, dictionary in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                channel_lines = []
                for alias_type, alias_names in dictionary["aliases"].items():
                    channel_lines.append(f"{alias_type}")
                    for name in alias_names:
                        channel_lines.append(f"  ↳ {name}")
                if len(lines) >= chunk_size:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
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
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )

        await StateService.send_pages(obj=Alias, pages=pages, state=state)

    @app_commands.command(name="coords", description="Lists coords.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_coordinators_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Coordinator.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        coordinators = await Coordinator.select(**kwargs)

        for coordinator in coordinators:
            guild_dictionary.setdefault(coordinator.guild_snowflake, {"members": {}})
            guild_dictionary[coordinator.guild_snowflake]["members"].setdefault(
                coordinator.member_snowflake, {"coordinators": {}}
            )
            guild_dictionary[coordinator.guild_snowflake]["members"][
                coordinator.member_snowflake
            ]["coordinators"].setdefault(coordinator.channel_snowflake, {})
            guild_dictionary[coordinator.guild_snowflake]["members"][
                coordinator.member_snowflake
            ]["coordinators"][coordinator.channel_snowflake].update(
                {"placeholder": "placeholder"}
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, coordinator_dictionary in guild_data.get(
                "members"
            ).items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in coordinator_dictionary.get(
                    "coordinators"
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Coordinator, pages=pages, state=state)

    # DONE
    @commands.command(name="coords", help="Lists coords.")
    @moderator_predicator()
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            description="Specify one of: `all`, channel ID/mention, or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Coordinator.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        coordinators = await Coordinator.select(**kwargs)

        for coordinator in coordinators:
            guild_dictionary.setdefault(coordinator.guild_snowflake, {"members": {}})
            guild_dictionary[coordinator.guild_snowflake]["members"].setdefault(
                coordinator.member_snowflake, {"coordinators": {}}
            )
            guild_dictionary[coordinator.guild_snowflake]["members"][
                coordinator.member_snowflake
            ]["coordinators"].setdefault(coordinator.channel_snowflake, {})
            guild_dictionary[coordinator.guild_snowflake]["members"][
                coordinator.member_snowflake
            ]["coordinators"][coordinator.channel_snowflake].update(
                {"placeholder": "placeholder"}
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, coordinator_dictionary in guild_data.get(
                "members"
            ).items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in coordinator_dictionary.get(
                    "coordinators"
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Coordinator, pages=pages, state=state)

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
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Developer.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        developers = await Developer.select(**kwargs)

        for developer in developers:
            guild_dictionary.setdefault(developer.guild_snowflake, {"members": {}})
            guild_dictionary[developer.guild_snowflake]["members"].setdefault(
                developer.member_snowflake, {"developers": {}}
            )
            guild_dictionary[developer.guild_snowflake]["members"][
                developer.member_snowflake
            ]["developers"].update({"placeholder": "placeholder"})

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            lines = []
            thumbnail = False
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, developer_dictionary in guild_data.get(
                "members"
            ).items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
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
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Channels in Server",
                )

        await StateService.send_pages(obj=Developer, pages=pages, state=state)

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
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Developer.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        developers = await Developer.select(**kwargs)

        for developer in developers:
            guild_dictionary.setdefault(developer.guild_snowflake, {"members": {}})
            guild_dictionary[developer.guild_snowflake]["members"].setdefault(
                developer.member_snowflake, {"developers": {}}
            )
            guild_dictionary[developer.guild_snowflake]["members"][
                developer.member_snowflake
            ]["developers"].update({"placeholder": "placeholder"})

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            lines = []
            thumbnail = False
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, developer_dictionary in guild_data.get(
                "members"
            ).items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
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
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Channels in Server",
                )

        await StateService.send_pages(obj=Developer, pages=pages, state=state)

    # DONE
    @app_commands.command(name="flags", description="List flags.")
    @moderator_predicator()
    async def list_flags_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Flag.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        flags = await Flag.select(**kwargs)

        for flag in flags:
            guild_dictionary.setdefault(flag.guild_snowflake, {"members": {}})
            guild_dictionary[flag.guild_snowflake]["members"].setdefault(
                flag.member_snowflake, {"flags": {}}
            )
            guild_dictionary[flag.guild_snowflake]["members"][flag.member_snowflake][
                "flags"
            ].setdefault(flag.channel_snowflake, {})
            guild_dictionary[flag.guild_snowflake]["members"][flag.member_snowflake][
                "flags"
            ][flag.channel_snowflake].update(
                {
                    "reason": flag.reason,
                }
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, flag_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in flag_dictionary.get(
                    "flags", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
                    if object_dict.get("type", None) == discord.Member:
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Flag, pages=pages, state=state)

    # DONE
    @commands.command(name="flags", help="List flags.")
    @moderator_predicator()
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Flag.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        flags = await Flag.select(**kwargs)

        for flag in flags:
            guild_dictionary.setdefault(flag.guild_snowflake, {"members": {}})
            guild_dictionary[flag.guild_snowflake]["members"].setdefault(
                flag.member_snowflake, {"flags": {}}
            )
            guild_dictionary[flag.guild_snowflake]["members"][flag.member_snowflake][
                "flags"
            ].setdefault(flag.channel_snowflake, {})
            guild_dictionary[flag.guild_snowflake]["members"][flag.member_snowflake][
                "flags"
            ][flag.channel_snowflake].update(
                {
                    "reason": flag.reason,
                }
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, flag_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in flag_dictionary.get(
                    "flags", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
                    if object_dict.get("type", None) == discord.Member:
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Flag, pages=pages, state=state)

    # DONE
    @app_commands.command(name="ls", description="List new vegans.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_new_vegans_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Vegan.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        vegans = await Vegan.select(**kwargs)

        for vegan in vegans:
            guild_dictionary.setdefault(vegan.guild_snowflake, {"members": {}})
            guild_dictionary[vegan.guild_snowflake]["members"].setdefault(
                vegan.member_snowflake, {"vegans": {}}
            )
            guild_dictionary[vegan.guild_snowflake]["members"][vegan.member_snowflake][
                "vegans"
            ].setdefault("placeholder", {})

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, vegan_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                else:
                    if not thumbnail:
                        embed.set_thumbnail(
                            url=object_dict.get("object", None).display_avatar.url
                        )
                        thumbnail = True
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Vegan, pages=pages, state=state)

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
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Vegan.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        vegans = await Vegan.select(**kwargs)

        for vegan in vegans:
            guild_dictionary.setdefault(vegan.guild_snowflake, {"members": {}})
            guild_dictionary[vegan.guild_snowflake]["members"].setdefault(
                vegan.member_snowflake, {"vegans": {}}
            )
            guild_dictionary[vegan.guild_snowflake]["members"][vegan.member_snowflake][
                "vegans"
            ].setdefault("placeholder", {})

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, vegan_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                else:
                    if not thumbnail:
                        embed.set_thumbnail(
                            url=object_dict.get("object", None).display_avatar.url
                        )
                        thumbnail = True
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Vegan, pages=pages, state=state)

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
                "room_name": channel_dict["name"],
            }
            where_kwargs = {
                "channel_snowflake": old_room.channel_snowflake,
                "guild_snowflake": interaction.guild.id,
            }
            await TemporaryRoom.update(
                set_kwargs=set_kwargs,
                where_kwargs=temp_where_kwargs,
            )
            await Alias.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Ban.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Cap.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Coordinator.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Flag.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Moderator.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Stage.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await TextMute.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Vegan.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await VoiceMute.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            return await state.end(
                success=f"Temporary room `{old_name}` migrated to {channel_dict['mention']}."
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
                "room_name": channel_dict["name"],
            }
            where_kwargs = {
                "channel_snowflake": old_room.channel_snowflake,
                "guild_snowflake": ctx.guild.id,
            }
            await TemporaryRoom.update(
                set_kwargs=set_kwargs,
                where_kwargs=temp_where_kwargs,
            )
            await Alias.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Ban.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Cap.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Coordinator.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Flag.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Moderator.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Stage.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await TextMute.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await Vegan.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            await VoiceMute.update(
                set_kwargs=set_kwargs,
                where_kwargs=where_kwargs,
            )
            return await state.end(
                success=f"Temporary room `{old_name}` migrated to {channel_dict['mention']}."
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
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Moderator.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        moderators = await Moderator.select(**kwargs)

        for moderator in moderators:
            guild_dictionary.setdefault(moderator.guild_snowflake, {"members": {}})
            guild_dictionary[moderator.guild_snowflake]["members"].setdefault(
                moderator.member_snowflake, {"moderators": {}}
            )
            guild_dictionary[moderator.guild_snowflake]["members"][
                moderator.member_snowflake
            ]["moderators"].setdefault(moderator.channel_snowflake, {})
            guild_dictionary[moderator.guild_snowflake]["members"][
                moderator.member_snowflake
            ]["moderators"][moderator.channel_snowflake].update(
                {"placeholder": "placeholder"}
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, member_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in member_dictionary.get(
                    "moderators", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Moderator, pages=pages, state=state)

    # DONE
    @commands.command(name="mods", help="Lists mods.")
    @moderator_predicator()
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {Moderator.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        moderators = await Moderator.select(**kwargs)

        for moderator in moderators:
            guild_dictionary.setdefault(moderator.guild_snowflake, {"members": {}})
            guild_dictionary[moderator.guild_snowflake]["members"].setdefault(
                moderator.member_snowflake, {"moderators": {}}
            )
            guild_dictionary[moderator.guild_snowflake]["members"][
                moderator.member_snowflake
            ]["moderators"].setdefault(moderator.channel_snowflake, {})
            guild_dictionary[moderator.guild_snowflake]["members"][
                moderator.member_snowflake
            ]["moderators"][moderator.channel_snowflake].update(
                {"placeholder": "placeholder"}
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, member_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in member_dictionary.get(
                    "moderators", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
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
        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=Moderator, pages=pages, state=state)

    # DONE
    @app_commands.command(name="mutes", description="List mutes.")
    @moderator_predicator()
    async def list_mutes_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {VoiceMute.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        voice_mutes = await VoiceMute.select(target="user", **kwargs)

        for voice_mute in voice_mutes:
            guild_dictionary.setdefault(voice_mute.guild_snowflake, {"members": {}})
            guild_dictionary[voice_mute.guild_snowflake]["members"].setdefault(
                voice_mute.member_snowflake, {"voice_mutes": {}}
            )
            guild_dictionary[voice_mute.guild_snowflake]["members"][
                voice_mute.member_snowflake
            ]["voice_mutes"].setdefault(voice_mute.channel_snowflake, {})
            guild_dictionary[voice_mute.guild_snowflake]["members"][
                voice_mute.member_snowflake
            ]["voice_mutes"][voice_mute.channel_snowflake].update(
                {
                    "reason": voice_mute.reason,
                    "expires_in": DurationObject.from_expires_in(voice_mute.expires_in),
                }
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, voice_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in voice_mute_dictionary.get(
                    "voice_mutes", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
                    lines.append(f"**Expires in:** {channel_dictionary['expires_in']}")
                    if object_dict.get("type", None) == discord.Member:
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=VoiceMute, pages=pages, state=state)

    # DONE
    @commands.command(name="mutes", help="List mutes.")
    @moderator_predicator()
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        title = f"{get_random_emoji()} {VoiceMute.PLURAL} {f'for {object_dict.get('name', None)}' if object_dict.get("type", None) == discord.Member else ''}"
        kwargs = object_dict["columns"]

        voice_mutes = await VoiceMute.select(target="user", **kwargs)

        for voice_mute in voice_mutes:
            guild_dictionary.setdefault(voice_mute.guild_snowflake, {"members": {}})
            guild_dictionary[voice_mute.guild_snowflake]["members"].setdefault(
                voice_mute.member_snowflake, {"voice_mutes": {}}
            )
            guild_dictionary[voice_mute.guild_snowflake]["members"][
                voice_mute.member_snowflake
            ]["voice_mutes"].setdefault(voice_mute.channel_snowflake, {})
            guild_dictionary[voice_mute.guild_snowflake]["members"][
                voice_mute.member_snowflake
            ]["voice_mutes"][voice_mute.channel_snowflake].update(
                {
                    "reason": voice_mute.reason,
                    "expires_in": DurationObject.from_expires_in(voice_mute.expires_in),
                }
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, voice_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in voice_mute_dictionary.get(
                    "voice_mutes", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
                    lines.append(f"**Expires in:** {channel_dictionary['expires_in']}")
                    if object_dict.get("type", None) == discord.Member:
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=VoiceMute, pages=pages, state=state)

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
        await has_equal_or_lower_role(
            source=interaction,
            member_snowflake=member_dict["id"],
            sender_snowflake=interaction.user.id,
        )
        kwargs = channel_dict["columns"]

        stage = await Stage.select(**kwargs)
        if not stage:
            return await state.end(
                warning=f"No active stage found in {channel_dict['mention']}."
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
        await has_equal_or_lower_role(
            source=ctx,
            member_snowflake=member_dict["id"],
            sender_snowflake=ctx.author.id,
        )
        kwargs = channel_dict["columns"]

        stage = await Stage.select(**kwargs)
        if not stage:
            return await state.end(
                warning=f"No active stage found in {channel_dict['mention']}."
            )
        try:
            await member_dict.get("object", None).edit(
                mute=not member_dict.get("object", None).voice.mute
            )
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())

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
    async def get_role_id(self, ctx: commands.Context, *, role_name: str):
        state = StateService(source=ctx)
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            return await state.end(success=f"Role `{role.name}` has ID `{role.id}`.")
        else:
            return await state.end(
                warning=f"No role named `{role_name}` found in this server."
            )

    # DONE
    @app_commands.command(name="stages", description="List stages.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_stages_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=interaction)
        title = f"{get_random_emoji()} Stages"

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        stages = await Stage.select(**kwargs)

        for stage in stages:
            guild_dictionary.setdefault(stage.guild_snowflake, {"channels": {}})
            guild_dictionary[stage.guild_snowflake]["channels"].setdefault(
                stage.channel_snowflake, {}
            )
            guild_dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ].setdefault("stages", {})
            guild_dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ]["stages"].update(
                {"expires_in": DurationObject.from_expires_in(stage.expires_in)}
            )

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
            for channel_snowflake, stage_dictionary in guild_data.get(
                "channels"
            ).items():
                channel = guild.get_channel(channel_snowflake)
                for channel_snowflake in stage_dictionary.get("stages", {}).items():
                    lines.append(f"**Expires in:** {stage_dictionary['expires_in']}")
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed, field_count = flush_page(embed, pages, title, guild.name)
                        lines = []
                    if lines:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
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
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )

        await StateService.send_pages(obj=Stage, pages=pages, state=state)

    # DONE
    @commands.command(name="stages", help="List stages.")
    @moderator_predicator()
    async def list_stages_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=ctx)
        title = f"{get_random_emoji()} Stages"

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        stages = await Stage.select(**kwargs)

        for stage in stages:
            guild_dictionary.setdefault(stage.guild_snowflake, {"channels": {}})
            guild_dictionary[stage.guild_snowflake]["channels"].setdefault(
                stage.channel_snowflake, {}
            )
            guild_dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ].setdefault("stages", {})
            guild_dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ]["stages"].update(
                {"expires_in": DurationObject.from_expires_in(stage.expires_in)}
            )

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
            for channel_snowflake, stage_dictionary in guild_data.get(
                "channels"
            ).items():
                channel = guild.get_channel(channel_snowflake)
                for channel_snowflake in stage_dictionary.get("stages", {}).items():
                    lines.append(f"**Expires in:** {stage_dictionary['expires_in']}")
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed, field_count = flush_page(embed, pages, title, guild.name)
                        lines = []
                    if lines:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
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
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )

        await StateService.send_pages(obj=Stage, pages=pages, state=state)

    @app_commands.command(name="summary", description="Moderation summary.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
        target="Specify 'all' or a channel ID/mention.",
    )
    @moderator_predicator()
    async def list_moderation_summary_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        target: Optional[str] = None,
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        objects_list = []

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        member_dict = await do.determine_from_target(target=member)
        await has_equal_or_lower_role(
            source=interaction,
            member_snowflake=member_dict["id"],
            sender_snowflake=interaction.user.id,
        )
        kwargs = member_dict["columns"]

        title = f"{get_random_emoji()} Summary for {member_dict.get('name', None)}"
        for name, object in member_relevant_objects_dict.items():
            if target and target.lower() == name or target.lower() == "all":
                objects = object.select(**kwargs)
                objects_list = [*objects]

        guild = self.bot.get_guild(interaction.guild.id)

        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for obj in objects_list:
            channel = guild.get_channel(obj.channel_snowflake)
            entries = [f"**User:** {member_dict['mention']}"]
            line = entries[0]
            if hasattr(obj, "expires_in"):
                entries.append(f"**Expires in:** {obj.expires_in}")
                line = "\n".join(entries)
            if hasattr(obj, "reason"):
                entries.append(f"**Reason:** {obj.reason}")
                line = "\n".join(entries)
            lines.append(line)
            field_count += 1
            if field_count >= chunk_size:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed.set_thumbnail(
                    url=member_dict.get("object", None).display_avatar.url
                )
                embed, field_count = flush_page(embed, pages, title, guild.name)
                lines = []
        if lines:
            embed.add_field(
                name=f"Channel: {channel.mention}",
                value="\n".join(lines),
                inline=False,
            )
            embed.set_thumbnail(url=member_dict.get("object", None).display_avatar.url)
            pages.append(embed)
            lines = []

        await StateService.send_pages(obj=Alias, pages=pages, state=state)

    @commands.command(name="summary", description="Moderation summary.")
    @moderator_predicator()
    async def list_moderation_summary_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            description="Specify a member ID/mention."
        ),
        target: Optional[str] = commands.parameter(
            description="Specify 'all' or a channel ID/mention."
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        member_dict = await do.determine_from_target(target=member)
        await has_equal_or_lower_role(
            source=ctx,
            member_snowflake=member_dict["id"],
            sender_snowflake=ctx.author.id,
        )
        kwargs = member_dict["columns"]

        title = f"{get_random_emoji()} Summary for {member_dict.get('name', None)}"
        for name, object in member_relevant_objects_dict.items():
            if target and target.lower() == name or target.lower() == "all":
                objects = object.select(**kwargs)
                objects_list = [*objects]

        guild = self.bot.get_guild(ctx.guild.id)

        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for obj in objects:
            channel = guild.get_channel(obj.channel_snowflake)
            entries = [f"**User:** {member_dict['mention']}"]
            line = entries[0]
            if hasattr(obj, "expires_in"):
                entries.append(f"**Expires in:** {obj.expires_in}")
                line = "\n".join(entries)
            if hasattr(obj, "reason"):
                entries.append(f"**Reason:** {obj.reason}")
                line = "\n".join(entries)
            lines.append(line)
            field_count += 1
            if field_count >= chunk_size:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed.set_thumbnail(
                    url=member_dict.get("object", None).display_avatar.url
                )
                embed, field_count = flush_page(embed, pages, title, guild.name)
                lines = []
        if lines:
            embed.add_field(
                name=f"Channel: {channel.mention}",
                value="\n".join(lines),
                inline=False,
            )
            embed.set_thumbnail(url=member_dict.get("object", None).display_avatar.url)
            pages.append(embed)
            lines = []

        await StateService.send_pages(obj=Alias, pages=pages, state=state)

    @app_commands.command(name="survey", description="Get all.")
    @app_commands.describe(channel="Tag a voice/stage channel")
    @moderator_predicator()
    async def stage_survey_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        chunk_size, pages = 7, []
        system_owners = developers = guild_owners = administrators = coordinators = (
            moderators
        ) = []

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        channel_dict = await do.determine_from_target(target=channel)

        for member in channel_dict.get("object", None).members:
            try:
                if await member_is_system_owner(member.id):
                    system_owners.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await member_is_developer(interaction.guild.id, member.id):
                    developers.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await member_is_guild_owner(interaction.guild.id, member.id):
                    guild_owners.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await member_is_administrator(interaction.guild.id, member.id):
                    administrators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await member_is_coordinator(
                    channel_dict.get("id", None), interaction.guild.id, member.id
                ):
                    coordinators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await member_is_moderator(
                    channel_dict.get("id", None), interaction.guild.id, member.id
                ):
                    moderators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
        system_owners_chunks = [
            system_owners[i : i + chunk_size]
            for i in range(0, len(system_owners), chunk_size)
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
            ("System Owners", system_owners, system_owners_chunks),
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

        await StateService.send_pages(obj=PermissionRole, pages=pages, state=state)

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
        system_owners = developers = guild_owners = administrators = coordinators = (
            moderators
        ) = []

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        channel_dict = await do.determine_from_target(target=channel)

        for member in channel_dict.get("object", None).members:
            try:
                if await member_is_system_owner(member.id):
                    system_owners.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await member_is_developer(ctx.guild.id, member.id):
                    developers.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await member_is_guild_owner(ctx.guild.id, member.id):
                    guild_owners.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await member_is_administrator(ctx.guild.id, member.id):
                    administrators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await member_is_coordinator(
                    channel_dict.get("id", None), ctx.guild.id, member.id
                ):
                    coordinators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
            try:
                if await member_is_moderator(
                    channel_dict.get("id", None), ctx.guild.id, member.id
                ):
                    moderators.append(member)
            except commands.CheckFailure as e:
                logger.warning(str(e).capitalize())
        system_owners_chunks = [
            system_owners[i : i + chunk_size]
            for i in range(0, len(system_owners), chunk_size)
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
            ("System Owners", system_owners, system_owners_chunks),
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

        await StateService.send_pages(obj=PermissionRole, pages=pages, state=state)

    # DONE
    @app_commands.command(name="tmutes", description="List text-mutes.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, member ID/mention, or server ID."
    )
    @moderator_predicator()
    async def list_text_mutes_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        title = f"{get_random_emoji()} {TextMute.PLURAL}"
        is_at_home = at_home(source=interaction)

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        text_mutes = await TextMute.select(**kwargs)

        for text_mute in text_mutes:
            guild_dictionary.setdefault(text_mute.guild_snowflake, {"members": {}})
            guild_dictionary[text_mute.guild_snowflake]["members"].setdefault(
                text_mute.member_snowflake, {"text_mutes": {}}
            )
            guild_dictionary[text_mute.guild_snowflake]["members"][
                text_mute.member_snowflake
            ]["text_mutes"].setdefault(text_mute.channel_snowflake, {})
            guild_dictionary[text_mute.guild_snowflake]["members"][
                text_mute.member_snowflake
            ]["text_mutes"][text_mute.channel_snowflake].update(
                {
                    "reason": text_mute.reason,
                    "expires_in": DurationObject.from_expires_in(text_mute.expires_in),
                }
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, text_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in text_mute_dictionary.get(
                    "text_mutes", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
                    lines.append(f"**Expires in:** {channel_dictionary['expires_in']}")
                    if object_dict.get("type", None) == discord.Member:
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=TextMute, pages=pages, state=state)

    # DONE
    @commands.command(name="tmutes", help="List text-mutes.")
    @moderator_predicator()
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            description="Specify one of: 'all', channel ID/mention, or server ID.",
        ),
    ):
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        title = f"{get_random_emoji()} {TextMute.PLURAL}"
        is_at_home = at_home(source=ctx)

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        object_dict = await do.determine_from_target(target=target)
        kwargs = object_dict["columns"]

        text_mutes = await TextMute.select(**kwargs)

        for text_mute in text_mutes:
            guild_dictionary.setdefault(text_mute.guild_snowflake, {"members": {}})
            guild_dictionary[text_mute.guild_snowflake]["members"].setdefault(
                text_mute.member_snowflake, {"text_mutes": {}}
            )
            guild_dictionary[text_mute.guild_snowflake]["members"][
                text_mute.member_snowflake
            ]["text_mutes"].setdefault(text_mute.channel_snowflake, {})
            guild_dictionary[text_mute.guild_snowflake]["members"][
                text_mute.member_snowflake
            ]["text_mutes"][text_mute.channel_snowflake].update(
                {
                    "reason": text_mute.reason,
                    "expires_in": DurationObject.from_expires_in(text_mute.expires_in),
                }
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, text_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if object_dict.get("type", None) != discord.Member:
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in text_mute_dictionary.get(
                    "text_mutes", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
                    lines.append(f"**Expires in:** {channel_dictionary['expires_in']}")
                    if object_dict.get("type", None) == discord.Member:
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
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

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )

        await StateService.send_pages(obj=TextMute, pages=pages, state=state)


async def setup(bot: DiscordBot):
    cog = ModeratorCommands(bot)
    await bot.add_cog(cog)
