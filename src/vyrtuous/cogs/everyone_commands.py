"""everyone_commands.py A discord.py cog containing everyone commands for the Vyrtuous bot.

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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.roles.administrator import Administrator
from vyrtuous.database.roles.coordinator import Coordinator
from vyrtuous.database.roles.developer import Developer
from vyrtuous.database.roles.moderator import Moderator
from vyrtuous.properties.snowflake import (
    AppChannelSnowflake,
    ChannelSnowflake,
)
from vyrtuous.service.check_service import (
    at_home,
    member_is_administrator,
    member_is_coordinator,
    member_is_developer,
    member_is_guild_owner,
    member_is_moderator,
    member_is_system_owner,
)
from vyrtuous.service.logging_service import logger
from vyrtuous.service.resolution.channel_service import resolve_channel
from vyrtuous.service.resolution.member_service import resolve_member
from vyrtuous.service.messaging.message_service import MessageService
from vyrtuous.service.scope_service import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    resolve_objects,
    generate_skipped_channels,
    generate_skipped_guilds,
    generate_skipped_roles,
    generate_skipped_members,
    clean_guild_dictionary,
    flush_page,
)
from vyrtuous.service.messaging.state_service import StateService
from vyrtuous.utils.emojis import get_random_emoji


class EveryoneCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.message_service = MessageService(self.bot, self.bot.db_pool)

    @app_commands.command(name="admins", description="Lists admins.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, " "server ID or empty."
    )
    async def list_administrators_app_command(
        self, interaction: discord.Interaction, target: str
    ):
        state = StateService(interaction)
        administrators, title = await resolve_objects(
            ctx_interaction_or_message=interaction,
            obj=Administrator,
            state=state,
            target=target,
        )

        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}

        for administrator in administrators:
            guild_dictionary.setdefault(administrator.guild_snowflake, {"channels": {}, "members": {}})
            guild_dictionary[administrator.guild_snowflake]["members"][administrator.member_snowflake] = administrator.role_snowflakes

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        skipped_roles = await generate_skipped_roles(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
            skipped_roles=skipped_roles,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, role_snowflakes in guild_data.get("members").items():
                if field_count >= chunk_size:
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                member = guild.get_member(member_snowflake)
                role_mentions = [guild.get_role(role_snowflake).mention for role_snowflake in role_snowflakes if guild.get_role(role_snowflake)]
                if role_mentions:
                    lines = f"**Roles:** {'\n'.join(role_mentions)}"
                else:
                    lines = "\u200b"
                embed.add_field(
                    name=f"{member.name} ({member_snowflake})",
                    value=lines,
                    inline=False,
                )
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
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )
            if skipped_roles:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_roles,
                    title="Skipped Roles in Server",
                )

        await StateService.send_pages(obj=Administrator, pages=pages, state=state)

    # DONE
    @commands.command(name="admins", help="Lists admins.")
    async def list_administrators_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: `all`, "
            "channel ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        administrators, title = await resolve_objects(
            ctx_interaction_or_message=ctx,
            obj=Administrator,
            state=state,
            target=target,
        )

        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}

        for administrator in administrators:
            guild_dictionary.setdefault(administrator.guild_snowflake, {"channels": {}, "members": {}})
            guild_dictionary[administrator.guild_snowflake]["members"][administrator.member_snowflake] = administrator.role_snowflakes

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        skipped_roles = await generate_skipped_roles(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
            skipped_roles=skipped_roles,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, role_snowflakes in guild_data.get("members").items():
                if field_count >= chunk_size:
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                member = guild.get_member(member_snowflake)
                role_mentions = [guild.get_role(role_snowflake).mention for role_snowflake in role_snowflakes if guild.get_role(role_snowflake)]
                if role_mentions:
                    lines = f"**Roles:** {'\n'.join(role_mentions)}"
                else:
                    lines = "\u200b"
                embed.add_field(
                    name=f"{member.name} ({member_snowflake})",
                    value=lines,
                    inline=False,
                )
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
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )
            if skipped_roles:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_roles,
                    title="Skipped Roles in Server",
                )

        await StateService.send_pages(obj=Administrator, pages=pages, state=state)

    # DONE
    @app_commands.command(name="coords", description="Lists coords.")
    @app_commands.describe(
        target="Specify one of: `all`, channel ID/mention, server ID or empty."
    )
    async def list_coordinators_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        coordinators, title = await resolve_objects(
            ctx_interaction_or_message=interaction,
            obj=Coordinator,
            state=state,
            target=target,
        )

        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(
                ctx_interaction_or_message=interaction, member_str=target
            )
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        channel_lines, chunk_size, field_count, lines, pages = [], 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False

        for coordinator in coordinators:
            guild_dictionary.setdefault(coordinator.guild_snowflake, {})
            guild_dictionary[coordinator.guild_snowflake].setdefault(
                coordinator.channel_snowflake, []
            )
            guild_dictionary[coordinator.guild_snowflake][
                coordinator.channel_snowflake
            ].append({"member_snowflake": coordinator.member_snowflake})

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                lines = []
                for member_data in members:
                    member = guild.get_member(member_data["member_snowflake"])
                    if not member_obj:
                        lines.append(f"**User**: {member.mention}")
                    else:
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
                field_count += 1
                if field_count == chunk_size:
                    if lines:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        lines = []
                    elif channel_lines:
                        embed.add_field(
                            name="Channels",
                            value="\n".join(channel_lines),
                            inline=False,
                        )
                        channel_lines = []
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                if lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                else:
                    channel_lines.append(channel.mention)
            if channel_lines:
                embed.add_field(
                    name="Channels", value="\n".join(channel_lines), inline=False
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
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: `all`, channel ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        coordinators, title = await resolve_objects(
            ctx_interaction_or_message=ctx, obj=Coordinator, state=state, target=target
        )

        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(
                ctx_interaction_or_message=ctx, member_str=target
            )
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        channel_lines, chunk_size, field_count, lines, pages = [], 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False

        for coordinator in coordinators:
            guild_dictionary.setdefault(coordinator.guild_snowflake, {})
            guild_dictionary[coordinator.guild_snowflake].setdefault(
                coordinator.channel_snowflake, []
            )
            guild_dictionary[coordinator.guild_snowflake][
                coordinator.channel_snowflake
            ].append({"member_snowflake": coordinator.member_snowflake})

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, channels in guild_dictionary.items():
            channel_lines = []
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                lines = []
                for member_data in members:
                    member = guild.get_member(member_data["member_snowflake"])
                    if member_obj:
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
                    else:
                        lines.append(f"**User**: {member.mention}")
                if lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    field_count += 1
                else:
                    channel_lines.append(channel.mention)
                if field_count == chunk_size:
                    if channel_lines:
                        embed.add_field(
                            name="Channels",
                            value="\n".join(channel_lines),
                            inline=False,
                        )
                        channel_lines.clear()
                    embed, field_count = flush_page(embed, pages, title, guild.name)
            if channel_lines:
                embed.add_field(
                    name="Channels",
                    value="\n".join(channel_lines),
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
    @app_commands.command(name="devs", description="List devs.")
    @app_commands.describe(target="Specify one of: 'all', server ID or empty.")
    async def list_developers_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        developers, title = await resolve_objects(
            ctx_interaction_or_message=interaction,
            obj=Developer,
            state=state,
            target=target,
        )

        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(
                ctx_interaction_or_message=interaction, member_str=target
            )
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}
        thumbnail = False

        for developer in developers:
            guild_dictionary.setdefault(developer.guild_snowflake, {})
            guild_dictionary[developer.guild_snowflake].setdefault(
                developer.member_snowflake, []
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, member_ids in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            if member_obj:
                focused_member = guild.get_member(member_obj.id)
            else:
                focused_member = None
            thumbnail = False
            members = (
                [focused_member]
                if focused_member
                else [guild.get_member(mid) for mid in member_ids]
            )
            for member in members:
                if not member:
                    continue
                if field_count >= chunk_size:
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                if focused_member and not thumbnail:
                    embed.set_thumbnail(url=member.display_avatar.url)
                    thumbnail = True
                embed.add_field(
                    name=f"{member.name} ({member.id})",
                    value=member.mention,
                    inline=False,
                )
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
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None, description="'all', a specific server or user mention/ID"
        ),
    ):
        state = StateService(ctx)
        developers, title = await resolve_objects(
            ctx_interaction_or_message=ctx, obj=Developer, state=state, target=target
        )

        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(
                ctx_interaction_or_message=ctx, member_str=target
            )
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}

        for developer in developers:
            guild_dictionary.setdefault(developer.guild_snowflake, {})
            guild_dictionary[developer.guild_snowflake].setdefault(
                developer.member_snowflake, []
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, member_ids in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            if member_obj:
                focused_member = guild.get_member(member_obj.id)
            else:
                focused_member = None
            thumbnail = False
            members = (
                [focused_member]
                if focused_member
                else [guild.get_member(mid) for mid in member_ids]
            )
            for member in members:
                if not member:
                    continue
                if field_count >= chunk_size:
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                if focused_member and not thumbnail:
                    embed.set_thumbnail(url=member.display_avatar.url)
                    thumbnail = True
                embed.add_field(
                    name=f"{member.name} ({member.id})",
                    value=member.mention,
                    inline=False,
                )
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
    @app_commands.command(name="mods", description="Lists mods.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, server ID or empty."
    )
    async def list_moderators_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        moderators, title = await resolve_objects(
            ctx_interaction_or_message=interaction,
            obj=Moderator,
            state=state,
            target=target,
        )

        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(
                ctx_interaction_or_message=interaction, member_str=target
            )
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        channel_lines, chunk_size, field_count, lines, pages = [], 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False

        for moderator in moderators:
            guild_dictionary.setdefault(moderator.guild_snowflake, {})
            guild_dictionary[moderator.guild_snowflake].setdefault(
                moderator.channel_snowflake, []
            )
            guild_dictionary[moderator.guild_snowflake][
                moderator.channel_snowflake
            ].append({"member_snowflake": moderator.member_snowflake})

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                lines = []
                for member_data in members:
                    member = guild.get_member(member_data["member_snowflake"])
                    if not member_obj:
                        lines.append(f"**User**: {member.mention}")
                    else:
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
                field_count += 1
                if field_count == chunk_size:
                    if lines:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        lines = []
                    elif channel_lines:
                        embed.add_field(
                            name="Channels",
                            value="\n".join(channel_lines),
                            inline=False,
                        )
                        channel_lines = []
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                if lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                else:
                    channel_lines.append(f"{channel.mention}")
            if channel_lines:
                embed.add_field(
                    name="Channels", value="\n".join(channel_lines), inline=False
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
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        moderators, title = await resolve_objects(
            ctx_interaction_or_message=ctx, obj=Moderator, state=state, target=target
        )

        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(
                ctx_interaction_or_message=ctx, member_str=target
            )
        except Exception as e:
            logger.warning(str(e).capitalize())
            pass

        channel_lines, chunk_size, field_count, lines, pages = [], 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False

        for moderator in moderators:
            guild_dictionary.setdefault(moderator.guild_snowflake, {})
            guild_dictionary[moderator.guild_snowflake].setdefault(
                moderator.channel_snowflake, []
            )
            guild_dictionary[moderator.guild_snowflake][
                moderator.channel_snowflake
            ].append({"member_snowflake": moderator.member_snowflake})

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                lines = []
                for member_data in members:
                    member = guild.get_member(member_data["member_snowflake"])
                    if not member_obj:
                        lines.append(f"**User**: {member.mention}")
                    else:
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
                field_count += 1
                if field_count == chunk_size:
                    if lines:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        lines = []
                    elif channel_lines:
                        embed.add_field(
                            name="Channels",
                            value="\n".join(channel_lines),
                            inline=False,
                        )
                        channel_lines = []
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                if lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                else:
                    channel_lines.append(f"{channel.mention}")
            if channel_lines:
                embed.add_field(
                    name="Channels", value="\n".join(channel_lines), inline=False
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
    @app_commands.command(name="roleid", description="Get role by name.")
    @app_commands.describe(role_name="The name of the role to look up")
    async def get_role_id_app_command(
        self, interaction: discord.Interaction, role_name: str
    ):
        state = StateService(interaction)
        role = discord.utils.get(interaction.guild.roles, name=role_name)
        if role:
            try:
                return await state.end(
                    success=f"{get_random_emoji()} Role `{role.name}` has ID `{role.id}`."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        else:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f No role named `{role_name}` found in this server."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @commands.command(name="roleid", help="Get role by name.")
    async def get_role_id(self, ctx: commands.Context, *, role_name: str):
        state = StateService(ctx)
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            try:
                return await state.end(
                    success=f"{get_random_emoji()} Role `{role.name}` has ID `{role.id}`."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        else:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f No role named `{role_name}` found in this server."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @app_commands.command(name="survey", description="Get all.")
    @app_commands.describe(channel="Tag a voice/stage channel")
    async def stage_survey_app_command(
        self, interaction: discord.Interaction, channel: AppChannelSnowflake
    ):
        state = StateService(interaction)
        channel_obj = None
        chunk_size = 7
        pages = []
        (
            system_owners,
            developers,
            guild_owners,
            administrators,
            moderators,
            coordinators,
        ) = ([], [], [], [], [], [])
        try:
            channel_obj = await resolve_channel(interaction, channel)
        except Exception as e:
            channel_obj = interaction.channel
            logger.warning(str(e).capitalize())
        for member in channel_obj.members:
            try:
                if await member_is_system_owner(member.id):
                    system_owners.append(member)
            except commands.CheckFailure:
                pass
            try:
                if await member_is_developer(interaction.guild.id, member.id):
                    developers.append(member)
            except commands.CheckFailure:
                pass
            try:
                if await member_is_guild_owner(interaction.guild.id, member.id):
                    guild_owners.append(member)
            except commands.CheckFailure:
                pass
            try:
                if await member_is_administrator(interaction.guild.id, member.id):
                    administrators.append(member)
            except commands.CheckFailure:
                pass
            try:
                if await member_is_coordinator(
                    channel_obj.id, interaction.guild.id, member.id
                ):
                    coordinators.append(member)
            except commands.CheckFailure:
                pass
            try:
                if await member_is_moderator(
                    channel_obj.id, interaction.guild.id, member.id
                ):
                    moderators.append(member)
            except commands.CheckFailure:
                pass
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
                title=f"{get_random_emoji()} Survey results for {channel_obj.name}",
                description=f"Total surveyed: {len(channel_obj.members)}",
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
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    logger.warning(str(e).capitalize())
                    return await state.end(
                        warning="\U000026a0\U0000fe0f Embed size is too large. Limit the scope."
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
        else:
            try:
                return await state.end(
                    warning="\U000026a0\U0000fe0f No authorized roles found."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @commands.command(name="survey", help="Get all.")
    async def stage_survey_text_command(
        self,
        ctx: commands.Context,
        *,
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
    ):
        state = StateService(ctx)
        channel_obj = None
        chunk_size = 7
        pages = []
        (
            system_owners,
            developers,
            guild_owners,
            administrators,
            moderators,
            coordinators,
        ) = ([], [], [], [], [], [])
        try:
            channel_obj = await resolve_channel(ctx, channel)
        except Exception as e:
            channel_obj = ctx.channel
            logger.warning(str(e).capitalize())
        for member in channel_obj.members:
            try:
                if await member_is_system_owner(member.id):
                    system_owners.append(member)
            except commands.CheckFailure:
                pass
            try:
                if await member_is_developer(ctx.guild.id, member.id):
                    developers.append(member)
            except commands.CheckFailure:
                pass
            try:
                if await member_is_guild_owner(ctx.guild.id, member.id):
                    guild_owners.append(member)
            except commands.CheckFailure:
                pass
            try:
                if await member_is_administrator(ctx.guild.id, member.id):
                    administrators.append(member)
            except commands.CheckFailure:
                pass
            try:
                if await member_is_coordinator(channel_obj.id, ctx.guild.id, member.id):
                    coordinators.append(member)
            except commands.CheckFailure:
                pass
            try:
                if await member_is_moderator(channel_obj.id, ctx.guild.id, member.id):
                    moderators.append(member)
            except commands.CheckFailure:
                pass
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
                title=f"{get_random_emoji()} Survey results for {channel_obj.name}",
                description=f"Total surveyed: {len(channel_obj.members)}",
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
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    logger.warning(str(e).capitalize())
                    return await state.end(
                        warning="\U000026a0\U0000fe0f Embed size is too large. Limit the scope."
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
        else:
            try:
                return await state.end(
                    warning="\U000026a0\U0000fe0f No authorized roles found."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")


async def setup(bot: DiscordBot):
    cog = EveryoneCommands(bot)
    await bot.add_cog(cog)
