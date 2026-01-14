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
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import resolve_channel
from vyrtuous.service.member_service import resolve_member
from vyrtuous.service.message_service import MessageService
from vyrtuous.database.actions.action import Action
from vyrtuous.database.actions.ban import Ban
from vyrtuous.database.settings.cap import Cap
from vyrtuous.utils.properties.duration import DurationObject
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.database.actions.flag import Flag
from vyrtuous.utils.setup_logging import logger
from vyrtuous.utils.properties.snowflake import *
from vyrtuous.database.rooms.stage import Stage
from vyrtuous.service.state_service import StateService
from vyrtuous.database.rooms.temporary_room import TemporaryRoom
from vyrtuous.database.actions.text_mute import TextMute
from vyrtuous.database.roles.vegan import Vegan
from vyrtuous.database.actions.voice_mute import VoiceMute
from vyrtuous.service.scope_service import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    resolve_objects,
    generate_skipped_channels,
    generate_skipped_guilds,
    generate_skipped_members,
    clean_guild_dictionary,
    flush_page
)
import discord
from vyrtuous.service.check_service import (
    moderator_predicator,
    at_home,
    has_equal_or_higher_role,
    not_bot,
    role_check_without_specifics
)
from vyrtuous.utils.properties.snowflake import (
    AppChannelSnowflake,
    AppMemberSnowflake,
    ChannelSnowflake,
    MemberSnowflake
)

class ModeratorCommands(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.message_service = MessageService(self.bot, self.bot.db_pool)

    @app_commands.command(name="bans", description="List bans.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, server ID or empty."
    )
    @moderator_predicator()
    async def list_bans_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        bans, title = await resolve_objects(ctx_interaction_or_message=interaction, obj=Ban, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass
        
        channel_lines, chunk_size, field_count, pages = [], 7, 0, []
        thumbnail = False
        guild_dictionary = {}
    
        for ban in bans:
            guild_dictionary.setdefault(ban.guild_snowflake, {})
            guild_dictionary[ban.guild_snowflake].setdefault(ban.channel_snowflake, [])
            guild_dictionary[ban.guild_snowflake][ban.channel_snowflake].append(
                {
                    "member_snowflake": ban.member_snowflake,
                    "reason": ban.reason,
                    "expires_in": DurationObject.from_expires_in(ban.expires_in),
                }
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_members=skipped_members)

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
                    channel_lines.append(
                        f"**Channel**: {channel.mention}\n**Expires in**: {member_data['expires_in']}\n**Reason**: {member_data['reason']}"
                    )
                    if not thumbnail:
                        embed.set_thumbnail(url=member.display_avatar.url)
                        thumbnail = True
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
                            name=f"Information",
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
            if channel_lines:
                embed.add_field(
                    name=f"Information", value="\n".join(channel_lines), inline=False
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
        await StateService.send_pages(obj=Ban, pages=pages, state=state)
    
    # DONE
    @commands.command(name="bans", description="List bans.")
    @moderator_predicator()
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        bans, title = await resolve_objects(ctx_interaction_or_message=ctx, obj=Ban, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass
        
        channel_lines, chunk_size, field_count, pages = [], 7, 0, []
        thumbnail = False
        guild_dictionary = {}
    
        for ban in bans:
            guild_dictionary.setdefault(ban.guild_snowflake, {})
            guild_dictionary[ban.guild_snowflake].setdefault(ban.channel_snowflake, [])
            guild_dictionary[ban.guild_snowflake][ban.channel_snowflake].append(
                {
                    "member_snowflake": ban.member_snowflake,
                    "reason": ban.reason,
                    "expires_in": DurationObject.from_expires_in(ban.expires_in),
                }
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_members=skipped_members)

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
                    channel_lines.append(
                        f"**Channel**: {channel.mention}\n**Expires in**: {member_data['expires_in']}\n**Reason**: {member_data['reason']}"
                    )
                    if not thumbnail:
                        embed.set_thumbnail(url=member.display_avatar.url)
                        thumbnail = True
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
                            name=f"Information",
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
            if channel_lines:
                embed.add_field(
                    name=f"Information", value="\n".join(channel_lines), inline=False
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
        await StateService.send_pages(obj=Ban, pages=pages, state=state)
    
    # DONE
    @app_commands.command(name="caps", description="List caps.")
    @moderator_predicator()
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, server ID or empty."
    )
    async def list_caps_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        caps, title = await resolve_objects(ctx_interaction_or_message=interaction, obj=Cap, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass
        
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
    
        for cap in caps:
            guild_dictionary.setdefault(cap.guild_snowflake, {})
            guild_dictionary[cap.guild_snowflake].setdefault(cap.channel_snowflake, [])
            channel_entries = guild_dictionary[cap.guild_snowflake][cap.channel_snowflake]
            entry_found = False
            for entry in channel_entries:
                if entry["moderation_type"] == cap.moderation_type:
                    entry["durations"].append(cap.duration)
                    entry_found = True
                    break
            if not entry_found:
                channel_entries.append(
                    {"moderation_type": cap.moderation_type, "durations": [cap.duration]}
                )

        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds)
    
        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            channel_header_added = False
            for channel_snowflake, channel_entries in channels.items():
                lines = []
                channel = guild.get_channel(channel_snowflake)
                for entry in channel_entries:
                    for duration in entry["durations"]:
                        lines.append(
                            f'  ↳ {entry["moderation_type"]} ({DurationObject.from_seconds(duration)})'
                        )
                if not channel_header_added:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    channel_header_added = True
                if field_count + 1 >= chunk_size:
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                embed.add_field(
                    name=f"Channel: {channel.mention}", value="\n".join(lines), inline=False
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
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )

        await StateService.send_pages(obj=Cap, pages=pages, state=state)
    
    # DONE
    @commands.command(name="caps", help="List caps.")
    @moderator_predicator()
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        caps, title = await resolve_objects(ctx_interaction_or_message=ctx, obj=Cap, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass
        
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
    
        for cap in caps:
            guild_dictionary.setdefault(cap.guild_snowflake, {})
            guild_dictionary[cap.guild_snowflake].setdefault(cap.channel_snowflake, [])
            channel_entries = guild_dictionary[cap.guild_snowflake][cap.channel_snowflake]
            entry_found = False
            for entry in channel_entries:
                if entry["moderation_type"] == cap.moderation_type:
                    entry["durations"].append(cap.duration)
                    entry_found = True
                    break
            if not entry_found:
                channel_entries.append(
                    {"moderation_type": cap.moderation_type, "durations": [cap.duration]}
                )
                
        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds)
    
        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            channel_header_added = False
            for channel_snowflake, channel_entries in channels.items():
                lines = []
                channel = guild.get_channel(channel_snowflake)
                for entry in channel_entries:
                    for duration in entry["durations"]:
                        lines.append(
                            f'  ↳ {entry["moderation_type"]} ({DurationObject.from_seconds(duration)})'
                        )
                if not channel_header_added:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    channel_header_added = True
                if field_count + 1 >= chunk_size:
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                embed.add_field(
                    name=f"Channel: {channel.mention}", value="\n".join(lines), inline=False
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
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )

        await StateService.send_pages(obj=Cap, pages=pages, state=state)
    
    # DONE
    @app_commands.command(name="cmds", description="List aliases.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, server ID or empty."
    )
    @moderator_predicator()
    async def list_commands_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        aliases, title = await resolve_objects(ctx_interaction_or_message=interaction, obj=Action, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass
        
        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}
    
        for alias in aliases:
            guild_dictionary.setdefault(alias.guild_snowflake, {})
            guild_dictionary[alias.guild_snowflake].setdefault(alias.channel_snowflake, [])
            channel_entries = guild_dictionary[alias.guild_snowflake][
                alias.channel_snowflake
            ]
            entry_found = False
            for entry in channel_entries:
                if alias.alias_type in entry:
                    entry[alias.alias_type].append(alias.alias_name)
                    entry_found = True
                    break
            if not entry_found:
                channel_entries.append({alias.alias_type: [alias.alias_name]})

        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds)

        for guild_snowflake, channels in guild_dictionary.items():
            current_channel = None
            lines = []
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_entries in channels.items():
                channel = guild.get_channel(channel_snowflake)
                channel_lines = []
                for entry in channel_entries:
                    for alias_type, alias_names in entry.items():
                        channel_lines.append(f"{alias_type}")
                        for name in alias_names:
                            channel_lines.append(f"  ↳ {name}")
                i = 0
                while i < len(channel_lines):
                    remaining_space = chunk_size - len(lines)
                    chunk = channel_lines[i : i + remaining_space]
                    if not lines:
                        current_channel = channel
                    lines.extend(chunk)
                    i += remaining_space
                    if len(lines) >= chunk_size:
                        embed.add_field(
                            name=f"Channel: {current_channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed, field_count = flush_page(embed, pages, title, guild.name)
                        lines = []
                        current_channel = None
            if lines:
                embed.add_field(
                    name=f"Channel: {current_channel.mention}",
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

        await StateService.send_pages(obj=Action, pages=pages, state=state)
    
    # DONE
    @commands.command(name="cmds", help="List aliases.")
    @moderator_predicator()
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        aliases, title = await resolve_objects(ctx_interaction_or_message=ctx, obj=Action, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass
        
        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}
    
        for alias in aliases:
            guild_dictionary.setdefault(alias.guild_snowflake, {})
            guild_dictionary[alias.guild_snowflake].setdefault(alias.channel_snowflake, [])
            channel_entries = guild_dictionary[alias.guild_snowflake][
                alias.channel_snowflake
            ]
            entry_found = False
            for entry in channel_entries:
                if alias.alias_type in entry:
                    entry[alias.alias_type].append(alias.alias_name)
                    entry_found = True
                    break
            if not entry_found:
                channel_entries.append({alias.alias_type: [alias.alias_name]})


        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds)

        for guild_snowflake, channels in guild_dictionary.items():
            current_channel = None
            lines = []
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_entries in channels.items():
                channel = guild.get_channel(channel_snowflake)
                channel_lines = []
                for entry in channel_entries:
                    for alias_type, alias_names in entry.items():
                        channel_lines.append(f"{alias_type}")
                        for name in alias_names:
                            channel_lines.append(f"  ↳ {name}")
                i = 0
                while i < len(channel_lines):
                    remaining_space = chunk_size - len(lines)
                    chunk = channel_lines[i : i + remaining_space]
                    if not lines:
                        current_channel = channel
                    lines.extend(chunk)
                    i += remaining_space
                    if len(lines) >= chunk_size:
                        embed.add_field(
                            name=f"Channel: {current_channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed, field_count = flush_page(embed, pages, title, guild.name)
                        lines = []
                        current_channel = None
            if lines:
                embed.add_field(
                    name=f"Channel: {current_channel.mention}",
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

        await StateService.send_pages(obj=Action, pages=pages, state=state)

    # DONE
    @app_commands.command(name="del", description="Delete message.")
    @app_commands.describe(
        message="Message ID", channel="Tag a channel or include its ID"
    )
    @moderator_predicator()
    async def delete_message_app_command(
        self,
        interaction: discord.Interaction,
        message: AppMessageSnowflake,
        channel: AppChannelSnowflake = None,
    ):
        state = StateService(interaction)
        channel_obj = None
        try:
            channel_obj = await resolve_channel(
                interaction, channel
            )
        except:
            channel_obj = interaction.channel
        msg = await channel_obj.fetch_message(message)
        if not msg:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f Message `{message}` does not exist."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            await has_equal_or_higher_role(
                interaction,
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=msg.author.id,
                sender_snowflake=interaction.user.id,
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                success=f"{get_random_emoji()} Message `{message}` deleted successfully."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @commands.command(name="del", help="Delete message.")
    @moderator_predicator()
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        message: MessageSnowflake = commands.parameter(
            default=None, description="Message snowflake"
        ),
        *,
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Channel or snowflake"
        ),
    ):
        state = StateService(ctx)
        channel_obj = None
        try:
            channel_obj = await resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
        msg = await channel_obj.fetch_message(message)
        if not msg:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f Message `{message}` does not exist."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            await has_equal_or_higher_role(
                ctx,
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=msg.author.id,
                sender_snowflake=ctx.author.id,
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            await msg.delete()
        except discord.Forbidden as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                success=f"{get_random_emoji()} Message `{message}` deleted successfully."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @app_commands.command(name="flags", description="List flags.")
    @moderator_predicator()
    async def list_flags_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        flags, title = await resolve_objects(ctx_interaction_or_message=interaction, obj=Flag, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(ctx_interaction_or_message=interaction, member_str=target)
        except Exception as e:
            pass

        channel_lines, chunk_size, field_count, pages = [], 7, 0, []
        thumbnail = False
        guild_dictionary = {}

        for flag in flags:
            guild_dictionary.setdefault(flag.guild_snowflake, {})
            guild_dictionary[flag.guild_snowflake].setdefault(
                flag.channel_snowflake, []
            )
            guild_dictionary[flag.guild_snowflake][flag.channel_snowflake].append(
                {"member_snowflake": flag.member_snowflake, "reason": flag.reason}
            )

        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_members=skipped_members)

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
                        channel_lines.append(channel.mention)
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
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
                            name=f"Channels",
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
                elif channel_lines:
                    channel_lines.append(channel.mention)
            if channel_lines:
                embed.add_field(
                    name=f"Channels", value="\n".join(channel_lines), inline=False
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

        await StateService.send_pages(obj=Flag, pages=pages, state=state)

    # DONE
    @commands.command(name="flags", help="List flags.")
    @moderator_predicator()
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        flags, title = await resolve_objects(ctx_interaction_or_message=ctx, obj=Flag, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(ctx_interaction_or_message=ctx, member_str=target)
        except Exception as e:
            pass
                
        channel_lines, chunk_size, field_count, pages = [], 7, 0, []
        thumbnail = False
        guild_dictionary = {}

        for flag in flags:
            guild_dictionary.setdefault(flag.guild_snowflake, {})
            guild_dictionary[flag.guild_snowflake].setdefault(
                flag.channel_snowflake, []
            )
            guild_dictionary[flag.guild_snowflake][flag.channel_snowflake].append(
                {"member_snowflake": flag.member_snowflake, "reason": flag.reason}
            )

        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_members=skipped_members)

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
                        channel_lines.append(channel.mention)
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
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
                            name=f"Channels",
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
            if channel_lines:
                embed.add_field(
                    name=f"Channels", value="\n".join(channel_lines), inline=False
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

        await StateService.send_pages(obj=Flag, pages=pages, state=state)

    # DONE
    @app_commands.command(name="ls", description="List new vegans.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, member ID/mention, server ID or empty."
    )
    async def list_new_vegans_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        vegans, title = await resolve_objects(ctx_interaction_or_message=interaction, obj=Vegan, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(ctx_interaction_or_message=interaction, member_str=target)
        except Exception as e:
            pass
                
        channel_lines, chunk_size, field_count, pages = [], 7, 0, []
        thumbnail = False
        guild_dictionary = {}

        for vegan in vegans:
            guild_dictionary.setdefault(vegan.guild_snowflake, {})
            guild_dictionary[vegan.guild_snowflake].setdefault(
                vegan.channel_snowflake, []
            )
            guild_dictionary[vegan.guild_snowflake][vegan.channel_snowflake].append(
                {"member_snowflake": vegan.member_snowflake}
            )

        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_members=skipped_members)

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
                        channel_lines.append(channel.mention)
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
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
                            name=f"Channels",
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
            if channel_lines:
                embed.add_field(
                    name=f"Channels", value="\n".join(channel_lines), inline=False
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

        await StateService.send_pages(obj=Vegan, pages=pages, state=state)

    # DONE
    @commands.command(name="ls", help="List new vegans.")
    async def list_new_vegans_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        vegans, title = await resolve_objects(ctx_interaction_or_message=ctx, obj=Vegan, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(ctx_interaction_or_message=ctx, member_str=target)
        except Exception as e:
            pass
                
        channel_lines, chunk_size, field_count, pages = [], 7, 0, []
        thumbnail = False
        guild_dictionary = {}

        for vegan in vegans:
            guild_dictionary.setdefault(vegan.guild_snowflake, {})
            guild_dictionary[vegan.guild_snowflake].setdefault(
                vegan.channel_snowflake, []
            )
            guild_dictionary[vegan.guild_snowflake][vegan.channel_snowflake].append(
                {"member_snowflake": vegan.member_snowflake}
            )

        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_members=skipped_members)

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
                        channel_lines.append(channel.mention)
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
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
                            name=f"Channels",
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
            if channel_lines:
                embed.add_field(
                    name=f"Channels", value="\n".join(channel_lines), inline=False
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
        state = StateService(interaction)
        channel_obj = None
        old_room = await TemporaryRoom.select_and_room_name(
            guild_snowflake=interaction.guild.id, room_name=old_name
        )
        if old_room:
            try:
                channel_obj = await resolve_channel(
                    interaction, channel
                )
            except Exception as e:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
            is_owner = old_room.member_snowflake == interaction.user.id
            highest_role = await role_check_without_specifics(interaction)
            if (
                highest_role
                not in ("System Owner", "Developer", "Guild Owner", "Administrator")
                and not is_owner
            ):
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f Only owners, developers and administrators can migrate rooms."
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
            await TemporaryRoom.update_channel(
                guild_snowflake=interaction.guild.id,
                room_name=channel_obj.name,
                source_channel_snowflake=old_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            new_room = await TemporaryRoom.select_and_room_name(
                guild_snowflake=interaction.guild.id, room_name=channel_obj.name
            )
            await Action.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Ban.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Cap.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Coordinator.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Flag.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Moderator.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Stage.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await TextMute.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Vegan.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await VoiceMute.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            try:
                return await state.end(
                    success=f"{get_random_emoji()} Temporary room `{old_name}` migrated to {channel_obj.mention}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                warning=f"\U000026a0\U0000fe0f No temporary rooms found called `{old_name}` in {interaction.guild.name}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
        old_name: str = commands.parameter(
            default=None, description="Provide a channel name"
        ),
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include its ID"
        ),
    ):
        state = StateService(ctx)
        channel_obj = None
        old_room = await TemporaryRoom.select_and_room_name(
            guild_snowflake=ctx.guild.id, room_name=old_name
        )
        if old_room:
            try:
                channel_obj = await resolve_channel(ctx, channel)
            except Exception as e:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
            is_owner = old_room.member_snowflake == ctx.author.id
            highest_role = await role_check_without_specifics(ctx)
            if (
                highest_role
                not in ("System Owner", "Developer", "Guild Owner", "Administrator")
                and not is_owner
            ):
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f Only owners, developers and administrators can migrate rooms."
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
            await TemporaryRoom.update_channel(
                guild_snowflake=ctx.guild.id,
                room_name=channel_obj.name,
                source_channel_snowflake=old_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            new_room = await TemporaryRoom.select_and_room_name(
                guild_snowflake=ctx.guild.id, room_name=channel_obj.name
            )
            await Action.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Ban.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Cap.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Coordinator.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Flag.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Moderator.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Stage.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await TextMute.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await Vegan.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            await VoiceMute.update_channel(
                source_channel_snowflake=new_room.channel_snowflake,
                target_channel_snowflake=channel_obj.id,
            )
            try:
                return await state.end(
                    success=f"{get_random_emoji()} Temporary room `{old_name}` migrated to {channel_obj.mention}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                warning=f"\U000026a0\U0000fe0f No temporary rooms found called `{old_name}` in {ctx.guild.name}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @app_commands.command(name="mutes", description="List mutes.")
    @moderator_predicator()
    async def list_mutes_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        voice_mutes, title = await resolve_objects(ctx_interaction_or_message=interaction, obj=VoiceMute, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(ctx_interaction_or_message=interaction, member_str=target)
        except Exception as e:
            pass
                
        channel_lines, chunk_size, field_count, pages = [], 7, 0, []
        thumbnail = False
        guild_dictionary = {}

        for voice_mute in voice_mutes:
            guild_dictionary.setdefault(voice_mute.guild_snowflake, {})
            guild_dictionary[voice_mute.guild_snowflake].setdefault(
                voice_mute.channel_snowflake, []
            )
            guild_dictionary[voice_mute.guild_snowflake][
                voice_mute.channel_snowflake
            ].append(
                {
                    "member_snowflake": voice_mute.member_snowflake,
                    "reason": voice_mute.reason,
                    "expires_in": DurationObject.from_expires_in(voice_mute.expires_in),
                }
            )

        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_members=skipped_members)

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
                        lines.append(
                            f"**User**: {member.mention}\n**Expires in**: {member_data['expires_in']}"
                        )
                    else:
                        channel_lines.append(
                            f"**Channel**: {channel.mention}\n**Expires in**: {member_data['expires_in']}\n**Reason**: {member_data['reason']}"
                        )
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
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
                            name=f"Information",
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
            if channel_lines:
                embed.add_field(
                    name=f"Information", value="\n".join(channel_lines), inline=False
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

        await StateService.send_pages(obj=VoiceMute, pages=pages, state=state)

    # DONE
    @commands.command(name="mutes", help="List mutes.")
    @moderator_predicator()
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, member ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        voice_mutes, title = await resolve_objects(ctx_interaction_or_message=ctx, obj=VoiceMute, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(ctx_interaction_or_message=ctx, member_str=target)
        except Exception as e:
            pass
                
        channel_lines, chunk_size, field_count, pages = [], 7, 0, []
        thumbnail = False
        guild_dictionary = {}

        for voice_mute in voice_mutes:
            guild_dictionary.setdefault(voice_mute.guild_snowflake, {})
            guild_dictionary[voice_mute.guild_snowflake].setdefault(
                voice_mute.channel_snowflake, []
            )
            guild_dictionary[voice_mute.guild_snowflake][
                voice_mute.channel_snowflake
            ].append(
                {
                    "member_snowflake": voice_mute.member_snowflake,
                    "reason": voice_mute.reason,
                    "expires_in": DurationObject.from_expires_in(voice_mute.expires_in),
                }
            )

        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_members=skipped_members)

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
                        lines.append(
                            f"**User**: {member.mention}\n**Expires in**: {member_data['expires_in']}"
                        )
                    else:
                        channel_lines.append(
                            f"**Channel**: {channel.mention}\n**Expires in**: {member_data['expires_in']}\n**Reason**: {member_data['reason']}"
                        )
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
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
                            name=f"Information",
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
            if channel_lines:
                embed.add_field(
                    name=f"Information", value="\n".join(channel_lines), inline=False
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
        state = StateService(interaction)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await resolve_channel(
                interaction, channel
            )
        except:
            channel_obj = interaction.channel
        try:
            member_obj = await resolve_member(interaction, member)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f Could not resolve a valid member `{member}`."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            not_bot(interaction, member_snowflake=member_obj.id)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f You are not authorized to affect {interaction.guild.me.mention}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            highest_role = await has_equal_or_higher_role(
                interaction,
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id,
                sender_snowflake=interaction.user.id,
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        stage = await Stage.select(
            channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id
        )
        if not stage:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f No active stage found in {channel_obj.mention}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            await member_obj.edit(mute=not member_obj.voice.mute)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                success=f"{get_random_emoji()} {member_obj.mention} has been {'muted' if member_obj.voice.mute else 'unmuted'}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @commands.command(name="mstage", help="Stage mute/unmute.")
    @moderator_predicator()
    async def stage_mute_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            default=None, description="Tag a member or include their ID"
        ),
        channel: ChannelSnowflake = commands.parameter(
            default=None, description="Tag a channel or include it's ID"
        ),
    ):
        state = StateService(ctx)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await resolve_channel(ctx, channel)
        except Exception as e:
            channel_obj = ctx.channel
        try:
            member_obj = await resolve_member(ctx, member)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f Could not resolve a valid member `{member}`."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            not_bot(ctx, member_snowflake=member_obj.id)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f You are not authorized to affect {ctx.guild.me.mention}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            highest_role = await has_equal_or_higher_role(
                ctx,
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id,
                sender_snowflake=ctx.author.id,
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        stage = await Stage.select(
            channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id
        )
        if not stage:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f No active stage found in {channel_obj.mention}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            await member_obj.edit(mute=not member_obj.voice.mute)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                success=f"{get_random_emoji()} {member_obj.mention} has been {'muted' if member_obj.voice.mute else 'unmuted'}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @app_commands.command(name="stages", description="List stages.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, server ID or empty."
    )
    @moderator_predicator()
    async def list_stages_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        stages, title = await resolve_objects(ctx_interaction_or_message=interaction, obj=Stage, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass
                
        chunk_size, field_count, pages = 7, 0, []
        guild_dictionary = {}

        for stage in stages:
            guild_dictionary.setdefault(stage.guild_snowflake, {})
            guild_dictionary[stage.guild_snowflake].setdefault(
                stage.channel_snowflake, []
            )
            guild_dictionary[stage.guild_snowflake][stage.channel_snowflake][
                "expires_in"
            ] = DurationObject.from_expires_in(stage.expires_in)

        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds)

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                lines = []
                lines.append(f"**Expires in**: {channel_data['expires_in']}")
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
            default=None,
            description="Specify one of: 'all', channel ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        stages, title = await resolve_objects(ctx_interaction_or_message=ctx, obj=Stage, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass
                
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}

        for stage in stages:
            guild_dictionary.setdefault(stage.guild_snowflake, {})
            guild_dictionary[stage.guild_snowflake].setdefault(
                stage.channel_snowflake, []
            )
            guild_dictionary[stage.guild_snowflake][stage.channel_snowflake][
                "expires_in"
            ] = DurationObject.from_expires_in(stage.expires_in)

        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds)

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                lines = []
                lines.append(f"**Expires in**: {channel_data['expires_in']}")
                field_count += 1
                if field_count == chunk_size:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
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
        state = StateService(interaction)

        chunk_size, field_count, lines, pages = 7, 0, [], []

        try:
            member_obj = await resolve_member(interaction, member)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
            
        bans, ban_title = await resolve_objects(ctx_interaction_or_message=interaction, obj=Ban, state=state, target=target)
        flags, flag_title = await resolve_objects(ctx_interaction_or_message=interaction, obj=Flag, state=state, target=target)
        text_mutes, text_title = await resolve_objects(ctx_interaction_or_message=interaction, obj=TextMute, state=state, target=target)
        vegans, vegan_title = await resolve_objects(ctx_interaction_or_message=interaction, obj=Vegan, state=state, target=target)
        voice_mutes, voice_mute_title = await resolve_objects(ctx_interaction_or_message=interaction, obj=Ban, state=state, target=target)
        
        guild = self.bot.get_guild(interaction.guild.id)

        embed = discord.Embed(
            title=ban_title, description=guild.name, color=discord.Color.blue()
        )
        if bans:
            for ban in bans:
                if ban:
                    channel = guild.get_channel(ban.channel_snowflake)
                    lines.append(
                        f"**User**: {member_obj.mention}\n**Expires in**: {ban.expires_in}\n**Reason**: {ban.reason}"
                    )
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed.set_thumbnail(url=member_obj.display_avatar.url)
                        embed, field_count = flush_page(embed, pages, ban_title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed.set_thumbnail(url=member_obj.display_avatar.url)
                pages.append(embed)
                lines = []

        embed = discord.Embed(
            title=flag_title, description=guild.name, color=discord.Color.blue()
        )
        if flags:
            for flag in flags:
                if flag:
                    channel = guild.get_channel(flag.channel_snowflake)
                    lines.append(
                        f"**User**: {member_obj.mention}\n**Reason**: {flag.reason}"
                    )
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed.set_thumbnail(url=member_obj.display_avatar.url)
                        embed, field_count = flush_page(embed, pages, flag_title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed.set_thumbnail(url=member_obj.display_avatar.url)
                pages.append(embed)
                lines = []

        embed = discord.Embed(
            title=text_title, description=guild.name, color=discord.Color.blue()
        )
        if text_mutes:
            for text_mute in text_mutes:
                if text_mute:
                    channel = guild.get_channel(text_mute.channel_snowflake)
                    lines.append(
                        f"**User**: {member_obj.mention}\n**Expires in**: {text_mute.expires_in}\n**Reason**: {text_mute.reason}"
                    )
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed.set_thumbnail(url=member_obj.display_avatar.url)
                        embed, field_count = flush_page(embed, pages, ban_title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed.set_thumbnail(url=member_obj.display_avatar.url)
                pages.append(embed)
                lines = []

        embed = discord.Embed(
            title=vegan_title, description=guild.name, color=discord.Color.blue()
        )
        if vegans:
            for vegan in vegans:
                if vegan:
                    channel = guild.get_channel(vegan.channel_snowflake)
                    lines.append(f"**User**: {member_obj.mention}")
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed.set_thumbnail(url=member_obj.display_avatar.url)
                        embed, field_count = flush_page(embed, pages, vegan_title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed.set_thumbnail(url=member_obj.display_avatar.url)
                pages.append(embed)
                lines = []

        embed = discord.Embed(
            title=voice_mute_title, description=guild.name, color=discord.Color.blue()
        )
        if voice_mutes:
            for voice_mute in voice_mutes:
                if voice_mute:
                    channel = guild.get_channel(voice_mute.channel_snowflake)
                    lines.append(
                        f"**User**: {member_obj.mention}\n**Expires in**: {voice_mute.expires_in}\n**Reason**: {voice_mute.reason}"
                    )
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed.set_thumbnail(url=member_obj.display_avatar.url)
                        embed, field_count = flush_page(embed, pages, voice_mute_title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                pages.append(embed)
        
        await StateService.send_pages(obj=Action, pages=pages, state=state)

    @commands.command(name="summary", description="Moderation summary.")
    @moderator_predicator()
    async def list_moderation_summary_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            default=None, description="Specify a member ID/mention."
        ),
        *,
        target: Optional[str] = commands.parameter(
            default=None, description="Specify 'all' or a channel ID/mention."
        ),
    ):
        state = StateService(ctx)

        chunk_size, field_count, lines, pages = 7, 0, [], []

        try:
            member_obj = await resolve_member(ctx, member)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
            
        bans, ban_title = await resolve_objects(ctx_interaction_or_message=ctx, obj=Ban, state=state, target=target)
        flags, flag_title = await resolve_objects(ctx_interaction_or_message=ctx, obj=Flag, state=state, target=target)
        text_mutes, text_title = await resolve_objects(ctx_interaction_or_message=ctx, obj=TextMute, state=state, target=target)
        vegans, vegan_title = await resolve_objects(ctx_interaction_or_message=ctx, obj=Vegan, state=state, target=target)
        voice_mutes, voice_mute_title = await resolve_objects(ctx_interaction_or_message=ctx, obj=Ban, state=state, target=target)
        
        guild = self.bot.get_guild(ctx.guild.id)
        
        embed = discord.Embed(
            title=ban_title, description=guild.name, color=discord.Color.blue()
        )
        if bans:
            for ban in bans:
                if ban:
                    channel = guild.get_channel(ban.channel_snowflake)
                    lines.append(
                        f"**User**: {member_obj.mention}\n**Expires in**: {ban.expires_in}\n**Reason**: {ban.reason}"
                    )
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed.set_thumbnail(url=member_obj.display_avatar.url)
                        embed, field_count = flush_page(embed, pages, ban_title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed.set_thumbnail(url=member_obj.display_avatar.url)
                pages.append(embed)
                lines = []

        embed = discord.Embed(
            title=flag_title, description=guild.name, color=discord.Color.blue()
        )
        if flags:
            for flag in flags:
                if flag:
                    channel = guild.get_channel(flag.channel_snowflake)
                    lines.append(
                        f"**User**: {member_obj.mention}\n**Reason**: {flag.reason}"
                    )
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed.set_thumbnail(url=member_obj.display_avatar.url)
                        embed, field_count = flush_page(embed, pages, flag_title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed.set_thumbnail(url=member_obj.display_avatar.url)
                pages.append(embed)
                lines = []

        embed = discord.Embed(
            title=text_title, description=guild.name, color=discord.Color.blue()
        )
        if text_mutes:
            for text_mute in text_mutes:
                if text_mute:
                    channel = guild.get_channel(text_mute.channel_snowflake)
                    lines.append(
                        f"**User**: {member_obj.mention}\n**Expires in**: {text_mute.expires_in}\n**Reason**: {text_mute.reason}"
                    )
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed.set_thumbnail(url=member_obj.display_avatar.url)
                        embed, field_count = flush_page(embed, pages, ban_title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed.set_thumbnail(url=member_obj.display_avatar.url)
                pages.append(embed)
                lines = []

        embed = discord.Embed(
            title=vegan_title, description=guild.name, color=discord.Color.blue()
        )
        if vegans:
            for vegan in vegans:
                if vegan:
                    channel = guild.get_channel(vegan.channel_snowflake)
                    lines.append(f"**User**: {member_obj.mention}")
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed.set_thumbnail(url=member_obj.display_avatar.url)
                        embed, field_count = flush_page(embed, pages, vegan_title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed.set_thumbnail(url=member_obj.display_avatar.url)
                pages.append(embed)
                lines = []

        embed = discord.Embed(
            title=voice_mute_title, description=guild.name, color=discord.Color.blue()
        )
        if voice_mutes:
            for voice_mute in voice_mutes:
                if voice_mute:
                    channel = guild.get_channel(voice_mute.channel_snowflake)
                    lines.append(
                        f"**User**: {member_obj.mention}\n**Expires in**: {voice_mute.expires_in}\n**Reason**: {voice_mute.reason}"
                    )
                    field_count += 1
                    if field_count == chunk_size:
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed.set_thumbnail(url=member_obj.display_avatar.url)
                        embed, field_count = flush_page(embed, pages, voice_mute_title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                pages.append(embed)
        
        await StateService.send_pages(obj=Action, pages=pages, state=state)

    # DONE
    @app_commands.command(name="tmutes", description="List text-mutes.")
    @app_commands.describe(
        target="Specify one of: 'all', channel ID/mention, member ID/mention, server ID or empty."
    )
    @moderator_predicator()
    async def list_text_mutes_app_command(
        self, interaction: discord.Interaction, target: Optional[str] = None
    ):
        state = StateService(interaction)
        text_mutes, title = await resolve_objects(ctx_interaction_or_message=interaction, obj=TextMute, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(ctx_interaction_or_message=interaction, member_str=target)
        except Exception as e:
            pass
                
        channel_lines, chunk_size, field_count, lines, pages = [], 7, 0, [], []
        thumbnail = False
        guild_dictionary = {}

        for text_mute in text_mutes:
            guild_dictionary.setdefault(text_mute.guild_snowflake, {})
            guild_dictionary[text_mute.guild_snowflake].setdefault(
                text_mute.channel_snowflake, []
            )
            guild_dictionary[text_mute.guild_snowflake][
                text_mute.channel_snowflake
            ].append(
                {
                    "member_snowflake": text_mute.member_snowflake,
                    "reason": text_mute.reason,
                    "expires_in": DurationObject.from_expires_in(text_mute.expires_in),
                }
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_members=skipped_members)

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
                    if not member:
                        skipped_members.setdefault(
                            guild_snowflake, []
                        ).append(member_data["member_snowflake"])
                        continue
                    if not member_obj:
                        lines.append(
                            f"**User**: {member.mention}\n**Expires in**: {member_data['expires_in']}"
                        )
                    else:
                        channel_lines.append(
                            f"**Channel**: {channel.mention}\n**Expires in**: {member_data['expires_in']}\n**Reason**: {member_data['reason']}"
                        )
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
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
                            name=f"Information",
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
            if channel_lines:
                embed.add_field(
                    name=f"Information", value="\n".join(channel_lines), inline=False
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

        await StateService.send_pages(obj=TextMute, pages=pages, state=state)

    # DONE
    @commands.command(name="tmutes", help="List text-mutes.")
    @moderator_predicator()
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: 'all', channel ID/mention, server ID or empty.",
        ),
    ):
        state = StateService(ctx)
        text_mutes, title = await resolve_objects(ctx_interaction_or_message=ctx, obj=TextMute, state=state, target=target)
    
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass

        member_obj = None
        try:
            member_obj = await resolve_member(ctx_interaction_or_message=ctx, member_str=target)
        except Exception as e:
            pass
                
        channel_lines, chunk_size, field_count, lines, pages = [], 7, 0, [], []
        thumbnail = False
        guild_dictionary = {}

        for text_mute in text_mutes:
            guild_dictionary.setdefault(text_mute.guild_snowflake, {})
            guild_dictionary[text_mute.guild_snowflake].setdefault(
                text_mute.channel_snowflake, []
            )
            guild_dictionary[text_mute.guild_snowflake][
                text_mute.channel_snowflake
            ].append(
                {
                    "member_snowflake": text_mute.member_snowflake,
                    "reason": text_mute.reason,
                    "expires_in": DurationObject.from_expires_in(text_mute.expires_in),
                }
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_channels = await generate_skipped_channels(guild_dictionary)
        skipped_members = await generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_members=skipped_members)

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
                        lines.append(
                            f"**User**: {member.mention}\n**Expires in**: {member_data['expires_in']}"
                        )
                    else:
                        channel_lines.append(
                            f"**Channel**: {channel.mention}\n**Expires in**: {member_data['expires_in']}\n**Reason**: {member_data['reason']}"
                        )
                        if not thumbnail:
                            embed.set_thumbnail(url=member.display_avatar.url)
                            thumbnail = True
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
                            name=f"Information",
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
            if channel_lines:
                embed.add_field(
                    name=f"Information", value="\n".join(channel_lines), inline=False
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
                
        await StateService.send_pages(obj=TextMute, pages=pages, state=state)

async def setup(bot: DiscordBot):
    cog = ModeratorCommands(bot)
    await bot.add_cog(cog)
