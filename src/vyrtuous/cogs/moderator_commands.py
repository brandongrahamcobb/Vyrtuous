''' moderator_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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
'''
from collections import defaultdict
from typing import Optional
from discord import app_commands
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.message_service import MessageService
from vyrtuous.utils.paginator import Paginator
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.all import All
from vyrtuous.utils.ban import Ban
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.duration import Duration, DurationObject
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.flag import Flag
from vyrtuous.utils.setup_logging import logger
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.stage import Stage
from vyrtuous.utils.state import State
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.vegan import Vegan
from vyrtuous.utils.voice_mute import VoiceMute
   
class ModeratorCommands(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.emoji = Emojis()
        self.message_service = MessageService(self.bot, self.bot.db_pool)
        self.channel_service = ChannelService()
        self.member_service = MemberService()
    
    # DONE
    @app_commands.command(name='bans', description='Lists ban statistics.')
    @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, server ID or empty.")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_bans_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        bans, pages = [], []
        channel_obj = None
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        skipped_channel_snowflakes_by_guild_snowflake, skipped_member_snowflakes_by_guild_snowflake = {}, {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Ban(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list bans across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            bans = await Ban.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope) 
                bans = await Ban.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(interaction, scope)
                    bans = await Ban.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
                except Exception as e:
                    if highest_role not in ('Owner', 'Developer', 'Administrator'):
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list bans for specific servers.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    bans = await Ban.fetch_by_guild(guild_snowflake=int(scope))
        else:
            bans = await Ban.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
            channel_obj = interaction.channel 

        if not bans:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No ban(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for ban in bans:
            guild_dictionary.setdefault(ban.guild_snowflake, {})
            guild_dictionary[ban.guild_snowflake].setdefault(ban.channel_snowflake, [])
            guild_dictionary[ban.guild_snowflake][ban.channel_snowflake].append({
                'member_snowflake': ban.member_snowflake,
                'reason': ban.reason,
                'expires_in': DurationObject.from_expires_in(ban.expires_in)
            })
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for member_data in members:
                    member = guild.get_member(member_data['member_snowflake'])
                    if not member:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_data['member_snowflake'])
                        continue
                    entry_lines = [f'**User**: {member.mention}', f'**Expires in**: {member_data["expires_in"]}']
                    if member_obj:
                        entry_lines.append(f'**Reason**: {member_data["reason"]}')
                    if field_count >= chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(entry_lines), inline=False)
                    field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
        if skipped_member_snowflakes_by_guild_snowflake:
            for guild_snowflake, member_list in skipped_member_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for member_snowflake in member_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(member_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No bans found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='bans', description='Lists ban statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
       scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', channel ID/mention, server ID or empty.")
    ):
        state = State(ctx)
        bans, pages = [], []
        channel_obj = None
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        skipped_channel_snowflakes_by_guild_snowflake, skipped_member_snowflakes_by_guild_snowflake = {}, {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Ban(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list bans across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            bans = await Ban.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
                bans = await Ban.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(ctx, scope)
                    bans = await Ban.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                except Exception as e:
                    if highest_role not in ('Owner', 'Developer', 'Administrator'):
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list bans for specific servers.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    bans = await Ban.fetch_by_guild(guild_snowflake=int(scope))
        else:
            bans = await Ban.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            channel_obj = ctx.channel

        if not bans:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No ban(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        guild_dictionary = {}
        for ban in bans:
            guild_dictionary.setdefault(ban.guild_snowflake, {})
            guild_dictionary[ban.guild_snowflake].setdefault(ban.channel_snowflake, [])
            guild_dictionary[ban.guild_snowflake][ban.channel_snowflake].append({
                'member_snowflake': ban.member_snowflake,
                'reason': ban.reason,
                'expires_in': DurationObject.from_expires_in(ban.expires_in)
            })
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for member_data in members:
                    member = guild.get_member(member_data['member_snowflake'])
                    if not member:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_data['member_snowflake'])
                        continue
                    entry_lines = [f'**User**: {member.mention}', f'**Expires in**: {member_data["expires_in"]}']
                    if member_obj:
                        entry_lines.append(f'**Reason**: {member_data["reason"]}')
                    if field_count >= chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(entry_lines), inline=False)
                    field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
        if skipped_member_snowflakes_by_guild_snowflake:
            for guild_snowflake, member_list in skipped_member_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for member_snowflake in member_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(member_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No bans found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @app_commands.command(name='caps', description="List active caps for a channel or all channels if 'all' is provided.")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, server ID or empty.")
    async def list_caps_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        field_count = 0
        caps, lines, pages = [], [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Cap(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list caps across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            caps = await Cap.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
                caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list caps for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, server ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                caps = await Cap.fetch_by_guild(guild_snowflake=int(scope))
        else:
            caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
            channel_obj = interaction.channel

        if not caps:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No caps setup for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for cap in caps:
            guild_dictionary.setdefault(cap.guild_snowflake, {})
            guild_dictionary[cap.guild_snowflake].setdefault(cap.channel_snowflake, {})
            guild_dictionary[cap.guild_snowflake][cap.channel_snowflake].setdefault(cap.moderation_type, [])
            guild_dictionary[cap.guild_snowflake][cap.channel_snowflake][cap.moderation_type].append(cap.duration)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))
        
        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                lines = []
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for moderation_type, durations in channel_data.items():
                    for duration in durations:
                        lines.append(f'  ↳ {moderation_type} ({DurationObject.from_seconds(duration)})')
                if field_count + 1 >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(lines), inline=False)
                field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No caps found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @commands.command(name='caps', help="List active caps for a channel or all channels if 'all' is provided.")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', channel ID/mention, server ID or empty.")
    ):
        state = State(ctx)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        field_count = 0
        caps, lines, pages = [], [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Cap(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list caps across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            caps = await Cap.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
                caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list caps for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, server ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                caps = await Cap.fetch_by_guild(guild_snowflake=int(scope))
        else:
            caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            channel_obj = ctx.channel
        
        if not caps:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No caps setup for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for cap in caps:
            guild_dictionary.setdefault(cap.guild_snowflake, {})
            guild_dictionary[cap.guild_snowflake].setdefault(cap.channel_snowflake, {})
            guild_dictionary[cap.guild_snowflake][cap.channel_snowflake].setdefault(cap.moderation_type, [])
            guild_dictionary[cap.guild_snowflake][cap.channel_snowflake][cap.moderation_type].append(cap.duration)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))
        
        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                lines = []
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for moderation_type, durations in channel_data.items():
                    for duration in durations:
                        lines.append(f'  ↳ {moderation_type} ({DurationObject.from_seconds(duration)})')
                if field_count + 1 >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(lines), inline=False)
                field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No caps found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
 
    # DONE
    @app_commands.command(name='cmds', description="List command aliases routed to a specific channel, temp room, or all channels if 'all' is provided.")
    @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, server ID or empty.")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_commands_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        field_count = 0
        aliases, lines, pages = [], [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Alias(es)'
        
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, scope)
        except:
            channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all aliases in all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            aliases = await Alias.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope) 
                aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list aliases for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, server ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                aliases = await Alias.fetch_by_guild(guild_snowflake=int(scope))
        else:
            aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
            channel_obj = interaction.channel

        if not aliases:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No aliases found for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

        guild_dictionary = {}
        for alias in aliases:
            guild_dictionary.setdefault(alias.guild_snowflake, {})
            guild_dictionary[alias.guild_snowflake].setdefault(alias.channel_snowflake, {})
            guild_dictionary[alias.guild_snowflake][alias.channel_snowflake].setdefault(alias.alias_type, [])
            guild_dictionary[alias.guild_snowflake][alias.channel_snowflake][alias.alias_type].append(alias.alias_name)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            current_channel = None
            lines = []
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                channel_lines = []
                for alias_type, alias_names in channel_data.items():
                    channel_lines.append(f"{alias_type}")
                    for name in alias_names:
                        channel_lines.append(f"  ↳ {name}")
                i = 0
                while i < len(channel_lines):
                    remaining_space = chunk_size - len(lines)
                    chunk = channel_lines[i:i + remaining_space]
                    if not lines:
                        current_channel = channel
                    lines.extend(chunk)
                    i += remaining_space
                    if len(lines) >= chunk_size:
                        embed.add_field(name=f'Channel: {current_channel.mention}', value='\n'.join(lines), inline=False)
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        lines = []
                        current_channel = None
            if lines:
                embed.add_field(name=f'Channel: {current_channel.mention}', value='\n'.join(lines), inline=False)
                pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No aliases found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='cmds', help="List command aliases routed to a specific channel, temp room, or all channels if 'all' is provided.")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', channel ID/mention, server ID or empty.")
    ):
        state = State(ctx)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        field_count = 0
        aliases, lines, pages = [], [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Alias(es)'

        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, scope)
        except:
            channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all aliases in all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            aliases = await Alias.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope) 
                aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list aliases for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, server ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                aliases = await Alias.fetch_by_guild(guild_snowflake=int(scope))
        else:
            aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            channel_obj = ctx.channel

        if not aliases:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No aliases found for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

        guild_dictionary = {}
        for alias in aliases:
            guild_dictionary.setdefault(alias.guild_snowflake, {})
            guild_dictionary[alias.guild_snowflake].setdefault(alias.channel_snowflake, {})
            guild_dictionary[alias.guild_snowflake][alias.channel_snowflake].setdefault(alias.alias_type, [])
            guild_dictionary[alias.guild_snowflake][alias.channel_snowflake][alias.alias_type].append(alias.alias_name)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            current_channel = None
            lines = []
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                channel_lines = []
                for alias_type, alias_names in channel_data.items():
                    channel_lines.append(f"{alias_type}")
                    for name in alias_names:
                        channel_lines.append(f"  ↳ {name}")
                i = 0
                while i < len(channel_lines):
                    remaining_space = chunk_size - len(lines)
                    chunk = channel_lines[i:i + remaining_space]
                    if not lines:
                        current_channel = channel
                    lines.extend(chunk)
                    i += remaining_space
                    if len(lines) >= chunk_size:
                        embed.add_field(name=f'Channel: {current_channel.mention}', value='\n'.join(lines), inline=False)
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        lines = []
                        current_channel = None
            if lines:
                embed.add_field(name=f'Channel: {current_channel.mention}', value='\n'.join(lines), inline=False)
                pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No aliases found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
 
    # DONE
    @app_commands.command(name='del', description='Delete a message by ID (only if you are coordinator/moderator of that temp room).')
    @app_commands.describe(
        message='Message snowflake ID',
        channel='Tag a channel or include its snowflake ID'
    )
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def delete_message_app_command(
        self,
        interaction: discord.Interaction,
        message: AppMessageSnowflake,
        channel: AppChannelSnowflake = None
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
        msg = await channel_obj.fetch_message(message)
        if not msg:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F Message `{message}` does not exist.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            await has_equal_or_higher_role(interaction, channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=msg.author.id, sender_snowflake=interaction.user.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            await msg.delete()
        except discord.Forbidden:
            try:
                return await state.end(warning='\U000026A0\U0000FE0F Missing permissions to delete the message.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Message `{message}` deleted successfully.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @commands.command(name='del', help='Delete a message by ID (only if you are coordinator/moderator of that temp room).')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        message: MessageSnowflake = commands.parameter(default=None, description='Message snowflake'),
        *,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Channel or snowflake')
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
        msg = await channel_obj.fetch_message(message)
        if not msg:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F Message `{message}` does not exist.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            await has_equal_or_higher_role(ctx, channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=msg.author.id, sender_snowflake=ctx.author.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            await msg.delete()
        except discord.Forbidden:
            try:
                return await state.end(warning='\U000026A0\U0000FE0F Missing permissions to delete the message.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Message `{message}` deleted successfully.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @app_commands.command(name='flags', description='List flag statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_flags_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        flages, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake, skipped_member_snowflakes_by_guild_snowflake = {}, {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Flag(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list flags across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            flags = await Flag.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope) 
                title = f'{self.emoji.get_random_emoji()} Flags'
                flags = await Flag.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(interaction, scope) 
                    title = f'{self.emoji.get_random_emoji()} Flags'
                    flags = await Flag.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
                except Exception as e:
                    if highest_role not in ('Owner', 'Developer', 'Administrator'):
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list flags for specific servers.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    title = f'{self.emoji.get_random_emoji()} Flags'
                    flags = await Flag.fetch_by_guild(guild_snowflake=int(scope))
        else:
            title = f'{self.emoji.get_random_emoji()} Flags'
            flags = await Flag.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
        
        if not flags:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No flag(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for flag in flags:
            guild_dictionary.setdefault(flag.guild_snowflake, {})
            guild_dictionary[flag.guild_snowflake].setdefault(flag.channel_snowflake, [])
            guild_dictionary[flag.guild_snowflake][flag.channel_snowflake].append({
                'member_snowflake': flag.member_snowflake,
                'reason': flag.reason
            })
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for member_data in members:
                    member = guild.get_member(member_data['member_snowflake'])
                    if not member:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_data['member_snowflake'])
                        continue
                    entry_lines = [f'**User**: {member.mention}']
                    if member_obj:
                        entry_lines.append(f'**Reason**: {member_data["reason"]}')
                    if field_count >= chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(entry_lines), inline=False)
                    field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No flags found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @commands.command(name='flags', help='List flag statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', channel ID/mention, member ID/mention, server ID or empty.")
    ):
        state = State(ctx)
        channel_obj = None
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        flags, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake, skipped_member_snowflakes_by_guild_snowflake = {}, {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Flag(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list flags across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            flags = await Flag.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope) 
                flags = await Flag.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(ctx, scope)
                    flags = await Flag.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                except Exception as e:
                    if highest_role not in ('Owner', 'Developer', 'Administrator'):
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list flags for specific servers.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    flags = await Flag.fetch_by_guild(guild_snowflake=int(scope))
        else:
            flags = await Flag.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            channel_obj = ctx.channel
        
        if not flags:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No flag(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for flag in flags:
            guild_dictionary.setdefault(flag.guild_snowflake, {})
            guild_dictionary[flag.guild_snowflake].setdefault(flag.channel_snowflake, [])
            guild_dictionary[flag.guild_snowflake][flag.channel_snowflake].append({
                'member_snowflake': flag.member_snowflake,
                'reason': flag.reason
            })
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for member_data in members:
                    member = guild.get_member(member_data['member_snowflake'])
                    if not member:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_data['member_snowflake'])
                        continue
                    entry_lines = [f'**User**: {member.mention}']
                    if member_obj:
                        entry_lines.append(f'**Reason**: {member_data["reason"]}')
                    if field_count >= chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(entry_lines), inline=False)
                    field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
        if skipped_member_snowflakes_by_guild_snowflake:
            for guild_snowflake, member_list in skipped_member_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for member_snowflake in member_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(member_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No flags found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @app_commands.command(name='ls', description='List users veganed as going vegan in this server.')
    @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, member ID/mention, server ID or empty.")
    async def list_vegans_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        field_count = 0
        pages = []
        skipped_channel_snowflakes_by_guild_snowflake, skipped_member_snowflakes_by_guild_snowflake = {}, {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Vegan(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list vegans across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            vegans = await Vegan.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
                vegans = await Vegan.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(interaction, scope)
                    vegans = await Vegan.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
                except Exception as e:
                    if highest_role not in ('Owner', 'Developer', 'Administrator'):
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list vegans for specific servers.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    vegans = await Vegan.fetch_by_guild(guild_snowflake=int(scope))
        else:
            vegans = await Vegan.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
            channel_obj = interaction.channel
        
        if not vegans:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No vegan(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for vegan in vegans:
            guild_dictionary.setdefault(vegan.guild_snowflake, {})
            guild_dictionary[vegan.guild_snowflake].setdefault(vegan.channel_snowflake, [])
            guild_dictionary[vegan.guild_snowflake][vegan.channel_snowflake].append({
                'member_snowflake': vegan.member_snowflake
            })
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for entry in channel_data:
                    member = guild.get_member(entry['member_snowflake'])
                    if not member:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(entry['member_snowflake'])
                        continue
                    entry_lines = [f'**User**: {member.mention}']
                    if field_count >= chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(entry_lines), inline=False)
                    field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
        if skipped_member_snowflakes_by_guild_snowflake:
            for guild_snowflake, member_list in skipped_member_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for member_snowflake in member_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(member_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No vegans found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
                            
    # DONE
    @commands.command(name='ls', help='List users veganed as going vegan in this server.')
    async def list_vegans_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', channel ID/mention, member ID/mention, server ID or empty.")
    ):
        state = State(ctx)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        field_count = 0
        pages = []
        skipped_channel_snowflakes_by_guild_snowflake, skipped_member_snowflakes_by_guild_snowflake = {}, {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Vegan(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list vegans across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            vegans = await Vegan.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope) 
                vegans = await Vegan.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(ctx, scope)
                    vegans = await Vegan.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                except Exception as e:
                    if highest_role not in ('Owner', 'Developer', 'Administrator'):
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list vegans for specific servers.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    vegans = await Vegan.fetch_by_guild(guild_snowflake=int(scope))
        else:
            vegans = await Vegan.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            channel_obj = ctx.channel
        
        if not vegans:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No vegan(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for vegan in vegans:
            guild_dictionary.setdefault(vegan.guild_snowflake, {})
            guild_dictionary[vegan.guild_snowflake].setdefault(vegan.channel_snowflake, [])
            guild_dictionary[vegan.guild_snowflake][vegan.channel_snowflake].append({
                'member_snowflake': vegan.member_snowflake
            })
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            field_count = 0
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for entry in channel_data:
                    member = guild.get_member(entry['member_snowflake'])
                    if not member:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(entry['member_snowflake'])
                        continue
                    entry_lines = [f'**User**: {member.mention}']
                    if field_count >= chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(entry_lines), inline=False)
                    field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
        if skipped_member_snowflakes_by_guild_snowflake:
            for guild_snowflake, member_list in skipped_member_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for member_snowflake in member_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(member_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No vegans found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

 # DONE
    @app_commands.command(name='migrate', description='Migrate a temporary room to a new channel.')
    @app_commands.describe(
        old_name='Old temporary room name',
        channel='New channel to migrate to'
    )
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def migrate_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        old_name: str,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        old_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=interaction.guild.id, room_name=old_name)
        if old_room:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, channel)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            is_owner = old_room.member_snowflake == interaction.user.id
            highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
            if highest_role not in ('Owner', 'Developer', 'Administrator') and not is_owner:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners, developers and administrators can migrate rooms.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            await TemporaryRoom.update_by_source_and_target(guild_snowflake=interaction.guild.id, room_name=channel_obj.name, source_channel_snowflake=old_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            new_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=interaction.guild.id, room_name=channel_obj.name)
            await Alias.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Ban.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Cap.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Coordinator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Flag.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Moderator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Stage.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await TextMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Vegan.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await VoiceMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Temporary room `{old_name}` migrated to {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(warning=f'\U000026A0\U0000FE0F No temporary rooms found called `{old_name}` in {interaction.guild.name}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @commands.command(name='migrate', help='Migrate a temporary room to a new channel by snowflake.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def migrate_temp_room_text_command(
        self,
        ctx: commands.Context,
        old_name: Optional[str] = commands.parameter(default=None, description='Provide a channel name'),
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        channel_obj = None
        old_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=ctx.guild.id, room_name=old_name)
        if old_room:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, channel)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            is_owner = old_room.member_snowflake == ctx.author.id
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
            if highest_role not in ('Owner', 'Developer', 'Administrator') and not is_owner:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners, developers and administrators can migrate rooms.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            await TemporaryRoom.update_by_source_and_target(guild_snowflake=ctx.guild.id, room_name=channel_obj.name, source_channel_snowflake=old_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            new_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=ctx.guild.id, room_name=channel_obj.name)
            await Alias.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Ban.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Cap.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Coordinator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Flag.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Moderator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Stage.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await TextMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Vegan.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await VoiceMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Temporary room `{old_name}` migrated to {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(warning=f'\U000026A0\U0000FE0F No temporary rooms found called `{old_name}` in {ctx.guild.name}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @app_commands.command(name='mutes', description='Lists mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_mutes_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        voice_mutes, pages = [], []
        channel_obj = None
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        skipped_channel_snowflakes_by_guild_snowflake, skipped_member_snowflakes_by_guild_snowflake = {}, {}
        skipped_guild_snowflakes = set()
        target = 'user'
        title = f'{self.emoji.get_random_emoji()} Voice Mute(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list voice mutes across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            voice_mutes = await VoiceMute.fetch_all_by_target(target=target)
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope) 
                voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, target=target)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(interaction, scope) 
                    voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, target=target)
                except Exception as e:
                    if highest_role not in ('Owner', 'Developer', 'Administrator'):
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list voice mutes for specific servers.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=scope, target=target)
        else:
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id, target=target)
            channel_obj = interaction.channel

        if not voice_mutes:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No voice mute(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for voice_mute in voice_mutes:
            guild_dictionary.setdefault(voice_mute.guild_snowflake, {})
            guild_dictionary[voice_mute.guild_snowflake].setdefault(voice_mute.channel_snowflake, [])
            guild_dictionary[voice_mute.guild_snowflake][voice_mute.channel_snowflake].append({
                'member_snowflake': voice_mute.member_snowflake,
                'reason': voice_mute.reason,
                'expires_in': DurationObject.from_expires_in(voice_mute.expires_in)
            })
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for member_data in members:
                    member = guild.get_member(member_data['member_snowflake'])
                    if not member:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_data['member_snowflake'])
                        continue
                    entry_lines = [f'**User**: {member.mention}', f'**Expires in**: {member_data["expires_in"]}']
                    if member_obj:
                        entry_lines.append(f'**Reason**: {member_data["reason"]}')
                    if field_count >= chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(entry_lines), inline=False)
                    field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
        if skipped_member_snowflakes_by_guild_snowflake:
            for guild_snowflake, member_list in skipped_member_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for member_snowflake in member_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(member_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No voice mutes found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(name='mutes', help='Lists mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', channel ID/mention, member ID/mention, server ID or empty.")
    ):
        state = State(ctx)
        voice_mutes, pages = [], []
        channel_obj = None
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        skipped_channel_snowflakes_by_guild_snowflake, skipped_member_snowflakes_by_guild_snowflake = {}, {}
        skipped_guild_snowflakes = set()
        target = 'user'
        title = f'{self.emoji.get_random_emoji()} Voice Mute(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list voice mutes across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            voice_mutes = await VoiceMute.fetch_all_by_target(target=target)
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope) 
                voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, target=target)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(ctx, scope) 
                    voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, target=target)
                except Exception as e:
                    if highest_role not in ('Owner', 'Developer', 'Administrator'):
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list voice mutes for specific servers.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=scope, target=target)
        else:
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id, target=target)
            channel_obj = ctx.channel

        if not voice_mutes:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No voice mute(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for voice_mute in voice_mutes:
            guild_dictionary.setdefault(voice_mute.guild_snowflake, {})
            guild_dictionary[voice_mute.guild_snowflake].setdefault(voice_mute.channel_snowflake, [])
            guild_dictionary[voice_mute.guild_snowflake][voice_mute.channel_snowflake].append({
                'member_snowflake': voice_mute.member_snowflake,
                'reason': voice_mute.reason,
                'expires_in': DurationObject.from_expires_in(voice_mute.expires_in)
            })
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for member_data in members:
                    member = guild.get_member(member_data['member_snowflake'])
                    if not member:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_data['member_snowflake'])
                        continue
                    entry_lines = [f'**User**: {member.mention}', f'**Expires in**: {member_data["expires_in"]}']
                    if member_obj:
                        entry_lines.append(f'**Reason**: {member_data["reason"]}')
                    if field_count >= chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(entry_lines), inline=False)
                    field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
        if skipped_member_snowflakes_by_guild_snowflake:
            for guild_snowflake, member_list in skipped_member_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for member_snowflake in member_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(member_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No voice mutes found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @app_commands.command(name='mstage', description='Mute/unmute a member in the active stage.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    async def stage_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
        try:
            member_obj = await self.member_service.resolve_member(interaction, member)
            if not member_obj.voice:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')

            check_not_self(interaction, member_snowflake=member_obj.id)
            await has_equal_or_higher_role(interaction, channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, sender_snowflake=interaction.user.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
        if not stage:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No active stage found in {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            await member_obj.edit(mute=not member_obj.voice.mute)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been {'muted' if member_obj.voice.mute else 'unmuted'}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
                
    # DONE
    @commands.command(name='mstage', help='Mute/unmute a member in the active stage.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def stage_mute_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: ChannelSnowflake = commands.parameter(default=None, description="Tag a channel or include it's snowflake ID")
    ):
        state = State(ctx)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
        try:
            member_obj = await self.member_service.resolve_member(ctx, member)
            if not member_obj.voice:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            check_not_self(ctx, member_snowflake=member_obj.id)
            await has_equal_or_higher_role(ctx, channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, sender_snowflake=ctx.author.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
        if not stage:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No active stage found in {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            await member_obj.edit(mute=not member_obj.voice.mute)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been {'muted' if member_obj.voice.mute else 'unmuted'}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @app_commands.command(name='stages', description='Lists stage mute statistics.')
    @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, server ID or empty.")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_stages_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        field_count = 0
        pages = []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Stage(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list stages across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            stages = await Stage.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
                stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
                stages = [stage] if stage else []
            except Exception as e:
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all stages in a specific server.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, server ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                stages = await Stage.fetch_by_guild(guild_snowflake=int(scope))
        else:
            stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
            stages = [stage] if stage else []
            channel_obj = interaction.channel
        
        if not stages:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No stage(s) setup for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for stage in stages:
            guild_dictionary.setdefault(stage.guild_snowflake, {})
            guild_dictionary[stage.guild_snowflake].setdefault(stage.channel_snowflake, [])
            guild_dictionary[stage.guild_snowflake][stage.channel_snowflake]['expires_in'] = DurationObject.from_expires_in(stage.expires_in)
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                embed.add_field(name=f'Channel: {channel.mention}', value=f'**Expires in:** {channel_data['expires_in']}', inline=False)
                field_count += 1
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No stages found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
                 
    # DONE
    @commands.command(name='stages', help='Lists stage mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_stages_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', channel ID/mention, server ID or empty.")
    ):
        state = State(ctx)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        field_count = 0
        pages = []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Stage(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list stages across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            stages = await Stage.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
                stages = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all stages in a specific server.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, server ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                stages = await Stage.fetch_by_guild(guild_snowflake=int(scope))
        else:
            stages = await Stage.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            channel_obj = ctx.channel
        
        if not stages:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No stage(s) setup for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for stage in stages:
            guild_dictionary.setdefault(stage.guild_snowflake, {})
            guild_dictionary[stage.guild_snowflake].setdefault(stage.channel_snowflake, [])
            guild_dictionary[stage.guild_snowflake][stage.channel_snowflake]['expires_in'] = DurationObject.from_expires_in(stage.expires_in)
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                embed.add_field(name=f'Channel: {channel.mention}', value=f'**Expires in:** {channel_data['expires_in']}', inline=False)
                field_count += 1
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No stages found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @app_commands.command(name='tmutes', description='Lists text-mute statistics.')
    @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, member ID/mention, server ID or empty.")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_text_mutes_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        lines, pages, text_mutes = [], [], []
        skipped_channel_snowflakes_by_guild_snowflake, skipped_member_snowflakes_by_guild_snowflake = {}, {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Text Mute(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list text mutes across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            text_mutes = await TextMute.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope) 
                text_mutes = await TextMute.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(interaction, scope) 
                    text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
                except Exception as e:
                    if highest_role not in ('Owner', 'Developer', 'Administrator'):
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list text mutes for specific servers.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    text_mutes = await TextMute.fetch_by_guild(guild_snowflake=int(scope))
        else:
            text_mutes = await TextMute.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
        
        if not text_mutes:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No text mute(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for text_mute in text_mutes:
            guild_dictionary.setdefault(text_mute.guild_snowflake, {})
            guild_dictionary[text_mute.guild_snowflake].setdefault(text_mute.channel_snowflake, [])
            guild_dictionary[text_mute.guild_snowflake][text_mute.channel_snowflake].append({
                'member_snowflake': text_mute.member_snowflake,
                'reason': text_mute.reason,
                'expires_in': DurationObject.from_expires_in(text_mute.expires_in)
            })
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for member_data in members:
                    member = guild.get_member(member_data['member_snowflake'])
                    if not member:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_data['member_snowflake'])
                        continue
                    entry_lines = [f'**User**: {member.mention}', f'**Expires in**: {member_data["expires_in"]}']
                    if member_obj:
                        entry_lines.append(f'**Reason**: {member_data["reason"]}')
                    if field_count >= chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(entry_lines), inline=False)
                    field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
        if skipped_member_snowflakes_by_guild_snowflake:
            for guild_snowflake, member_list in skipped_member_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for member_snowflake in member_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(member_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No text mutes found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(name='tmutes', help='Lists text-mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', channel ID/mention, server ID or empty.")
    ):
        state = State(ctx)
        channel_obj = None
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        lines, pages, text_mutes = [], [], []
        skipped_channel_snowflakes_by_guild_snowflake, skipped_member_snowflakes_by_guild_snowflake = {}, {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Text Mute(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list text mutes across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            text_mutes = await TextMute.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope) 
                text_mutes = await TextMute.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(ctx, scope) 
                    text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                except Exception as e:
                    if highest_role not in ('Owner', 'Developer', 'Administrator'):
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list text mutes for specific servers.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    text_mutes = await TextMute.fetch_by_guild(guild_snowflake=int(scope))
        else:
            text_mutes = await TextMute.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
        
        if not text_mutes:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No text mute(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for text_mute in text_mutes:
            guild_dictionary.setdefault(text_mute.guild_snowflake, {})
            guild_dictionary[text_mute.guild_snowflake].setdefault(text_mute.channel_snowflake, [])
            guild_dictionary[text_mute.guild_snowflake][text_mute.channel_snowflake].append({
                'member_snowflake': text_mute.member_snowflake,
                'reason': text_mute.reason,
                'expires_in': DurationObject.from_expires_in(text_mute.expires_in)
            })
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))
        
        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, members in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                for member_data in members:
                    member = guild.get_member(member_data['member_snowflake'])
                    if not member:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_data['member_snowflake'])
                        continue
                    entry_lines = [f'**User**: {member.mention}', f'**Expires in**: {member_data["expires_in"]}']
                    if member_obj:
                        entry_lines.append(f'**Reason**: {member_data["reason"]}')
                    if field_count >= chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(entry_lines), inline=False)
                    field_count += 1
            pages.append(embed)
        if skipped_guild_snowflakes:
            embed = discord.Embed(title=title, description='\u200b', color=discord.Color.blue())
            lines = []
            for guild_snowflake in skipped_guild_snowflakes:
                if field_count >= chunk_size:
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                    embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                    lines = []
                    field_count = 0
                lines.append(str(guild_snowflake))
                field_count += 1
            embed.description = '\n'.join(lines)
            pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake})')
                field_count = 0
                lines = []
                for channel_snowflake in channel_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(channel_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
        if skipped_member_snowflakes_by_guild_snowflake:
            for guild_snowflake, member_list in skipped_member_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for member_snowflake in member_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(member_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No text mutes found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
async def setup(bot: DiscordBot):
    cog = ModeratorCommands(bot)
    await bot.add_cog(cog)

