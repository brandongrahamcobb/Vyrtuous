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
    @app_commands.describe(scope="'all', channel name/ID/mention, or user mention/ID")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_bans_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, scope)
            check_not_self(interaction, member_snowflake=member_obj.id)
        except Exception as e:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            except:
                channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all bans in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            bans = await Ban.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active bans found in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=interaction.guild.id, moderations=bans, moderation_type=Ban)
        elif member_obj:
            bans = await Ban.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active bans found for user {member_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, moderations=bans, moderation_type=Ban)
        elif channel_obj:
            bans = await Ban.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active bans found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=interaction.guild.id, moderations=bans, moderation_type=Ban)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='bans', description='Lists ban statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="'all', channel name/ID/mention, or user mention/ID")
    ):
        state = State(ctx)
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, scope)
            check_not_self(ctx, member_snowflake=member_obj.id)
        except Exception as e:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            except:
                channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all bans in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            bans = await Ban.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active bans found in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=bans, moderation_type=Ban)
        elif member_obj:
            bans = await Ban.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not banned in any channels.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, moderations=bans, moderation_type=Ban)
        elif channel_obj:
            bans = await Ban.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active bans found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=ctx.guild.id, moderations=bans, moderation_type=Ban)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @app_commands.command(name='caps', description="List active caps for a channel or all channels if 'all' is provided.")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    @app_commands.describe(scope="'all', channel name/ID/mention")
    async def list_caps_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        lines, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            title = f'{self.emoji.get_random_emoji()} Caps in All Servers'
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list caps across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            caps = await Cap.fetch_all()
            title = f'{self.emoji.get_random_emoji()} Caps for All Servers'
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope) 
                title = f'{self.emoji.get_random_emoji()} Caps in {channel_obj.name}'
                caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list caps for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                try:
                    guild_obj = self.bot.get_guild(scope)
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, guild ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    title = f'{self.emoji.get_random_emoji()} Caps for {guild_obj.name}'
                    caps = await Cap.fetch_by_guild(guild_snowflake=scope)
                except:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, guild ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            title = f'{self.emoji.get_random_emoji()} Caps for {interaction.channel.mention}'
            caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
            channel_obj = interaction.channel

        if not caps:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No caps setup for {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for cap in caps:
            guild_dictionary.setdefault(cap.guild_snowflake, {})
            guild_dictionary[cap.guild_snowflake].setdefault(cap.channel_snowflake, {})
            guild_dictionary[cap.guild_snowflake][cap.channel_snowflake][cap.moderation_type] = cap.duration

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))
        
        for guild_snowflake, channels in guild_dictionary.items():
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            field_count = 0
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                if field_count == chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                field_count += 1
                for moderation_type, duration in channel_data.items():
                    duration_obj = DurationObject.from_seconds(duration)
                    lines.append(f'  ↳ {moderation_type} ({duration_obj})')
                embed.add_field(name=channel.mention, value=f"{'\n'.join(lines)}", inline=False)
                field_count += 1
                if field_count == chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
            pages.append(embed)
        if skipped_guild_snowflakes:
            for guild_snowflake in skipped_guild_snowflakes:
                embed = discord.Embed(title=title, description=f'Skipped Guild {guild_snowflake}', color=discord.Color.red())
                pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                guild_label = f'Skipped Guild {guild_snowflake}'
                field_count = 0
                embed = discord.Embed(title=title, description=guild_label, color=discord.Color.red())
                for channel_snowflake in channel_list:
                    if field_count == chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=guild_label + ' continued...', color=discord.Color.red())
                        field_count = 0
                    embed.add_field(name=f'Skipped Channel {channel_snowflake}', value='\u200b', inline=False)
                    field_count += 1
                pages.append(embed)

        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @commands.command(name='caps', help="List active caps for a channel or all channels if 'all' is provided.")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="'all', channel name/ID/mention")
    ):
        state = State(ctx)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        lines, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            title = f'{self.emoji.get_random_emoji()} Caps in All Servers'
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list caps across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            caps = await Cap.fetch_all()
            title = f'{self.emoji.get_random_emoji()} Caps for All Servers'
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope) 
                title = f'{self.emoji.get_random_emoji()} Caps in {channel_obj.name}'
                caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list caps for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                try:
                    guild_obj = self.bot.get_guild(scope)
                    if not guild_obj:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, guild ID or empty. Received: {scope}.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    title = f'{self.emoji.get_random_emoji()} Caps for {guild_obj.name}'
                    caps = await Cap.fetch_by_guild(guild_snowflake=scope)
                except:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, guild ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            title = f'{self.emoji.get_random_emoji()} Caps for {ctx.channel.mention}'
            caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            channel_obj = ctx.channel
        
        if not caps:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No caps setup for {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for cap in caps:
            guild_dictionary.setdefault(cap.guild_snowflake, {})
            guild_dictionary[cap.guild_snowflake].setdefault(cap.channel_snowflake, {})
            guild_dictionary[cap.guild_snowflake][cap.channel_snowflake][cap.moderation_type] = cap.duration

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))
        
        for guild_snowflake, channels in guild_dictionary.items():
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            field_count = 0
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                if field_count == chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                field_count += 1
                for moderation_type, duration in channel_data.items():
                    duration_obj = DurationObject.from_seconds(duration)
                    lines.append(f'  ↳ {moderation_type} ({duration_obj})')
                embed.add_field(name=channel.mention, value=f"{'\n'.join(lines)}", inline=False)
                field_count += 1
                if field_count == chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
            pages.append(embed)
        if skipped_guild_snowflakes:
            for guild_snowflake in skipped_guild_snowflakes:
                embed = discord.Embed(title=title, description=f'Skipped Guild {guild_snowflake}', color=discord.Color.red())
                pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                guild_label = f'Skipped Guild {guild_snowflake}'
                field_count = 0
                embed = discord.Embed(title=title, description=guild_label, color=discord.Color.red())
                for channel_snowflake in channel_list:
                    if field_count == chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=guild_label + ' continued...', color=discord.Color.red())
                        field_count = 0
                    embed.add_field(name=f'Skipped Channel {channel_snowflake}', value='\u200b', inline=False)
                    field_count += 1
                pages.append(embed)
                
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
 
    # DONE
    @app_commands.command(name='cmds', description="List command aliases routed to a specific channel, temp room, or all channels if 'all' is provided.")
    @app_commands.describe(scope="'all', channel name/ID/mention")
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
        lines, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, scope)
        except:
            channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            title = f'{self.emoji.get_random_emoji()} Aliases in All Servers'
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners and developers can list all aliases in all guilds.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            aliases = await Alias.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope) 
                title = f'{self.emoji.get_random_emoji()} Aliases in {channel_obj.name}'
                aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list aliases for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(scope)
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, guild ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                title = f'{self.emoji.get_random_emoji()} Aliases in {guild_obj.name}'
                aliases = await Alias.fetch_by_guild(guild_snowflake=scope)
        else:
            title = f'{self.emoji.get_random_emoji()} Aliases in {interaction.channel.mention}'
            aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
            channel_obj = interaction.channel

        if not aliases:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No aliases found for {scope}.')
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
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            field_count = 0
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                if field_count == chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                field_count += 1
                for alias_type, alias_names in channel_data.items():
                    lines.append(f"{alias_type}\n" + '\n'.join(f"  ↳ {name}" for name in alias_names))
                    embed.add_field(name=channel.mention, value='\n'.join(lines), inline=False)
                    field_count += 1
                    if field_count == chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    lines = []
            pages.append(embed)
        if skipped_guild_snowflakes:
            for guild_snowflake in skipped_guild_snowflakes:
                embed = discord.Embed(title=title, description=f'Skipped Guild {guild_snowflake}', color=discord.Color.red())
                pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                guild_label = f'Skipped Guild {guild_snowflake}'
                field_count = 0
                embed = discord.Embed(title=title, description=guild_label, color=discord.Color.red())
                for channel_snowflake in channel_list:
                    if field_count == chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=guild_label + ' continued...', color=discord.Color.red())
                        field_count = 0
                    embed.add_field(name=f'Skipped Channel {channel_snowflake}', value='\u200b', inline=False)
                    field_count += 1
                pages.append(embed)

        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='cmds', help="List command aliases routed to a specific channel, temp room, or all channels if 'all' is provided.")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="'all', channel name/ID/mention, or temp room name")
    ):
        state = State(ctx)
        channel_obj = None
        guild_obj = None
        chunk_size = 18
        lines, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, scope)
        except:
            channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            title = f'{self.emoji.get_random_emoji()} Aliases in All Servers'
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners and developers can list all aliases in all guilds.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            aliases = await Alias.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope) 
                title = f'{self.emoji.get_random_emoji()} Aliases in {channel_obj.name}'
                aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list aliases for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(scope)
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, guild ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                title = f'{self.emoji.get_random_emoji()} Aliases in {guild_obj.name}'
                aliases = await Alias.fetch_by_guild(guild_snowflake=scope)
        else:
            title = f'{self.emoji.get_random_emoji()} Aliases in {ctx.channel.mention}'
            aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            channel_obj = ctx.channel

        if not aliases:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No aliases found for {scope}.')
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
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            field_count = 0
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                if field_count == chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                field_count += 1
                for alias_type, alias_names in channel_data.items():
                    lines.append(f"{alias_type}\n" + '\n'.join(f"  ↳ {name}" for name in alias_names))
                    embed.add_field(name=channel.mention, value='\n'.join(lines), inline=False)
                    field_count += 1
                    if field_count == chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                        field_count = 0
                    lines = []
            pages.append(embed)
        if skipped_guild_snowflakes:
            for guild_snowflake in skipped_guild_snowflakes:
                embed = discord.Embed(title=title, description=f'Skipped Guild {guild_snowflake}', color=discord.Color.red())
                pages.append(embed)
        if skipped_channel_snowflakes_by_guild_snowflake:
            for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                guild_label = f'Skipped Guild {guild_snowflake}'
                field_count = 0
                embed = discord.Embed(title=title, description=guild_label, color=discord.Color.red())
                for channel_snowflake in channel_list:
                    if field_count == chunk_size:
                        pages.append(embed)
                        embed = discord.Embed(title=title, description=guild_label + ' continued...', color=discord.Color.red())
                        field_count = 0
                    embed.add_field(name=f'Skipped Channel {channel_snowflake}', value='\u200b', inline=False)
                    field_count += 1
                pages.append(embed)

        try:
            return await state.end(success=pages)
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
        channel: AppChannelSnowflake
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
            has_equal_or_higher_role(interaction, channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=message.author.id, sender_snowflake=interaction.user.id)
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
            has_equal_or_higher_role(ctx, channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=message.author.id, sender_snowflake=ctx.author.id)
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
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, scope)
            try:
                check_not_self(interaction, member_snowflake=member_obj.id)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            except:
                channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all flags in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            flags = await Flag.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active flags found in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=interaction.guild.id, moderations=flags, moderation_type=Flag)
        elif member_obj:
            flags = await Flag.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not flagged in any channels.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, moderations=flags, moderation_type=Flag)
        elif channel_obj:
            flags = await Flag.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active flags found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=interaction.guild.id, moderations=flags, moderation_type=Flag)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @commands.command(name='flags', help='List flag statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="'all', channel name/ID/mention, or user mention/ID")
    ):
        state = State(ctx)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, scope)
            try:
                check_not_self(ctx, member_snowflake=member_obj.id)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            except:
                channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all flags in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            flags = await Flag.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active flags found in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=flags, moderation_type=Flag)
        elif member_obj:
            flags = await Flag.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not flagged in any channels.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, moderations=flags, moderation_type=Flag)
        elif channel_obj:
            flags = await Flag.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active flags found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=ctx.guild.id, moderations=flags, moderation_type=Flag)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @app_commands.command(name='ls', description='List users veganed as going vegan in this server.')
    @app_commands.describe(scope='A member or channel snowflake ID/mention')
    async def list_vegans_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, scope)
            try:
                check_not_self(interaction, member_snowflake=member_obj.id)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            except:
                channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all new vegans in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            vegans = await Vegan.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No new vegans found in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=interaction.guild.id, moderations=vegans, moderation_type=Vegan)
        elif member_obj:
            vegans = await Vegan.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No new vegans found for user {member_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, moderations=vegans, moderation_type=Vegan)
        elif channel_obj:
            vegans = await Vegan.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active vegans found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=interaction.guild.id, moderations=vegans, moderation_type=Vegan)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
                            
    # DONE
    @commands.command(name='ls', help='List users veganed as going vegan in this server.')
    async def list_members_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, scope)
            try:
                check_not_self(ctx, member_snowflake=member_obj.id)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            except:
                channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all new vegans in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            vegans = await Vegan.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No new vegans found in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=vegans, moderation_type=Vegan)
        elif member_obj:
            vegans = await Vegan.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No new vegans found for user {member_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, moderations=vegans, moderation_type=Vegan)
        elif channel_obj:
            vegans = await Vegan.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active vegans found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=ctx.guild.id, moderations=vegans, moderation_type=Vegan)
        try:
            return await state.end(success=pages)
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
                return await state.end(success=f'{self.emoji.get_random_emoji()} Temporary room {old_name} migrated to {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(warning=f'\U000026A0\U0000FE0F No temporary rooms found called {old_name} in {interaction.guild.name}.')
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
                return await state.end(success=f'{self.emoji.get_random_emoji()} Temporary room {old_name} migrated to {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(warning=f'\U000026A0\U0000FE0F No temporary rooms found called {old_name} in {ctx.guild.name}.')
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
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, scope)
            try:
                check_not_self(interaction, member_snowflake=member_obj.id)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            except:
                channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all voice mutes in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=interaction.guild.id, target='user')
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=interaction.guild.id, moderations=voice_mutes, moderation_type=VoiceMute)
        elif member_obj:
            voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, target='user')
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found for user {member_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, moderations=voice_mutes, moderation_type=VoiceMute)
        elif channel_obj:
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, target='user')
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=interaction.guild.id, moderations=voice_mutes, moderation_type=VoiceMute)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(name='mutes', help='Lists mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="'all', channel name/ID/mention, or user mention/ID")
    ):
        state = State(ctx)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, scope)
            try:
                check_not_self(ctx, member_snowflake=member_obj.id)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            except:
                channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all voice mutes in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=ctx.guild.id, target='user')
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=voice_mutes, moderation_type=VoiceMute)
        elif member_obj:
            voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, target='user')
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found for user {member_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, moderations=voice_mutes, moderation_type=VoiceMute)
        elif channel_obj:
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, target='user')
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=ctx.guild.id, moderations=voice_mutes, moderation_type=VoiceMute)
        try:
            return await state.end(success=pages)
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
            await self.message_service.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        try:
            member_obj = await self.member_service.resolve_member(interaction, member)
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
                return await state.end(warning=f'\U000026A0\U0000FE0F No active stage found.')
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
            await self.message_service.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        try:
            member_obj = await self.member_service.resolve_member(ctx, member)
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
                return await state.end(warning=f'\U000026A0\U0000FE0F No active stage found.')
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
    @app_commands.describe(scope="'all', channel name/ID/mention")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_stages_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        chunk_size = 18
        lines, pages = [], []
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, scope)
        except:
            channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners/devs can view all stages in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            stages = await Stage.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not stages:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active stages in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            for i in range(0, len(stages), chunk_size):
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Active Stages in {interaction.guild.name}',
                    color=discord.Color.purple()
                )
                for s in stages[i:i+chunk_size]:
                    ch = interaction.guild.get_channel(s.channel_snowflake)
                    ch_name = ch.mention
                    voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=interaction.guild.id, target='room')
                    for voice_mute in voice_mutes:
                        embed.add_field(name=ch_name, value=f'Active stage mutes: {voice_mute.member_snowflake}', inline=False)
                pages.append(embed)
        else:
            stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not stage:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active stages in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, target='room')
            initiator = interaction.guild.get_member(stage.member_snowflake)
            initiator_name = initiator.mention
            expires = DurationObject.from_expires_at(stage.expires_at) if stage.expires_at else 'No expiration'
            for m in voice_mutes:
                user = interaction.guild.get_member(m.member_snowflake)
                duration_str = DurationObject.from_expires_at(m.expires_at) if m.expires_at else 'No expiration'
                reason = m.reason or 'No reason provided'
                lines.append(f'• {user.mention} — {reason} — {duration_str}')
            description = f'Initiator: {initiator_name}\nStage Expires: {expires}\nActive stage mutes: {len(lines)}\n\n'+'\n'.join(lines)
            for i in range(0, len(description.splitlines()), chunk_size):
                embed=discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Stage info for {channel_obj.mention}',
                    description='\n'.join(description.splitlines()[i:i+chunk_size]),
                    color=discord.Color.purple()
                )
                pages.append(embed)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')          
    # DONE
    @commands.command(name='stages', help='Lists stage mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_stages_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="'all', channel name/ID/mention")
    ):
        state = State(ctx)
        channel_obj = None
        chunk_size = 18
        lines, pages = [], []
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, scope)
        except:
            channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners/devs can view all stages in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            stages = await Stage.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not stages:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active stages in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            for i in range(0, len(stages), chunk_size):
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Active Stages in {ctx.guild.name}',
                    color=discord.Color.purple()
                )
                for s in stages[i:i+chunk_size]:
                    ch = ctx.guild.get_channel(s.channel_snowflake)
                    ch_name = ch.mention
                    voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=ctx.guild.id, target='room')
                    for voice_mute in voice_mutes:
                        embed.add_field(name=ch_name, value=f'Active stage mutes: {voice_mute.member_snowflake}', inline=False)
                pages.append(embed)
        else:
            stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not stage:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active stages in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, target='room')
            initiator = ctx.guild.get_member(stage.member_snowflake)
            initiator_name = initiator.mention
            expires = DurationObject.from_expires_at(stage.expires_at) if stage.expires_at else 'No expiration'
            for m in voice_mutes:
                user = ctx.guild.get_member(m.member_snowflake)
                duration_str = DurationObject.from_expires_at(m.expires_at) if m.expires_at else 'No expiration'
                reason = m.reason or 'No reason provided'
                lines.append(f'• {user.mention} — {reason} — {duration_str}')
            description = f'Initiator: {initiator_name}\nStage Expires: {expires}\nActive stage mutes: {len(lines)}\n\n'+'\n'.join(lines)
            for i in range(0, len(description.splitlines()), chunk_size):
                embed=discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Stage info for {channel_obj.mention}',
                    description='\n'.join(description.splitlines()[i:i+chunk_size]),
                    color=discord.Color.purple()
                )
                pages.append(embed)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @app_commands.command(name='tmutes', description='Lists text-mute statistics.')
    @app_commands.describe(scope="'all', channel name/ID/mention, or user mention/ID")
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_text_mutes_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, scope)
            try:
                check_not_self(interaction, member_snowflake=member_obj.id)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            except:
                channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners, developers and administrators can list all text-mutes in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            text_mutes = await TextMute.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No users are currently text-muted in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=interaction.guild.id, moderations=text_mutes, moderation_type=TextMute)
        elif member_obj:
            text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not text-muted in any channels.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=interaction.guild.id, moderations=text_mutes, moderation_type=TextMute)
        elif channel_obj:
            text_mutes = await TextMute.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No text-muted users currently in {channel_obj.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, moderations=text_mutes, moderation_type=TextMute)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(name='tmutes', help='Lists text-mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description="'all', channel name/ID/mention, or user mention/ID")
    ):
        state = State(ctx)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, scope)
            try:
                check_not_self(ctx, member_snowflake=member_obj.id)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            except:
                channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners, developers and administrators can list all text-mutes in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            text_mutes = await TextMute.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No users are currently text-muted in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=text_mutes, moderation_type=TextMute)
        elif member_obj:
            text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not text-muted in any channels.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, moderations=text_mutes, moderation_type=TextMute)
        elif channel_obj:
            text_mutes = await TextMute.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No text-muted users currently in {channel_obj.name}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, moderations=text_mutes, moderation_type=TextMute)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
            
async def setup(bot: DiscordBot):
    cog = ModeratorCommands(bot)
    await bot.add_cog(cog)

