''' commands.py

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
from discord import app_commands
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.role_service import RoleService
from vyrtuous.utils.administrator import Administrator
from vyrtuous.utils.all import All
from vyrtuous.utils.developer import Developer
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.moderator import Moderator
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.state import State
from vyrtuous.utils.temporary_room import TemporaryRoom
   
import discord

class EveryoneCommands(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.emoji = Emojis()
        self.message_service = MessageService(self.bot, self.bot.db_pool)
        self.channel_service = ChannelService()
        self.member_service = MemberService()
        self.role_service = RoleService()
        
    # DONE
    @app_commands.command(name='admins', description='Lists all members with server mute privileges in this server.')
    @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, server ID or empty.")
    async def list_administrators_app_command(
        self,
        interaction: discord.Interaction,
        scope:str
    ):
        state = State(interaction)
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        pages = []
        skipped_guild_snowflakes = set()
        skipped_member_snowflakes_by_guild_snowflake = {}
        skipped_role_snowflakes_by_guild_snowflake = {}
        title = f'{self.emoji.get_random_emoji()} Administrator(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list administrators across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            administrators = await Administrator.fetch_all()
        elif scope:
            try:
                member_obj = await self.member_service.resolve_member(interaction, scope) 
                administrator = await Administrator.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
                administrators = [administrator] if administrator else None
            except Exception as e:
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list administrators for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, role ID/mention, server ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                administrators = await Administrator.fetch_by_guild(guild_snowflake=int(scope))
        else:
            administrators = await Administrator.fetch_by_guild(guild_snowflake=interaction.guild.id)
            guild_obj = interaction.guild
        
        if not administrators:
            try:
                if guild_obj:
                    scope = guild_obj.name
                return await state.end(warning=f'\U000026A0\U0000FE0F No administrator(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for administrator in administrators:
            guild_dictionary.setdefault(administrator.guild_snowflake, {})
            guild_dictionary[administrator.guild_snowflake].setdefault(administrator.member_snowflake, [])
            guild_dictionary[administrator.guild_snowflake][administrator.member_snowflake].extend(administrator.role_snowflakes)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake in guild_dictionary.keys():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for member_snowflake, role_snowflakes in guild_dictionary[guild_snowflake].items():
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                member_obj = guild.get_member(member_snowflake)
                if not member_obj:
                    skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_snowflake)
                    continue
                role_mentions = []
                for role_snowflake in role_snowflakes:
                    role_obj = guild.get_role(role_snowflake)
                    if role_obj:
                        role_mentions.append(role_obj.mention)
                    else:
                        skipped_role_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(role_snowflake)
                if role_mentions:
                    lines = f"**Roles:** {'\n'.join(role_mentions)}"
                else:
                    lines = '\u200b'
                embed.add_field(name=f'{member_obj.name} ({member_snowflake})', value=lines, inline=False)
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
        if skipped_role_snowflakes_by_guild_snowflake:
            for guild_snowflake, role_list in skipped_role_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for role_snowflake in role_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(role_snowflake))
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
                return await state.end(warning=f'\U000026A0\U0000FE0F No administrators found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(name='admins', help='Lists all members with server mute privileges in this server.')
    async def list_administrators_text_command(
        self,
        ctx: commands.Context,
        *,
        scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', channel ID/mention, server ID or empty.")
    ):
        state = State(ctx)
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        pages = []
        skipped_guild_snowflakes = set()
        skipped_member_snowflakes_by_guild_snowflake = {}
        skipped_role_snowflakes_by_guild_snowflake = {}
        title = f'{self.emoji.get_random_emoji()} Administrator(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list administrators across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            administrators = await Administrator.fetch_all()
        elif scope:
            try:
                member_obj = await self.member_service.resolve_member(ctx, scope)
                administrator = await Administrator.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                administrators = [administrator] if administrator else None
            except Exception as e:
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list administrators for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                administrators = await Administrator.fetch_by_guild(guild_snowflake=int(scope))
        else:
            administrators = await Administrator.fetch_by_guild(guild_snowflake=ctx.guild.id)
            guild_obj = ctx.guild
        
        if not administrators:
            try:
                if guild_obj:
                    scope = guild_obj.name
                return await state.end(warning=f'\U000026A0\U0000FE0F No administrator(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for administrator in administrators:
            guild_dictionary.setdefault(administrator.guild_snowflake, {})
            guild_dictionary[administrator.guild_snowflake].setdefault(administrator.member_snowflake, [])
            guild_dictionary[administrator.guild_snowflake][administrator.member_snowflake].extend(administrator.role_snowflakes)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake in guild_dictionary.keys():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for member_snowflake, role_snowflakes in guild_dictionary[guild_snowflake].items():
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                member_obj = guild.get_member(member_snowflake)
                if not member_obj:
                    skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_snowflake)
                    continue
                role_mentions = []
                for role_snowflake in role_snowflakes:
                    role_obj = guild.get_role(role_snowflake)
                    if role_obj:
                        role_mentions.append(role_obj.mention)
                    else:
                        skipped_role_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(role_snowflake)
                if role_mentions:
                    lines = f"**Roles:** {'\n'.join(role_mentions)}"
                else:
                    lines = '\u200b'
                embed.add_field(name=f'{member_obj.name} ({member_snowflake})', value=lines, inline=False)
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
        if skipped_role_snowflakes_by_guild_snowflake:
            for guild_snowflake, role_list in skipped_role_snowflakes_by_guild_snowflake.items():
                embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake})')
                field_count = 0
                lines = []
                for role_snowflake in role_list:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(color=discord.Color.red(), title=f'Server ({guild_snowflake}) continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(role_snowflake))
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
                return await state.end(warning=f'\U000026A0\U0000FE0F No administrators found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @app_commands.command(name='coords', description='Lists coordinators for a specific voice channel, all, or a member.')
    @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, server ID or empty.")
    async def list_coordinators_app_command(
        self,
        interaction : discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        lines, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_member_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Coordinator(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list text mutes across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            coordinators = await Coordinator.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope) 
                coordinators = await Coordinator.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(interaction, scope)
                    coordinators = await Coordinator.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
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
                    coordinators = await Coordinator.fetch_by_guild(guild_snowflake=int(scope))
        else:
            coordinators = await Coordinator.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
            channel_obj = interaction.channel
        
        if not coordinators:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No coordinator(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for coordinator in coordinators:
            guild_dictionary.setdefault(coordinator.guild_snowflake, {})
            guild_dictionary[coordinator.guild_snowflake].setdefault(coordinator.channel_snowflake, [])
            guild_dictionary[coordinator.guild_snowflake][coordinator.channel_snowflake].append(coordinator.member_snowflake)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, member_snowflakes in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                member_lines = []
                for member_snowflake in member_snowflakes:
                    member_obj = guild.get_member(member_snowflake)
                    if not member_obj:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_snowflake)
                        continue
                    member_lines.append(f'{member_obj.mention} ({member_snowflake})')
                    if member_lines:
                        lines = '\n'.join(member_lines)
                    else:
                        lines = '\u200b'
                field_count += 1
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                embed.add_field(name=f'Channel: {channel.mention}', value=lines, inline=False)
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
                return await state.end(warning=f'\U000026A0\U0000FE0F No coordinators found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @commands.command(name='coords', help='Lists coordinators for a specific voice channel, all, or a member.')
    async def list_coordinators_text_command(
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
        lines, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_member_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Coordinator(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list text mutes across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            coordinators = await Coordinator.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope) 
                coordinators = await Coordinator.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(ctx, scope)
                    coordinators = await Coordinator.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
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
                    coordinators = await Coordinator.fetch_by_guild(guild_snowflake=int(scope))
        else:
            coordinators = await Coordinator.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            channel_obj = ctx.channel
        
        if not coordinators:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No coordinator(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for coordinator in coordinators:
            guild_dictionary.setdefault(coordinator.guild_snowflake, {})
            guild_dictionary[coordinator.guild_snowflake].setdefault(coordinator.channel_snowflake, [])
            guild_dictionary[coordinator.guild_snowflake][coordinator.channel_snowflake].append(coordinator.member_snowflake)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, member_snowflakes in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                member_lines = []
                for member_snowflake in member_snowflakes:
                    member_obj = guild.get_member(member_snowflake)
                    if not member_obj:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_snowflake)
                        continue
                    member_lines.append(f'{member_obj.mention} ({member_snowflake})')
                    if member_lines:
                        lines = '\n'.join(member_lines)
                    else:
                        lines = '\u200b'
                field_count += 1
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                embed.add_field(name=f'Channel: {channel.mention}', value=lines, inline=False)
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
                return await state.end(warning=f'\U000026A0\U0000FE0F No coordinators found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @app_commands.command(name='devs', description='Lists developers.')
    @app_commands.describe(scope="Specify one of: 'all', server ID or empty.")
    async def list_developers_app_command(
        self,
        interaction : discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        pages = []
        skipped_guild_snowflakes = set()
        skipped_member_snowflakes_by_guild_snowflake = {}
        title = f'{self.emoji.get_random_emoji()} Developer(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list developers across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            developers = await Developer.fetch_all()
        elif scope:
            try:
                member_obj = await self.member_service.resolve_member(interaction, scope) 
                developers = await Developer.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer', 'Developer'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list developers for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, server ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                developers = await Developer.fetch_by_guild(guild_snowflake=int(scope))
        else:
            developers = await Developer.fetch_by_guild(guild_snowflake=interaction.guild.id)
            guild_obj = interaction.guild
        
        if not developers:
            try:
                if guild_obj:
                    scope = guild_obj.name
                return await state.end(warning=f'\U000026A0\U0000FE0F No developer(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for developer in developers:
            guild_dictionary.setdefault(developer.guild_snowflake, {})
            guild_dictionary[developer.guild_snowflake].setdefault(developer.member_snowflake, [])

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake in guild_dictionary.keys():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for member_snowflake in guild_dictionary[guild_snowflake]:
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                member_obj = guild.get_member(member_snowflake)
                if not member_obj:
                    skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_snowflake)
                    continue
                embed.add_field(name=f'{member_obj.name} ({member_snowflake})', value=member_obj.mention, inline=False)
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
                return await state.end(warning=f'\U000026A0\U0000FE0F No developers found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='devs', help='Lists developers.')
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        *,
        scope: Optional[str] = commands.parameter(default=None, description="'all', a specific server or user mention/ID")
    ):
        state = State(ctx)
        guild_obj = None
        member_obj = None
        chunk_size = 18
        field_count = 0
        pages = []
        skipped_guild_snowflakes = set()
        skipped_member_snowflakes_by_guild_snowflake = {}
        title = f'{self.emoji.get_random_emoji()} Developer(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list developers across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            developers = await Developer.fetch_all()
        elif scope:
            try:
                member_obj = await self.member_service.resolve_member(ctx, scope) 
                developers = await Developer.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            except Exception as e:
                if highest_role not in ('Owner', 'Developer', 'Developer'):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list developers for specific servers.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f"\U000026A0\U0000FE0F Scope must be one of: 'all', channel ID/mention, member ID/mention, role ID/mention, server ID or empty. Received: {scope}.")
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                developers = await Developer.fetch_by_guild(guild_snowflake=int(scope))
        else:
            developers = await Developer.fetch_by_guild(guild_snowflake=ctx.guild.id)
            guild_obj = ctx.guild
        
        if not developers:
            try:
                if guild_obj:
                    scope = guild_obj.name
                return await state.end(warning=f'\U000026A0\U0000FE0F No developer(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for developer in developers:
            guild_dictionary.setdefault(developer.guild_snowflake, {})
            guild_dictionary[developer.guild_snowflake].setdefault(developer.member_snowflake, [])

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake in guild_dictionary.keys():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for member_snowflake in guild_dictionary[guild_snowflake]:
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                member_obj = guild.get_member(member_snowflake)
                if not member_obj:
                    skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_snowflake)
                    continue
                embed.add_field(name=f'{member_obj.name} ({member_snowflake})', value=member_obj.mention, inline=False)
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
                return await state.end(warning=f'\U000026A0\U0000FE0F No developers found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @app_commands.command(name='mods', description='Lists moderator statistics.')
    @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, server ID or empty.")
    async def list_moderators_app_command(
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
        lines, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_member_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Moderator(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list text mutes across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            moderators = await Moderator.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope) 
                moderators = await Moderator.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(interaction, scope)
                    moderators = await Moderator.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
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
                    moderators = await Moderator.fetch_by_guild(guild_snowflake=int(scope))
        else:
            moderators = await Moderator.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
            channel_obj = interaction.channel
        
        if not moderators:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No moderator(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for moderator in moderators:
            guild_dictionary.setdefault(moderator.guild_snowflake, {})
            guild_dictionary[moderator.guild_snowflake].setdefault(moderator.channel_snowflake, [])
            guild_dictionary[moderator.guild_snowflake][moderator.channel_snowflake].append(moderator.member_snowflake)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, member_snowflakes in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                member_lines = []
                for member_snowflake in member_snowflakes:
                    member_obj = guild.get_member(member_snowflake)
                    if not member_obj:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_snowflake)
                        continue
                    member_lines.append(f'{member_obj.mention} ({member_snowflake})')
                    if member_lines:
                        lines = '\n'.join(member_lines)
                    else:
                        lines = '\u200b'
                field_count += 1
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                embed.add_field(name=f'Channel: {channel.mention}', value=lines, inline=False)
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
                return await state.end(warning=f'\U000026A0\U0000FE0F No moderators found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='mods',help='Lists moderator statistics.')
    async def list_moderators_text_command(
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
        lines, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_member_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        title = f'{self.emoji.get_random_emoji()} Moderator(s)'

        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list text mutes across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            moderators = await Moderator.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope) 
                moderators = await Moderator.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                try:
                    member_obj = await self.member_service.resolve_member(ctx, scope)
                    moderators = await Moderator.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
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
                    moderators = await Moderator.fetch_by_guild(guild_snowflake=int(scope))
        else:
            moderators = await Moderator.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            channel_obj = ctx.channel
        
        if not moderators:
            try:
                if guild_obj:
                    scope = guild_obj.name
                elif channel_obj:
                    scope = channel_obj.mention
                return await state.end(warning=f'\U000026A0\U0000FE0F No moderator(s) exist for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        guild_dictionary = {}
        for moderator in moderators:
            guild_dictionary.setdefault(moderator.guild_snowflake, {})
            guild_dictionary[moderator.guild_snowflake].setdefault(moderator.channel_snowflake, [])
            guild_dictionary[moderator.guild_snowflake][moderator.channel_snowflake].append(moderator.member_snowflake)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, member_snowflakes in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                member_lines = []
                for member_snowflake in member_snowflakes:
                    member_obj = guild.get_member(member_snowflake)
                    if not member_obj:
                        skipped_member_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_snowflake)
                        continue
                    member_lines.append(f'{member_obj.mention} ({member_snowflake})')
                    if member_lines:
                        lines = '\n'.join(member_lines)
                    else:
                        lines = '\u200b'
                field_count += 1
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
                embed.add_field(name=f'Channel: {channel.mention}', value=lines, inline=False)
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
                return await state.end(warning=f'\U000026A0\U0000FE0F No moderators found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(name='roleid', description='Get the ID of a role by name in this server.')
    @app_commands.describe(role_name='The name of the role to look up')
    async def get_role_id_app_command(
        self,
        interaction: discord.Interaction,
        role_name: str
    ):
        state = State(interaction)
        role = discord.utils.get(interaction.guild.roles, name=role_name)
        if role:
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Role `{role.name}` has ID `{role.id}`.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No role named `{role_name}` found in this server.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
    # DONE
    @commands.command(name='roleid', help='Get the ID of a role by name in this server.')
    async def get_role_id(self, ctx: commands.Context, *, role_name: str):
        state = State(ctx)
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Role `{role.name}` has ID `{role.id}`.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No role named `{role_name}` found in this server.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
    # DONE
    @app_commands.command(name='survey', description='Survey moderators, developers, owners, and coordinators in the current or specified channel.')
    @app_commands.describe(channel='Tag a voice/stage channel')
    async def stage_survey_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        chunk_size = 18
        pages = []
        owners, developers, administrators, moderators, coordinators = [], [], [], [], []
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
        for member in channel_obj.members:
            try:
                if await member_is_owner(interaction.guild.id, member.id):
                    owners.append(member)
                    break  # stop after first match if desired
            except commands.CheckFailure:
                pass
            try:
                if await member_is_developer(interaction.guild.id, member.id):
                    developers.append(member)
                    break
            except commands.CheckFailure:
                pass
            try:
                if await member_is_administrator(interaction.guild.id, member.id):
                    administrators.append(member)
                    break
            except commands.CheckFailure:
                pass
            try:
                if await member_is_coordinator(channel_obj.id, interaction.guild.id, member.id):
                    coordinators.append(member)
                    break
            except commands.CheckFailure:
                pass
            try:
                if await member_is_moderator(channel_obj.id, interaction.guild.id, member.id):
                    moderators.append(member)
                    break
            except commands.CheckFailure:
                pass
        owners_chunks = [owners[i:i + chunk_size] for i in range(0, len(owners), chunk_size)]
        developers_chunks = [developers[i:i + chunk_size] for i in range(0, len(developers), chunk_size)]
        administrators_chunks = [administrators[i:i + chunk_size] for i in range(0, len(administrators), chunk_size)]
        coordinators_chunks = [coordinators[i:i + chunk_size] for i in range(0, len(coordinators), chunk_size)]
        moderators_chunks = [moderators[i:i + chunk_size] for i in range(0, len(moderators), chunk_size)]
        roles_chunks = [
            ('Owners', owners, owners_chunks),
            ('Developers', developers, developers_chunks),
            ('Administrators', administrators, administrators_chunks),
            ('Coordinators', coordinators, coordinators_chunks),
            ('Moderators', moderators, moderators_chunks),
        ]
        max_pages = max(len(c[2]) for c in roles_chunks)
        for page in range(max_pages):
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Survey results for {channel_obj.name}',
                description=f'Total surveyed: {len(channel_obj.members)}',
                color=discord.Color.blurple()
            )
            for role_name, role_list, chunks in roles_chunks:
                chunk = chunks[page] if page < len(chunks) else []
                embed.add_field(
                    name=f'{role_name} ({len(chunk)}/{len(role_list)})',
                    value=', '.join(u.mention for u in chunk) if chunk else '*None*',
                    inline=False
                )
            pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No authorized roles found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(name='survey', help='Survey moderators, developers, owners, and coordinators in the current or specified channel.')
    async def stage_survey_text_command(
        self,
        ctx: commands.Context,
        *,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        channel_obj = None
        chunk_size = 18
        pages = []
        owners, developers, administrators, moderators, coordinators = [], [], [], [], []
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
        for member in channel_obj.members:
            try:
                if await member_is_owner(ctx.guild.id, member.id):
                    owners.append(member)
                    break
            except commands.CheckFailure:
                pass
            try:
                if await member_is_developer(ctx.guild.id, member.id):
                    developers.append(member)
                    break
            except commands.CheckFailure:
                pass
            try:
                if await member_is_administrator(ctx.guild.id, member.id):
                    administrators.append(member)
                    break
            except commands.CheckFailure:
                pass
            try:
                if await member_is_coordinator(channel_obj.id, ctx.guild.id, member.id):
                    coordinators.append(member)
                    break
            except commands.CheckFailure:
                pass
            try:
                if await member_is_moderator(channel_obj.id, ctx.guild.id, member.id):
                    moderators.append(member)
                    break
            except commands.CheckFailure:
                pass
        owners_chunks = [owners[i:i + chunk_size] for i in range(0, len(owners), chunk_size)]
        developers_chunks = [developers[i:i + chunk_size] for i in range(0, len(developers), chunk_size)]
        administrators_chunks = [administrators[i:i + chunk_size] for i in range(0, len(administrators), chunk_size)]
        coordinators_chunks = [coordinators[i:i + chunk_size] for i in range(0, len(coordinators), chunk_size)]
        moderators_chunks = [moderators[i:i + chunk_size] for i in range(0, len(moderators), chunk_size)]
        roles_chunks = [
            ('Owners', owners, owners_chunks),
            ('Developers', developers, developers_chunks),
            ('Administrators', administrators, administrators_chunks),
            ('Coordinators', coordinators, coordinators_chunks),
            ('Moderators', moderators, moderators_chunks),
        ]
        max_pages = max(len(c[2]) for c in roles_chunks)
        for page in range(max_pages):
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Survey results for {channel_obj.name}',
                description=f'Total surveyed: {len(channel_obj.members)}',
                color=discord.Color.blurple()
            )
            for role_name, role_list, chunks in roles_chunks:
                chunk = chunks[page] if page < len(chunks) else []
                embed.add_field(
                    name=f'{role_name} ({len(chunk)}/{len(role_list)})',
                    value=', '.join(u.mention for u in chunk) if chunk else '*None*',
                    inline=False
                )
            pages.append(embed)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No authorized roles found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
async def setup(bot: DiscordBot):
    cog = EveryoneCommands(bot)
    await bot.add_cog(cog)

