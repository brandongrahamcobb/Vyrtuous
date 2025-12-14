''' commands.py

    Copyright (C) 2024  github.com/brandongrahamcobb

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
import os
import random
import subprocess
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from typing import Any, List, Optional, Union
import discord

from discord.ext.commands import Command
from discord import app_commands
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.discord_message_service import DiscordMessageService, ChannelPaginator, Paginator, UserPaginator
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.backup import Backup
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.duration import Duration
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.reason import Reason
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.utils.setup_logging import logger
   
class Hybrid(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.emoji = Emojis()
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
        self.channel_service = ChannelService()
        self.member_service = MemberService()
        
    # DONE
    @app_commands.command(name='admins', description='Lists all members with server mute privileges in this guild.')
    async def list_server_muters_app_command(
        self,
        interaction: discord.Interaction
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in a server.')
        async with self.bot.db_pool.acquire() as conn:
            records = await conn.fetch('''SELECT discord_snowflake FROM users WHERE $1=ANY(server_muter_guild_ids) ORDER BY discord_snowflake''', interaction.guild.id)
            if not records:
                return await interaction.response.send_message(content=f'\U0001F6AB No admins found in {guild.name}.', allowed_mentions=discord.AllowedMentions.none())
            description_lines = []
            for record in records:
                uid = record['discord_snowflake']
                member_obj = await self.member_service.resolve_member(interaction, uid)
                description_lines.append(f'• {member_obj.display_name} — {member_obj.mention}' if member_obj else f'• User ID `{uid}` (not in guild)')
            chunk_size = 18
            pages = []
            for i in range(0,len(description_lines),chunk_size):
                chunk = description_lines[i:i+chunk_size]
                embed = discord.Embed(
                    title=f'🔑 Administrators in {interaction.guild.name}',
                    color=discord.Color.blurple()
                )
                embed.add_field(name='Admins', value='\n'.join(chunk), inline=False)
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
            
    # DONE
    @commands.command(name='admins', help='Lists all members with server mute privileges in this guild.')
    async def list_server_muters_text_command(
        self,
        ctx: commands.Context
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        async with self.bot.db_pool.acquire() as conn:
            records = await conn.fetch('''
                SELECT discord_snowflake
                FROM users
                WHERE $1 = ANY(server_muter_guild_ids)
                ORDER BY discord_snowflake
            ''', ctx.guild.id)
            if not records:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No admins found in {ctx.guild.name}.', allowed_mentions=discord.AllowedMentions.none())
            description_lines = []
            for record in records:
                uid = record['discord_snowflake']
                member_obj = await self.member_service.resolve_member(ctx, uid)
                if member_obj:
                    description_lines.append(f'• {member_obj.display_name} — {member_obj.mention}')
                else:
                    description_lines.append(f'• User ID `{uid}` (not in guild)')
            chunk_size = 18
            pages = []
            for i in range(0, len(description_lines), chunk_size):
                chunk = description_lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'🔑 Administrators in {ctx.guild.name}',
                    color=discord.Color.blurple()
                )
                embed.add_field(name='Admins', value='\n'.join(chunk), inline=False)
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
    
    # DONE
    @app_commands.command(name='coords', description='Lists coordinators for a specific voice channel, all, or a member.')
    @app_commands.describe(target='"all", member or channel name/ID/mention')
    async def list_coordinators_app_command(
        self,
        interaction : discord.Interaction,
        target : Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content="This command must be used in a server.")
        member_obj = await self.member_service.resolve_member(interaction, target)
        if member_obj:
            target = None
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        async with self.bot.db_pool.acquire() as conn:
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
            if target and target.lower() == 'all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await interaction.response.send_message(content='\U0001F6AB You are not authorized to list all coordinators.')
                query = 'SELECT unnest(coordinator_channel_ids) AS channel_id, discord_snowflake FROM users WHERE coordinator_channel_ids IS NOT NULL'
                rows = await conn.fetch(query)
                if not rows:
                    return await interaction.response.send_message(content='\U0001F6AB No coordinators found in any voice channels.')
                channel_map = defaultdict(list)
                for row in rows:
                    channel_map[row['channel_id']].append(row['discord_snowflake'])
                pages = []
                for ch_id, user_ids in sorted(channel_map.items()):
                    vc = interaction.guild.get_channel(ch_id)
                    vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
                    embed = discord.Embed(
                        title=f'🧭 Coordinators for {vc_name}',
                        color=discord.Color.gold())
                    for uid in user_ids:
                        m = interaction.guild.get_member(uid)
                        name = m.display_name if m else f'User ID {uid}'
                        embed.add_field(name=f'{interaction.guild.name}', value=f'• {name} (<@{uid}>)', inline=False)
                    pages.append(embed)
                paginator = UserPaginator(self.bot, interaction, pages)
                return await paginator.start()
            elif member_obj:
                query = 'SELECT coordinator_channel_ids FROM users WHERE discord_snowflake=$1'
                row = await conn.fetchrow(query, member_obj.id)
                channels = []
                if row:
                    channels.extend(row.get('coordinator_channel_ids') or [])
                if not channels:
                    return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.display_name} is not a coordinator in any channels.')
                channel_mentions = []
                for ch_id in channels:
                    if not ch_id:
                        continue
                    vc = interaction.guild.get_channel(ch_id)
                    channel_mentions.append(vc.mention if vc else f'Unknown Channel ({ch_id})')
                embeds = []
                chunk_size = 18
                for i in range(0, len(channel_mentions), chunk_size):
                    chunk = channel_mentions[i:i+chunk_size]
                    embed = discord.Embed(
                        title=f'🧭 {member_obj.display_name} is a coordinator in:',
                        description = '\n'.join(f'• {ch}' for ch in chunk),
                        color = discord.Color.gold()
                    )
                    embeds.append(embed)
                paginator = UserPaginator(self.bot, interaction, embeds)
                return await paginator.start()
            elif channel_obj:
                if channel_obj.type != discord.ChannelType.voice:
                    return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
                query = 'SELECT discord_snowflake FROM users WHERE $1=ANY(coordinator_channel_ids)'
                rows = list({r['discord_snowflake']: r for r in rows}.values())
                if not rows:
                    return await interaction.response.send_message(content=f'\U0001F6AB No coordinators found for {channel_obj.mention}.')
                lines = []
                for r in rows:
                    uid = r['discord_snowflake']
                    m = interaction.guild.get_member(uid)
                    if m:
                        lines.append(f'• {m.display_name} — <@{uid}>')
                if not lines:
                    return await interaction.response.send_message(content=f'\U0001F6AB No coordinators currently in {guild.name}.')
                pages = []
                chunk_size = 18
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i+chunk_size]
                    embed = discord.Embed(
                        title=f'🧭 Coordinators for {channel_obj.name}',
                        description='\n'.join(chunk),
                        color=discord.Color.gold()
                    )
                    pages.append(embed)
                paginator = UserPaginator(self.bot, interaction, pages)
                return await paginator.start()
        
    # DONE
    @commands.command(name='coords', help='Lists coordinators for a specific voice channel, all, or a member.')
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Voice channel name, mention, ID, "all", or member ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, target)
        if member_obj:
            target=None
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        async with ctx.bot.db_pool.acquire() as conn:
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
            if target and target.lower() == 'all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all coordinators.')
                query = 'SELECT unnest(coordinator_channel_ids) AS channel_id, discord_snowflake FROM users WHERE coordinator_channel_ids IS NOT NULL'
                rows = await conn.fetch(query)
                if not rows:
                    return await self.handler.send_message(ctx, content='\U0001F6AB No coordinators found in any voice channels.')
                channel_map = defaultdict(list)
                for row in rows: channel_map[row['channel_id']].append(row['discord_snowflake'])
                pages = []
                for ch_id, user_ids in sorted(channel_map.items()):
                    vc = ctx.guild.get_channel(ch_id)
                    vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
                    embed = discord.Embed(
                        title = f'🧭 Coordinators for {vc_name}',
                        color = discord.Color.gold()
                    )
                    for uid in user_ids:
                        m = ctx.guild.get_member(uid)
                        name = m.display_name if m else f'User ID {uid}'
                        embed.add_field(name=f'{ctx.guild.name}', value=f'• {name} (<@{uid}>)', inline=False)
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            elif member_obj:
                query = 'SELECT coordinator_channel_ids, FROM users WHERE discord_snowflake=$1'
                row = await conn.fetchrow(query, member_obj.id)
                channels = []
                if row:
                    channels.extend(row.get('coordinator_channel_ids') or [])
                if not channels:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.display_name} is not a coordinator in any channels.')
                channel_mentions  =[]
                for ch_id in channels:
                    if not ch_id:
                        continue
                    vc = ctx.guild.get_channel(ch_id)
                    channel_mentions.append(vc.mention if vc else f'Unknown Channel ({ch_id})')
                embeds = []
                chunk_size = 18
                for i in range(0, len(channel_mentions), chunk_size):
                    chunk = channel_mentions[i:i+chunk_size]
                    embed = discord.Embed(
                        title=f'🧭 {member_obj.display_name} is a coordinator in:',
                        description='\n'.join(f'• {ch}' for ch in chunk),
                        color=discord.Color.gold()
                    )
                    embeds.append(embed)
                paginator = Paginator(self.bot, ctx, embeds)
                return await paginator.start()
            elif channel_obj:
                if channel_obj.type != discord.ChannelType.voice:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
                query = 'SELECT discord_snowflake FROM users WHERE $1=ANY(coordinator_channel_ids)'
                rows = await conn.fetch(query, channel_obj.id)
                rows = list({r['discord_snowflake']: r for r in rows}.values())
                if not rows:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No coordinators found for {channel_obj.mention}.')
                lines = []
                for r in rows:
                    uid = r['discord_snowflake']
                    m = ctx.guild.get_member(uid)
                    if m:
                        lines.append(f'• {m.display_name} — <@{uid}>')
                if not lines:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No coordinators currently in {ctx.guild.name}.')
                pages = []
                chunk_size = 18
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i+chunk_size]
                    embed = discord.Embed(
                        title=f'🧭 Coordinators for {channel_obj.name}',
                        description='\n'.join(chunk),
                        color=discord.Color.gold()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
    
    # DONE
    @app_commands.command(name='devs', description='Lists developers.')
    @app_commands.describe(target='"all" or a member snowflake ID/mention')
    async def list_developers_app_command(
        self,
        interaction : discord.Interaction,
        target : Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        async with self.bot.db_pool.acquire() as conn:
            pages = []
            member_obj = await self.member_service.resolve_member(interaction, target)
            highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
            if target is None:
                rows = await conn.fetch(
                    'SELECT discord_snowflake FROM users WHERE $1 = ANY(developer_guild_ids)',
                    interaction.guild.id
                )
                if not rows:
                    return await interaction.response.send_message(content=f'\U0001F6AB No developers are configured in {interaction.guild.name}.')
                for row in rows:
                    user = interaction.guild.get_member(row['discord_snowflake'])
                    name = user.display_name if user else f'User ID {row["discord_snowflake"]}'
                    embed = discord.Embed(
                        title = f'Developer: {name}',
                        color = discord.Color.blurple()
                    )
                    pages.append(embed)
            elif target.lower() == 'all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await interaction.response.send_message(ctx, content='\U0001F6AB You are not authorized to list all developers.')
                rows = await conn.fetch(
                    'SELECT discord_snowflake, developer_guild_ids FROM users WHERE array_length(developer_guild_ids, 1) > 0'
                )
                if not rows:
                    return await interaction.response.send_message(content='\U0001F6AB No developers are configured.')
                for row in rows:
                    user = self.bot.get_user(row['discord_snowflake'])
                    name = user.name if user else f'User ID {row["discord_snowflake"]}'
                    guilds = [self.bot.get_guild(gid).name for gid in row['developer_guild_ids'] if self.bot.get_guild(gid)]
                    embed = discord.Embed(
                        title = f'Developer: {name}',
                        description = ', '.join(guilds) if guilds else 'No known guilds',
                        color = discord.Color.blurple()
                    )
                    pages.append(embed)
            else:
                if not member_obj:
                    return await interaction.response.send_message(content=f'\U0001F6AB Could not resolve a valid member from input: {target}.')
                row = await conn.fetchrow(
                    'SELECT developer_guild_ids FROM users WHERE discord_snowflake=$1',
                    member_obj.id
                )
                if not row or not row['developer_guild_ids']:
                    return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.mention} is not a developer in any guilds.', allowed_mentions=discord.AllowedMentions.none())
                guilds = [self.bot.get_guild(gid).name for gid in row['developer_guild_ids'] if self.bot.get_guild(gid)]
                embed = discord.Embed(
                    title = f'Developer guilds for {member_obj.display_name}',
                    description = ', '.join(guilds) if guilds else 'No known guilds',
                    color = discord.Color.blurple()
                )
                pages.append(embed)
        paginator = UserPaginator(self.bot, interaction, pages)
        return await paginator.start()
        
    # DONE
    @commands.command(name='devs', help='Lists developers.')
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", or user mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        async with self.bot.db_pool.acquire() as conn:
            pages = []
            member_obj = await self.member_service.resolve_member(ctx, target)
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
            if member_obj and member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot list in the guilds the bot is a developer.')
            elif target and target.lower() == 'all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all developers.')
                rows = await conn.fetch('SELECT discord_snowflake, developer_guild_ids FROM users WHERE array_length(developer_guild_ids, 1) > 0')
                if not rows: return await self.handler.send_message(ctx, content='\U0001F6AB No developers are configured.')
                for row in rows:
                    user = self.bot.get_user(row['discord_snowflake'])
                    name = user.name if user else f'User ID {row["discord_snowflake"]}'
                    guilds = [self.bot.get_guild(gid).name for gid in row['developer_guild_ids'] if self.bot.get_guild(gid)]
                    embed = discord.Embed(title=f'Developer: {name}', description=', '.join(guilds) if guilds else 'No known guilds', color=discord.Color.blurple())
                    pages.append(embed)
            elif member_obj:
                row = await conn.fetchrow('SELECT developer_guild_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
                if not row or not row['developer_guild_ids']:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not a developer in any guilds.', allowed_mentions=discord.AllowedMentions.none())
                guilds = [self.bot.get_guild(gid).name for gid in row['developer_guild_ids'] if self.bot.get_guild(gid)]
                embed = discord.Embed(title=f'Developer guilds for {member_obj.display_name}', description=', '.join(guilds) if guilds else 'No known guilds', color=discord.Color.blurple())
                pages.append(embed)
            else:
                rows = await conn.fetch('SELECT discord_snowflake FROM users WHERE $1 = ANY(developer_guild_ids)', ctx.guild.id)
                if not rows: return await self.handler.send_message(ctx, content=f'\U0001F6AB No developers are configured in {ctx.guild.name}.')
                for row in rows:
                    user = ctx.guild.get_member(row['discord_snowflake'])
                    name = user.display_name if user else f'User ID {row["discord_snowflake"]}'
                    embed = discord.Embed(title=f'Developer: {name}', color=discord.Color.blurple())
                    pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        return await paginator.start()
        
    # DONE
    @app_commands.command(name='mods', description='Lists moderator statistics.')
    @app_commands.describe(target='A member or channel snowflake ID/mention')
    async def list_moderators_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        member_obj = await self.member_service.resolve_member(interaction, target)
        if member_obj:
            target = None
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await interaction.response.send_message(content='\U0001F6AB You are not authorized to list all moderators.')
            query = '''SELECT unnest(moderator_channel_ids) AS channel_id, discord_snowflake FROM users WHERE moderator_channel_ids IS NOT NULL'''
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(query)
            if not rows:
                return await interaction.response.send_message(content='\U0001F6AB No moderators found in any voice channels.')
            channel_map = defaultdict(list)
            for row in rows:
                channel_map[row['channel_id']].append(row['discord_snowflake'])
            pages = []
            for ch_id, user_ids in sorted(channel_map.items()):
                vc = interaction.guild.get_channel(ch_id)
                vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(title=f'\U0001F6E1 Moderators for {vc_name}', color=discord.Color.magenta())
                for uid in user_ids:
                    m = interaction.guild.get_member(uid)
                    name = m.display_name if m else f'User ID {uid}'
                    embed.add_field(name=f'{interaction.guild.name}', value=f'• {name} (<@{uid}>)', inline=False)
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
        if member_obj:
            query = '''SELECT moderator_channel_ids FROM users WHERE discord_snowflake=$1'''
            async with self.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(query, member_obj.id)
            if not row or not row['moderator_channel_ids']:
                return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.display_name} is not a moderator in any channels.')
            channel_mentions = [interaction.guild.get_channel(ch_id).mention if interaction.guild.get_channel(ch_id) else f'Unknown Channel ({ch_id})' for ch_id in row['moderator_channel_ids']]
            chunk_size = 18
            pages = []
            for i in range(0, len(channel_mentions), chunk_size):
                chunk = channel_mentions[i:i + chunk_size]
                embed = discord.Embed(title=f'🛡️ {member_obj.display_name} moderates:', description='\n'.join(f'• {ch}' for ch in chunk), color=discord.Color.magenta())
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
        if channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
            query = '''SELECT discord_snowflake FROM users WHERE $1=ANY(moderator_channel_ids)'''
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(query, channel_obj.id)
            if not rows:
                return await interaction.response.send_message(content=f'\U0001F6AB No moderators found for {channel_obj.mention}.')
            lines = []
            for row in rows:
                uid = row['discord_snowflake']
                m = interaction.guild.get_member(uid)
                if not m:
                    continue
                lines.append(f'• {m.display_name} — <@{uid}>')
            if not lines:
                return await interaction.response.send_message(content=f'\U0001F6AB No moderators currently in {interaction.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(title=f'\U0001F6E1 Moderators for {channel_obj.name}', description='\n'.join(chunk), color=discord.Color.magenta())
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
        
    # DONE
    @commands.command(name='mods',help='Lists moderator statistics.')
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Voice channel name/mention/ID, "all", or member mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all moderators.')
            async with ctx.bot.db_pool.acquire() as conn:
                query = 'SELECT unnest(moderator_channel_ids) AS channel_id, discord_snowflake FROM users WHERE moderator_channel_ids IS NOT NULL'
                rows = await conn.fetch(query)
                if not rows:
                    return await self.handler.send_message(ctx, content='\U0001F6AB No moderators found in any voice channels.')
                channel_map = defaultdict(list)
                for row in rows:
                    channel_map[row['channel_id']].append(row['discord_snowflake'])
                pages = []
                for ch_id,user_ids in sorted(channel_map.items()):
                    vc = ctx.guild.get_channel(ch_id)
                    vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
                    embed = discord.Embed(
                        title=f'\U0001F6E1 Moderators for {vc_name}',
                        color=discord.Color.magenta()
                    )
                    for uid in user_ids:
                        m = ctx.guild.get_member(uid)
                        name = m.display_name if m else f'User ID {uid}'
                        embed.add_field(name=f'{ctx.guild.name}',value=f'• {name} (<@{uid}>)',inline=False)
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
        elif member_obj:
            query = 'SELECT moderator_channel_ids FROM users WHERE discord_snowflake=$1'
            row = await conn.fetchrow(query, member_obj.id)
            channels = []
            if row:
                channels.extend(row.get('moderator_channel_ids') or [])
            if not channels: return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.display_name} is not a moderator in any channels.')
            channel_mentions = []
            for ch_id in channels:
                if not ch_id: continue
                vc = ctx.guild.get_channel(ch_id)
                channel_mentions.append(vc.mention if vc else f'Unknown Channel ({ch_id})')
            embeds = []
            chunk_size = 18
            for i in range(0, len(channel_mentions), chunk_size):
                chunk = channel_mentions[i:i+chunk_size]
                embed = discord.Embed(
                    title=f'🛡️ {member_obj.display_name} moderates:',
                    description='\n'.join(f'• {ch}' for ch in chunk),
                    color=discord.Color.magenta()
                )
                embeds.append(embed)
            paginator = Paginator(self.bot, ctx, embeds)
            return await paginator.start()
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            query = 'SELECT discord_snowflake FROM users WHERE $1=ANY(moderator_channel_ids)'
            rows = await conn.fetch(query,channel_obj.id)
            rows = list({r['discord_snowflake']:r for r in rows}.values())
            if not rows:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No moderators found for {channel_obj.mention}.')
            lines = []
            for r in rows:
                uid = r['discord_snowflake']
                m = ctx.guild.get_member(uid)
                if not m: continue
                lines.append(f'• {m.display_name} — <@{uid}>')
            if not lines:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No moderators currently in {ctx.guild.name}.')
            pages = []
            chunk_size = 18
            for i in range(0,len(lines),chunk_size):
                chunk = lines[i:i+chunk_size]
                embed = discord.Embed(
                    title=f'\U0001F6E1 Moderators for {channel_obj.name}',
                    description='\n'.join(chunk),
                    color=discord.Color.magenta()
                )
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()

    # DONE
    @app_commands.command(name='owners', description='Show temporary room stats for "all", a channel, or a member.')
    @app_commands.describe(target='"all", a channel mention/ID, or a member mention/ID')
    async def temp_room_stats_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        async with self.bot.db_pool.acquire() as conn:
            channel_obj = await self.channel_service.resolve_channel(interaction, target)
            member_obj = await self.member_service.resolve_member(interaction, target)
            highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
            if target and target.lower() == 'all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await interaction.response.send_message(content='\U0001F6AB You are not authorized to list all owners.')
                rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild(interaction.guild)
                if not rooms:
                    return await interaction.response.send_message(content='\U0001F6AB No temporary rooms exist.')
                pages = []
                chunk = 12
                for i in range(0, len(rooms), chunk):
                    subset = rooms[i:i+chunk]
                    embed = discord.Embed(
                        title='📊 Temporary Rooms',
                        color=discord.Color.blurple()
                    )
                    for room in subset:
                        embed.add_field(
                            name=room.room_name,
                            value=f'• Owner: {room.room_owner.display_name}\n• Channel: {room.channel.mention}',
                            inline=False
                        )
                    pages.append(embed)
                paginator = UserPaginator(self.bot, interaction, pages)
                return await paginator.start()
            if member_obj and is_owner_or_dev:
                rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild_and_member(interaction.guild, member_obj)
                if not rooms:
                    return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.display_name} does not own any temporary rooms.')
                pages = []
                chunk = 12
                for i in range(0, len(rooms), chunk):
                    subset = rooms[i:i+chunk]
                    embed = discord.Embed(
                        title=f'📊 Temporary Rooms Owned by {member_obj.display_name}',
                        color=discord.Color.blurple()
                    )
                    for room in subset:
                        embed.add_field(
                            name=room.room_name,
                            value=f'• Channel: {room.channel.mention}',
                            inline=False
                        )
                    pages.append(embed)
                paginator = UserPaginator(self.bot, interaction, pages)
                return await paginator.start()
            if channel_obj and is_owner_or_dev:
                room = await TemporaryRoom.fetch_temporary_room_by_channel(channel_obj)
                if not room:
                    return await interaction.response.send_message(content=f'\U0001F6AB {channel_obj.mention} is not a temporary room.')
                embed = discord.Embed(
                    title=f'📊 Temporary Room Info for {channel_obj.name}',
                    color=discord.Color.blurple()
                )
                embed.add_field(name='Room Name', value=room.room_name, inline=False)
                embed.add_field(name='Owner', value=f'{room.room_owner.display_name} ({room.room_owner.mention})', inline=False)
                embed.add_field(name='Channel', value=channel_obj.mention, inline=False)
                return await interaction.response.send_message(embed=embed)
            return await interaction.response.send_message(content='\U0001F6AB Could not interpret the target. Provide "all", a channel, or a member.')

    # DONE
    @commands.command(name='owners', help='Show temporary room stats for "all", a channel, or a member.')
    async def temp_room_stats_text_command(
        self,
        ctx,
        target: Optional[str] = commands.parameter(default=None, description='"all", a channel mention/ID, or a member mention/ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        async with self.bot.db_pool.acquire() as conn:
            channel_obj = await self.channel_service.resolve_channel(ctx, target)
            member_obj = await self.member_service.resolve_member(ctx, target)
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
            if target and target.lower() == 'all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all owners.')
                rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild(ctx.guild)
                if not rooms:
                    return await self.handler.send_message(ctx, content='\U0001F6AB No temporary rooms exist.')
                pages = []
                chunk = 12
                for i in range(0, len(rooms), chunk):
                    subset = rooms[i:i+chunk]
                    embed = discord.Embed(
                        title='📊 Temporary Rooms',
                        color=discord.Color.blurple()
                    )
                    for room in subset:
                        embed.add_field(
                            name=room.room_name,
                            value=f'• Owner: {room.room_owner.display_name}\n• Channel: {room.channel.mention}',
                            inline=False
                        )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            if member_obj and is_owner_or_dev:
                rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild_and_member(ctx.guild, member_obj)
                if not rooms:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.display_name} does not own any temporary rooms.')
                pages = []
                chunk = 12
                for i in range(0, len(rooms), chunk):
                    subset = rooms[i:i+chunk]
                    embed = discord.Embed(
                        title=f'📊 Temporary Rooms Owned by {member_obj.display_name}',
                        color=discord.Color.blurple()
                    )
                    for room in subset:
                        embed.add_field(
                            name=r['room_name'],
                            value=f'• Channel: {room.channel.mention}',
                            inline=False
                        )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            if channel_obj and is_owner_or_dev:
                room = await TemporaryRoom.fetch_temporary_room_by_channel(channel_obj)
                if not room:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {channel_obj.mention} is not a temporary room.')
                embed = discord.Embed(
                    title=f'📊 Temporary Room Info for {channel_obj.name}',
                    color=discord.Color.blurple()
                )
                embed.add_field(name='Room Name', value=room.room_name, inline=False)
                embed.add_field(name='Owner', value=f"{room.room_owner.display_name} ({room.room_owner.mention})", inline=False)
                embed.add_field(name='Channel', value=channel_obj.mention, inline=False)
                return await self.handler.send_message(ctx, embed=embed)
            return await self.handler.send_message(ctx, content='\U0001F6AB Could not interpret the target. Provide "all", a channel, or a member.')
    
    # DONE
    @app_commands.command(name='roleid', description='Get the ID of a role by name in this server.')
    @app_commands.describe(role_name='The name of the role to look up')
    async def get_role_id_app_command(
        self,
        interaction: discord.Interaction,
        role_name: str
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in a server.')
        role = discord.utils.get(interaction.guild.roles, name=role_name)
        if role:
            await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Role `{role.name}` has ID `{role.id}`.')
        else:
            await interaction.response.send_message(content=f'\U0001F6AB No role named "{role_name}" found in this server.')

    # DONE
    @commands.command(name='roleid', help='Get the ID of a role by name in this server.')
    async def get_role_id(self, ctx: commands.Context, *, role_name: str):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in a server.')
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Role `{role.name}` has ID `{role.id}`.')
        else:
            await self.handler.send_message(ctx, content=f'\U0001F6AB No role named "{role_name}" found in this server.')
    
    # DONE
    @app_commands.command(name='survey', description='Survey moderators, developers, owners, and coordinators in the current or specified channel.')
    @app_commands.describe(channel='Tag a voice/stage channel')
    async def stage_survey_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if not isinstance(channel_obj, discord.VoiceChannel):
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        owners, developers, administrators, moderators, coordinators = [], [], [], []
        for member in channel_obj.members:
            if await member_is_owner(member): owners.append(member)
            elif await member_is_developer(member): developers.append(member)
            elif await member_is_administrator(member): administrators.append(member)
            elif await member_is_coordinator(channel_obj, member): coordinators.append(member)
            elif await imember_is_moderator(channel_obj, member): moderators.append(member)
        def fmt(users): return ', '.join(u.mention for u in users) if users else '*None*'
        msg = (
            f'\U0001F50D **Survey results for {channel_obj.mention}:**\n'
            f'\n**Owners:** {fmt(owners)}'
            f'\n**Developers:** {fmt(developers)}'
            f'\n**Administrators:** {fmt(administrators)}'
            f'\n**Coordinators:** {fmt(coordinators)}'
            f'\n**Moderators:** {fmt(moderators)}'
            f'\n\nTotal surveyed: {len(channel_obj.members)}'
        )
        await interaction.response.send_message(content=msg, allowed_mentions=discord.AllowedMentions.none())
    
    # DONE
    @commands.command(name='survey', help='Survey moderators, developers, owners, and coordinators in the current or specified channel.')
    @is_owner_developer_predicator()
    async def stage_survey_text_command(
        self,
        ctx: commands.Context,
        *,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if not isinstance(channel_obj, (discord.VoiceChannel, discord.StageChannel)):
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid voice or stage channel.')
        owners, developers, administrators, moderators, coordinators = [], [], [], []
        for member in channel_obj.members:
            if await member_is_owner(member): owners.append(member)
            elif await member_is_developer(member): developers.append(member)
            elif await member_is_administrator(member): administrators.append(member)
            elif await member_is_coordinator(channel_obj, member): coordinators.append(member)
            elif await imember_is_moderator(channel_obj, member): moderators.append(member)
        def fmt(users): return ', '.join(u.mention for u in users) if users else '*None*'
        msg = (
            f'\U0001F50D **Survey results for {channel_obj.mention}:**\n'
            f'\n**Owners:** {fmt(owners)}'
            f'\n**Developers:** {fmt(developers)}'
            f'\n**Administrators:** {fmt(administrators)}'
            f'\n**Coordinators:** {fmt(coordinators)}'
            f'\n**Moderators:** {fmt(moderators)}'
            f'\n\nTotal surveyed: {len(channel_obj.members)}'
        )
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())
        
async def setup(bot: DiscordBot):
    cog = Hybrid(bot)
    await bot.add_cog(cog)

