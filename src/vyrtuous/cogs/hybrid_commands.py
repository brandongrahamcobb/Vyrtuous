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
import asyncio
import discord
import inspect
import logging
import re
from collections import defaultdict
from discord import app_commands
from discord.ext import commands
from discord.ext.commands import Command
from vyrtuous.service.discord_message_service import DiscordMessageService, Paginator
from vyrtuous.service.check_service import *
from vyrtuous.inc.helpers import *
from vyrtuous.utils.setup_logging import logger
from types import MethodType
from typing import List, Optional, Union

logger = logging.getLogger(__name__)
class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
  
    #
    #   Helper method for loading aliases at runtime.
    #
    async def cog_load(self) -> None:
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT guild_id, alias_type, alias_name, channel_id FROM command_aliases')
            for row in rows:
                guild_id = row['guild_id']
                alias_type = row['alias_type']
                alias_name = row['alias_name']
                channel_id = row['channel_id']
                self.bot.command_aliases[guild_id][alias_type][alias_name] = channel_id
                if alias_type == 'mute':
                    cmd = self.create_mute_alias(alias_name)
                elif alias_type == 'unmute':
                    cmd = self.create_unmute_alias(alias_name)
                elif alias_type == 'ban':
                    cmd = self.create_ban_alias(alias_name)
                elif alias_type == 'unban':
                    cmd = self.create_unban_alias(alias_name)
                elif alias_type == 'role':
                    cmd = self.create_role_alias(alias_name)
                elif alias_type == 'unrole':
                    cmd = self.create_unrole_alias(alias_name)
                elif alias_type == 'flag':
                    cmd = self.create_flag_alias(alias_name)
                self.bot.add_command(cmd)
        
    #
    #  Help Command: helper method for the help command.
    #
    async def get_available_commands(self, bot, ctx) -> list[commands.Command]:
        available_commands = []
        for command in bot.commands:
            try:
                if await command.can_run(ctx):
                    available_commands.append(command)
            except commands.CheckFailure:
                continue
            except Exception as e:
                logger.warning(f"❌ Exception while checking command '{command}': {e}")
        return available_commands
    
    @commands.hybrid_command(name='coord', help='Grants coordinator access to a user in this guild.')
    @commands.check(is_owner_developer)
    async def create_coordinator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ) -> None:
        guild_id = ctx.guild.id
        member_id = None
        member_object = None
    
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
    
        if member_id:
            member_object = ctx.guild.get_member(member_id)
    
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
    
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (user_id, coordinator_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[])
                ON CONFLICT (user_id) DO UPDATE
                SET coordinator_ids = (
                    SELECT ARRAY(
                        SELECT DISTINCT unnest(u.coordinator_ids || EXCLUDED.coordinator_ids)
                    )
                    FROM users u WHERE u.user_id = EXCLUDED.user_id
                ),
                updated_at = NOW()
            ''', member_object.id, guild_id)
    
        await self.handler.send_message(ctx, content=f'{member_object.mention} has been granted coordinator rights in this guild.')

    @commands.hybrid_command(name='xcoord', help='Revokes coordinator access from a user in this guild.')
    @commands.check(is_owner_developer)
    async def delete_coordinator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ) -> None:
        guild_id = ctx.guild.id
        member_id = None
        member_object = None
    
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
    
        if member_id:
            member_object = ctx.guild.get_member(member_id)
    
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
    
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET coordinator_ids = array_remove(coordinator_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member_object.id, guild_id)
    
        await self.handler.send_message(ctx, content=f'{member_object.mention}\'s coordinator access has been revoked in this guild.')

    @commands.hybrid_command(name='coords', help='Lists coordinators for a specific voice channel.')
    @commands.check(is_owner_developer_coordinator)
    async def list_coordinators(
        self,
        ctx,
        channel: str = commands.parameter(description='Voice channel ID, mention, or name.')
    ) -> None:
        guild = ctx.guild
        def resolve_channel(value: str) -> Optional[discord.VoiceChannel]:
            if value.isdigit():
                channel = guild.get_channel(int(value))
            elif value.startswith('<#') and value.endswith('>'):
                channel = guild.get_channel(int(value.strip('<#>')))
            else:
                channel = discord.utils.get(guild.voice_channels, name=value)
            return channel if isinstance(channel, discord.VoiceChannel) else None
        voice_channel = resolve_channel(channel)
        if not voice_channel:
            return await self.handler.send_message(ctx, content='❌ Could not resolve a valid voice channel.')
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT user_id FROM users
                WHERE $1 = ANY(coordinator_ids)
            ''', voice_channel.id)
        if not rows:
            return await self.handler.send_message(ctx, content=f'ℹ️ No coordinators found for <#{voice_channel.id}>.')
        pages = []
        for row in rows:
            user_id = row['user_id']
            user = guild.get_member(user_id)
            name = user.display_name if user else f'User ID {user_id}'
            embed = discord.Embed(
                title=f'Coordinator: {name}',
                description=f'Assigned to <#{voice_channel.id}>',
                color=discord.Color.gold()
            )
            embed.set_footer(text=f'User ID: {user_id}')
            pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()


    #
    # Developer Commands: creation
    #                     deletion
    #                     listing
    #
    @commands.hybrid_command(name='dev', help='Elevates a user\'s permissions to a bot developer.')
    @commands.check(is_owner)
    async def create_developer(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ) -> None:
        guild_id = ctx.guild.id
        member_id = None
        member_object = None
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (user_id, developer_guild_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[])
                ON CONFLICT (user_id) DO UPDATE
                SET developer_guild_ids = (
                    SELECT ARRAY(
                        SELECT DISTINCT unnest(u.developer_guild_ids || EXCLUDED.developer_guild_ids)
                    )
                    FROM users u WHERE u.user_id = EXCLUDED.user_id
                ),
                updated_at = NOW()
            ''', member_object.id, guild_id)
        await self.handler.send_message(ctx, content=f'{member_object.mention} has been granted developer rights in this guild.')
        
    @commands.hybrid_command(name='xdev', help='Removes a developer.')
    @commands.check(is_owner)
    async def delete_developer(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ) -> None:
        member_id = None
        member_object = None
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        guild_id = ctx.guild.id
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET developer_guild_ids = array_remove(developer_guild_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member_object.id, guild_id)
        await self.handler.send_message(ctx, content=f'{member_object.mention}\'s developer access has been revoked in this guild.')

    @commands.hybrid_command(name='devs', help='Lists all developers.')
    @commands.check(is_owner_developer)
    async def list_developers(self, ctx) -> None:
        guild = ctx.guild
        pages = []
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT user_id, developer_guild_ids
                FROM users
                WHERE $1 = ANY(developer_guild_ids)
            ''', guild.id)
        if not rows:
            await self.handler.send_message(ctx, content='No developers are configured in this guild.')
            return
        for row in rows:
            user_id = row['user_id']
            user = guild.get_member(user_id)
            name = user.display_name if user else f'User ID {user_id}'
            embed = discord.Embed(
                title=f'Developer: {name}',
                color=discord.Color.blue()
            )
            pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()
    
    #
    # Moderator Commands: creation
    #                 deletion
    #                 listing
    #
    @commands.hybrid_command(name='mod', help='Elevates a user\'s permission to VC moderator for a specific channel.')
    @commands.check(is_owner_developer_coordinator)
    async def create_moderator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel: str = commands.parameter(description='Tag a channel or include its snowflake ID.')
    ) -> None:
        member_id = None
        member_object = None
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        resolved_channel = None
        if channel.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel))
        elif channel.startswith('<#') and channel.endswith('>'):
            try:
                channel_id = int(channel.strip('<#>'))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel.lower():
                    resolved_channel = vc
                    break
        if not isinstance(resolved_channel, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content='Could not resolve a valid **voice** channel from your input.')
        async with self.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (user_id, moderator_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[])
                ON CONFLICT (user_id) DO UPDATE
                SET moderator_ids = (
                    SELECT ARRAY(
                        SELECT DISTINCT unnest(u.moderator_ids || EXCLUDED.moderator_ids)
                    )
                    FROM users u WHERE u.user_id = EXCLUDED.user_id
                ),
                updated_at = NOW()
            ''', member_object.id, resolved_channel.id)
        await self.handler.send_message(ctx, content=f'{member_object.mention} has been granted VC moderator access in {resolved_channel.name}.')
    
    @commands.hybrid_command(name='xmod', help='Revokes a member\'s VC moderator role for a given channel.')
    @commands.check(is_owner_developer_coordinator_moderator)
    async def delete_moderator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel: str = commands.parameter(description='Tag a VC or include its snowflake ID.')
    ) -> None:
        member_id = None
        member_object = None
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        resolved_channel = None
        if channel.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel))
        elif channel.startswith('<#') and channel.endswith('>'):
            try:
                channel_id = int(channel.strip('<#>'))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel.lower():
                    resolved_channel = vc
                    break
        if not isinstance(resolved_channel, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content='Could not resolve a valid **voice** channel from your input.')
        async with self.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET moderator_ids = array_remove(moderator_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member_object.id, resolved_channel.id)
        await self.handler.send_message(ctx, content=f'{member_object.mention} has been revoked moderator access in {resolved_channel.name}.')

    @commands.hybrid_command(name='mods', help='Lists moderators for a specific voice channel.')
    @commands.check(is_owner_developer_coordinator_moderator)
    async def list_moderators(
        self,
        ctx,
        channel: str = commands.parameter(description='Voice channel ID, mention, or name.')
    ) -> None:
        guild = ctx.guild
        def resolve_channel(value: str) -> Optional[discord.VoiceChannel]:
            if value.isdigit():
                channel = guild.get_channel(int(value))
            elif value.startswith('<#') and value.endswith('>'):
                channel = guild.get_channel(int(value.strip('<#>')))
            else:
                channel = discord.utils.get(guild.voice_channels, name=value)
            return channel if isinstance(channel, discord.VoiceChannel) else None
        voice_channel = resolve_channel(channel)
        if not voice_channel:
            return await self.handler.send_message(ctx, content='❌ Could not resolve a valid voice channel.')
        channel_id = voice_channel.id
    
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT user_id FROM users
                WHERE $1 = ANY(moderator_ids)
            ''', channel_id)
        if not rows:
            return await self.handler.send_message(ctx, content=f'ℹ️ No moderators found for <#{channel_id}>.')
        pages = []
        for row in rows:
            user_id = row['user_id']
            user = guild.get_member(user_id)
            display_name = user.display_name if user else f'User ID {user_id}'
    
            embed = discord.Embed(
                title=f'Moderator for {voice_channel.name}',
                description=f'{display_name} (<@{user_id}>)',
                color=discord.Color.green()
            )
            embed.set_footer(text=f'User ID: {user_id}')
            pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()

        
    @commands.hybrid_command(
        name='flags',
        help='List users flagged in the database for a specific voice channel.'
    )
    @commands.check(is_owner_developer_coordinator_moderator)
    async def flags(
        self,
        ctx,
        channel: str = commands.parameter(description='Voice channel ID, mention, or name.')
    ) -> None:
        guild = ctx.guild
        def resolve_channel(value: str) -> Optional[discord.VoiceChannel]:
            if value.isdigit():
                channel = guild.get_channel(int(value))
            elif value.startswith('<#') and value.endswith('>'):
                channel = guild.get_channel(int(value.strip('<#>')))
            else:
                channel = discord.utils.get(guild.voice_channels, name=value)
            return channel if isinstance(channel, discord.VoiceChannel) else None
        channel = resolve_channel(channel)
        if not channel:
            return await self.handler.send_message(ctx, content='❌ Could not resolve a valid voice channel.')
        channel_id = channel.id
        sql = '''
            SELECT user_id FROM users
            WHERE $1 = ANY(flagged_channel_ids)
        '''
        try:
            async with self.db_pool.acquire() as conn:
                rows = await conn.fetch(sql, channel_id)
            if not rows:
                return await self.handler.send_message(ctx, content=f'✅ No users are flagged for <#{channel_id}>.')
            mentions = [f'<@{row["user_id"]}>' for row in rows]
            message = f'🚩 Users flagged for <#{channel_id}>:\n' + '\n'.join(mentions)
            await self.handler.send_message(ctx, content=message)
        except Exception as e:
            await self.handler.send_message(ctx, content=f'❌ Database error: {e}')
            raise


    #
    # Alias Commands: creation
    #                 deletion
    #                 listing
    #
    @commands.hybrid_command(
        name='alias',
        help='Set an alias for a mute, unmute, ban, unban or flag action.'
    )
    @commands.check(is_owner_developer_coordinator)
    async def create_alias(
        self,
        ctx,
        alias_type: str = commands.parameter(description='One of: `mute`, `unmute`, `ban`, `unban`, `flag`'),
        alias_name: str = commands.parameter(description='Alias/Pseudonym'),
        channel_id: str = commands.parameter(description='Voice channel')
    ) -> None:
        alias_type = alias_type.lower()
        guild_id = ctx.guild.id
        valid_types = {'mute', 'unmute', 'ban', 'unban', 'flag'}
        if alias_type not in valid_types:
            return await self.handler.send_message(ctx, f'❌ Invalid alias type. Must be one of: {", ".join(valid_types)}', ephemeral=True)
        if not alias_name.strip():
            return await self.handler.send_message(ctx, '❌ Alias name cannot be empty.', ephemeral=True)
        def resolve_channel(value: str):
            if value.isdigit():
                return ctx.guild.get_channel(int(value))
            if value.startswith('<#') and value.endswith('>'):
                return ctx.guild.get_channel(int(value.strip('<#>')))
            return discord.utils.get(ctx.guild.voice_channels, name=value)
        channel = resolve_channel(channel_id)
        if not channel or channel.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, '❌ Could not resolve a valid voice channel.', ephemeral=True)
        async with self.db_pool.acquire() as conn:
            existing_alias = await conn.fetchrow(
                '''
                SELECT channel_id FROM command_aliases 
                WHERE guild_id = $1 AND alias_type = $2 AND alias_name = $3
                ''',
                guild_id, alias_type, alias_name
            )
            if existing_alias:
                existing_channel = ctx.guild.get_channel(existing_alias['channel_id'])
                channel_mention = existing_channel.mention if existing_channel else f"<#{existing_alias['channel_id']}>"
                return await self.handler.send_message(
                    ctx,
                    f'❌ Alias `{alias_name}` ({alias_type}) already exists and is set to {channel_mention}.',
                    ephemeral=True
                )
            if self.bot.get_command(alias_name):
                return await self.handler.send_message(
                    ctx,
                    f'❌ A command named `{alias_name}` already exists.',
                    ephemeral=True
                )
            await conn.execute(
                '''
                INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id)
                VALUES ($1, $2, $3, $4)
                ''',
                guild_id, alias_type, alias_name, channel.id
            )
        self.bot.command_aliases.setdefault(guild_id, {}).setdefault(alias_type, {})[alias_name] = channel.id
        if alias_type == 'mute':
            cmd = self.create_mute_alias(alias_name)
        elif alias_type == 'unmute':
            cmd = self.create_unmute_alias(alias_name)
        elif alias_type == 'ban':
            cmd = self.create_ban_alias(alias_name)
        elif alias_type == 'unban':
            cmd = self.create_unban_alias(alias_name)
        elif alias_type == 'flag':
            cmd = self.create_flag_alias(alias_name)
        self.bot.add_command(cmd)
        await self.handler.send_message(
            ctx,
            content=f'✅ Alias `{alias_name}` ({alias_type}) set to {channel.mention}.'
        )
            
    def create_ban_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(
            name=command_name,
            help='Ban a user from a voice channel. Use 0 hours for permanent ban (coordinators only).'
        )
        @commands.check(is_owner_developer_coordinator_moderator)
        async def ban_alias(
            ctx,
            member: str = commands.parameter(description='Mention or user ID of the member to ban.'),
            duration_hours: Optional[int] = commands.parameter(default=6, description='Duration of ban in hours (0 = permanent).'),
            *,
            reason: str = commands.parameter(default='', description='Reason for ban (required for permanent).')
        ) -> None:
            command_name = ctx.invoked_with
            guild_id = ctx.guild.id
            logger.info(f"[{guild_id}] Ban command invoked: {command_name} by {ctx.author} (User ID: {ctx.author.id})")
            member_id = None
            member_object = None
            if member.isdigit():
                member_id = int(member)
            elif member.startswith('<@') and member.endswith('>'):
                try:
                    member_id = int(member.strip('<@!>'))
                except ValueError:
                    pass
            member_object = ctx.guild.get_member(member_id) if member_id else None
            if not member_object:
                return await self.handler.send_message(ctx, content='❌ Could not resolve a valid member.')
            if duration_hours == 0:
                if not await is_owner_developer_coordinator(ctx):
                    return await self.handler.send_message(ctx, content='❌ Only coordinators can issue permanent bans.')
                if not reason.strip():
                    return await self.handler.send_message(ctx, content='❌ Reason is required for permanent bans.')
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('ban', {}).get(command_name)
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'❌ No channel alias mapping found for `{command_name}`.')
            role_id = None
            ban_channel_id = static_channel_id
            async with self.db_pool.acquire() as conn:
                role_id = await conn.fetchval(
                    '''
                    SELECT role_id FROM channel_roles
                    WHERE guild_id = $1 AND channel_id = $2
                    ''',
                    guild_id, static_channel_id
                )
                if not role_id:
                    return await self.handler.send_message(ctx, content=f'❌ No ban role alias mapping found for `{command_name}`.')
            async with self.db_pool.acquire() as conn:
                existing_ban = await conn.fetchval(
                    '''
                    SELECT 1 FROM active_bans
                    WHERE user_id = $1 AND channel_id = $2
                    ''',
                    member_object.id, static_channel_id
                )
                if existing_ban:
                    return await self.handler.send_message(
                        ctx,
                        content=f'ℹ️ {member_object.mention} is already banned from <#{static_channel_id}>.'
                    )
            role = ctx.guild.get_role(role_id)
            if not role:
                return await self.handler.send_message(ctx, content=f'❌ Could not resolve role ID `{role_id}`.')
            try:
                await member_object.add_roles(role, reason=f"Banned from <#{static_channel_id}>: {reason or 'No reason provided'}")
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content='❌ Missing permissions to assign ban role.')
            try:
                async with self.db_pool.acquire() as conn:
                    await conn.execute(
                        '''
                        INSERT INTO active_bans (user_id, channel_id, expires_at)
                        VALUES ($1, $2, $3)
                        ON CONFLICT (user_id, channel_id) DO UPDATE SET expires_at = EXCLUDED.expires_at
                        ''',
                        member_object.id, static_channel_id,
                        None if duration_hours == 0 else discord.utils.utcnow() + timedelta(hours=duration_hours)
                    )
                    await conn.execute(
                        '''
                        UPDATE users
                        SET ban_channel_ids = array_append(ban_channel_ids, $2),
                            updated_at = NOW()
                        WHERE user_id = $1
                        ''',
                        member_object.id, static_channel_id
                    )
                    await conn.execute(
                        '''
                        INSERT INTO ban_reasons (guild_id, user_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4)
                        ON CONFLICT (guild_id, user_id, channel_id)
                        DO UPDATE SET reason = EXCLUDED.reason
                        ''',
                        guild_id, member_object.id, static_channel_id, reason or 'No reason provided'
                    )
            except Exception as e:
                logger.warning(f"🔥 Database error occurred: {e}")
            await self.handler.send_message(
                ctx,
                content=f'🔇 {member_object.mention} has been banned from <#{static_channel_id}> {"permanently" if duration_hours == 0 else f"for {duration_hours} hour(s)"}.'
            )
        
            if duration_hours > 0:
                expires_at = datetime.datetime.utcnow() + datetime.timedelta(hours=duration_hours)
                async with self.db_pool.acquire() as conn:
                    await conn.execute(
                        '''
                        INSERT INTO ban_expirations (user_id, channel_id, expires_at)
                        VALUES ($1, $2, $3)
                        ON CONFLICT (user_id, channel_id)
                        DO UPDATE SET expires_at = EXCLUDED.expires_at
                        ''',
                        member_object.id, static_channel_id, expires_at
                    )
        return ban_alias
        
    @commands.hybrid_command(name='bans', help='Lists all users banned from a specific channel.')
    @commands.check(is_owner_developer)
    async def list_bans(
        self,
        ctx: commands.Context,
        channel: str = commands.parameter(description='Mention or ID of the channel to check bans for.')
    ) -> None:
        guild = ctx.guild
        channel_id = None
        if channel.isdigit():
            channel_id = int(channel)
        elif channel.startswith('<#') and channel.endswith('>'):
            try:
                channel_id = int(channel.strip('<#>'))
            except ValueError:
                pass
        if not channel_id or not guild.get_channel(channel_id):
            return await self.handler.send_message(ctx, content='❌ Could not resolve a valid channel.')
        async with self.db_pool.acquire() as conn:
            bans = await conn.fetch('''
                SELECT ab.user_id, ab.expires_at, br.reason
                FROM active_bans ab
                LEFT JOIN ban_reasons br ON ab.user_id = br.user_id AND ab.channel_id = br.channel_id AND br.guild_id = $1
                WHERE ab.channel_id = $2
                ORDER BY ab.expires_at NULLS LAST
            ''', guild.id, channel_id)
        if not bans:
            return await self.handler.send_message(ctx, content=f'✅ No active bans found for <#{channel_id}>.')
        embeds = []
        for record in bans:
            user = guild.get_member(record['user_id'])
            name = user.mention if user else f"User ID {record['user_id']}"
            reason = record['reason'] or 'No reason provided'
            expires_at = record['expires_at']
            duration_str = (
                "Permanent" if expires_at is None
                else discord.utils.format_dt(expires_at, style='R')  # e.g., "in 5 minutes"
            )
            embed = discord.Embed(
                title="⛔ Channel Ban",
                description=f"{name} is banned from <#{channel_id}>.",
                color=discord.Color.red()
            )
            embed.add_field(name="Reason", value=reason, inline=False)
            embed.add_field(name="Duration", value=duration_str, inline=True)
            embed.set_footer(text=f"Channel ID: {channel_id}")
            embeds.append(embed)
        paginator = Paginator(self.bot, ctx, embeds)
        await paginator.start()

    def create_mute_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(name=command_name, help='Mutes a member in a specific VC.')
        @commands.check(is_owner_developer_coordinator_moderator)
        async def mute_alias(
            ctx,
            member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Optionally include a reason for the mute.')
        ) -> None:
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member.isdigit():
                member_id = int(member)
            elif member.startswith('<@') and member.endswith('>'):
                try:
                    member_id = int(member.strip('<@!>'))
                except ValueError:
                    pass
            if member_id:
                member_object = ctx.guild.get_member(member_id)
            if not member_object:
                return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('mute', {}).get(command_name)
            is_owner = await self.bot.is_owner(ctx.author)
            mute_source = 'owner' if is_owner else 'bot'
            async with self.db_pool.acquire() as conn:
                await conn.execute('''
                    INSERT INTO active_mutes (user_id, channel_id, source, issuer_id)
                    VALUES ($1, $2, $3, $4)
                    ON CONFLICT (user_id, channel_id) DO UPDATE 
                    SET source = EXCLUDED.source, issuer_id = EXCLUDED.issuer_id
                ''', member_object.id, static_channel_id, mute_source, ctx.author.id)
                if is_owner:
                    await conn.execute('''
                        INSERT INTO users (user_id, server_mute_channel_ids)
                        VALUES ($1, ARRAY[$2]::BIGINT[])
                        ON CONFLICT (user_id) DO UPDATE
                        SET server_mute_channel_ids = (
                            SELECT ARRAY(
                                SELECT DISTINCT unnest(u.server_mute_channel_ids || ARRAY[$2])
                            )
                            FROM users u WHERE u.user_id = EXCLUDED.user_id
                        ),
                        updated_at = NOW()
                    ''', member_object.id, static_channel_id)
                else:
                    await conn.execute('''
                        INSERT INTO users (user_id, mute_channel_ids)
                        VALUES ($1, ARRAY[$2]::BIGINT[])
                        ON CONFLICT (user_id) DO UPDATE
                        SET mute_channel_ids = (
                            SELECT ARRAY(
                                SELECT DISTINCT unnest(u.mute_channel_ids || ARRAY[$2])
                            )
                            FROM users u WHERE u.user_id = EXCLUDED.user_id
                        ),
                        updated_at = NOW()
                    ''', member_object.id, static_channel_id)
                await conn.execute('''
                    INSERT INTO mute_reasons (guild_id, user_id, reason, channel_id)
                    VALUES ($1, $2, $3, $4)
                    ON CONFLICT (guild_id, user_id, channel_id)
                    DO UPDATE SET reason = EXCLUDED.reason
                ''', guild_id, member_object.id, reason, static_channel_id)
            if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                await member_object.edit(mute=True)
                mute_type = "server muted" if is_owner else "muted"
                await self.handler.send_message(ctx, content=f'{member_object.mention} has been {mute_type} in <#{static_channel_id}> with reason {reason}.')
            else:
                mute_type = "server muted" if is_owner else "muted"
                await self.handler.send_message(ctx, content=f'{member_object.mention} has been {mute_type} in <#{static_channel_id}> with reason {reason}.')
        return mute_alias
    
    @commands.hybrid_command(
        name='mutes',
        help='Lists all users currently muted in a specific voice channel.'
    )
    @commands.check(is_owner_developer_coordinator)
    async def list_mutes(
        ctx,
        channel: str = commands.parameter(description='Voice channel mention or ID.')
    ) -> None:
        guild_id = ctx.guild.id
        try:
            if isinstance(channel, discord.VoiceChannel):
                channel_id = channel.id
            elif channel.isdigit():
                channel_id = int(channel)
            elif channel.startswith('<#') and channel.endswith('>'):
                channel_id = int(channel.strip('<#>'))
            else:
                raise ValueError
        except ValueError:
            return await self.handler.send_message(ctx, content='❌ Invalid channel mention or ID.')
        async with self.db_pool.acquire() as conn:
            records = await conn.fetch('''
                SELECT u.user_id, u.mute_channel_ids, mr.reason
                FROM users u
                JOIN mute_reasons mr ON u.user_id = mr.user_id
                WHERE $1 = ANY(u.mute_channel_ids) AND mr.channel_id = $1 AND mr.guild_id = $2
            ''', channel_id, guild_id)
        if not records:
            return await self.handler.send_message(ctx, content=f'ℹ️ No users are currently muted in <#{channel_id}>.')
        description_lines = []
        for record in records:
            user = ctx.guild.get_member(record['user_id'])
            mention = user.mention if user else f'`{record["user_id"]}`'
            reason = record['reason']
            description_lines.append(f'• {mention} — {reason}')
        embed = discord.Embed(
            title=f'Muted Users in <#{channel_id}>',
            description='\n'.join(description_lines),
            color=discord.Color.orange()
        )
        await self.handler.send_message(ctx, embed=embed)

    def create_unban_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(
            name=command_name,
            help='Unban a user from a voice channel.'
        )
        @commands.check(is_owner_developer_coordinator_moderator)
        async def unban_alias(
            ctx,
            member: str = commands.parameter(description='Mention or user ID of the member to unban.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Optional reason for unbanning.')
        ) -> None:
            guild_id = ctx.guild.id
            member_id = None
            if member.isdigit():
                member_id = int(member)
            elif member.startswith('<@') and member.endswith('>'):
                try:
                    member_id = int(member.strip('<@!>'))
                except ValueError:
                    pass
            member_object = ctx.guild.get_member(member_id) if member_id else None
            if not member_object:
                return await self.handler.send_message(ctx, content='❌ Could not resolve a valid member.')
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('ban', {}).get(command_name)
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'❌ No channel alias mapping found for `{command_name}`.')
            role_id = None
            ban_channel_id = self.bot.command_aliases.get(guild_id, {}).get('ban', {}).get(command_name)
            role_aliases = self.bot.command_aliases.get(guild_id, {}).get('role', {})
            for alias_data in role_aliases.values():
                if alias_data.get('channel_id') == ban_channel_id:
                    role_id = alias_data.get('role_id')
                    break
            if not role_id:
                return await self.handler.send_message(ctx, content=f'❌ No ban role alias mapping found for `{command_name}`.')
            role = ctx.guild.get_role(role_id)
            if not role:
                return await self.handler.send_message(ctx, content=f'❌ Could not resolve role ID `{role_id}`.')
            try:
                await member_object.remove_roles(role, reason=f"Unbanned from <#{static_channel_id}>: {reason or 'No reason provided'}")
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content='❌ Missing permissions to remove ban role.')
            async with self.db_pool.acquire() as conn:
                await conn.execute('''
                    DELETE FROM active_bans
                    WHERE user_id = $1 AND channel_id = $2
                ''', member_object.id, static_channel_id)
                await conn.execute('''
                    DELETE FROM ban_expirations
                    WHERE user_id = $1 AND channel_id = $2
                ''', member_object.id, static_channel_id)
                await conn.execute('''
                    UPDATE users
                    SET ban_channel_ids = array_remove(ban_channel_ids, $2),
                        updated_at = NOW()
                    WHERE user_id = $1
                ''', member_object.id, static_channel_id)
            await self.handler.send_message(
                ctx,
                content=f'🔊 {member_object.mention} has been unbanned from <#{static_channel_id}>.'
            )
        return unban_alias
    
    def create_unmute_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(name=command_name, help='Unmutes a member in a specific VC.')
        @commands.check(is_owner_developer_coordinator_moderator)
        async def unmute_alias(
            ctx,
            member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Include a reason for the unmute.')
        ) -> None:
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member.isdigit():
                member_id = int(member)
            elif member.startswith('<@') and member.endswith('>'):
                try:
                    member_id = int(member.strip('<@!>'))
                except ValueError:
                    pass
            if member_id:
                member_object = ctx.guild.get_member(member_id)
            if not member_object:
                return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('unmute', {}).get(command_name)
            async with self.db_pool.acquire() as conn:
                row = await conn.fetchrow('''
                    SELECT source FROM active_mutes
                    WHERE user_id = $1 AND channel_id = $2
                ''', member_object.id, static_channel_id)
                if not row or row["source"] != "bot":
                    return await self.handler.send_message(ctx, content=f"❌ {member_object.mention} was not muted by the bot in <#{static_channel_id}>.")
                if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                    await member_object.edit(mute=False)
                await conn.execute('''
                    DELETE FROM active_mutes
                    WHERE user_id = $1 AND channel_id = $2 AND source = 'bot'
                ''', member_object.id, static_channel_id)
                await conn.execute('''
                    UPDATE users
                    SET mute_channel_ids = array_remove(mute_channel_ids, $2),
                        updated_at = NOW()
                    WHERE user_id = $1
                ''', member_object.id, static_channel_id)
                await conn.execute('''
                    INSERT INTO mute_reasons (guild_id, user_id, channel_id, reason)
                    VALUES ($1, $2, $3, $4)
                    ON CONFLICT (guild_id, user_id, channel_id)
                    DO UPDATE SET reason = EXCLUDED.reason
                ''', guild_id, member_object.id, static_channel_id, reason)
            if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                await member_object.edit(mute=False)
                await self.handler.send_message(ctx, content=f'{member_object.mention} has been unmuted in <#{static_channel_id}>.')
            else:
                await self.handler.send_message(ctx, content=f'{member_object.mention} is no longer marked as muted in <#{static_channel_id}>.')
        return unmute_alias
        
    @commands.hybrid_command(name='xalias', help='Deletes an alias.')
    @commands.check(is_owner_developer_coordinator)
    async def delete_alias(
        self,
        ctx,
        alias_name: str = commands.parameter(description='Includ an alias name')
    ) -> None:
        if not alias_name.strip():
            await self.handler.send_message(ctx, '❌ `alias_name` cannot be empty.', ephemeral=True)
            return
        guild_id = ctx.guild.id
        alias_type = None
        for candidate in ('mute', 'unmute', 'ban', 'unban', 'flag'):
            if alias_name in self.bot.command_aliases.get(guild_id, {}).get(candidate, {}):
                alias_type = candidate
                break
        if alias_type.lower() not in {'mute', 'unmute', 'ban', 'unban', 'flag'}:
            await self.handler.send_message(ctx, '❌ `alias_type` must be either `mute`, `unmute`, `ban`, `unban`, or `flag`.', ephemeral=True)
            return
        if not alias_type:
            await self.handler.send_message(ctx, f'❌ Alias `{alias_name}` not found.', ephemeral=True)
            return
        alias_map = self.bot.command_aliases.get(guild_id, {}).get(alias_type.lower(), {})
        if alias_name not in alias_map:
            await self.handler.send_message(ctx, f'❌ Alias `{alias_name}` not found in `{alias_type}` for guild `{guild_id}`.', ephemeral=True)
            return
        async with self.db_pool.acquire() as conn:
            await conn.execute(
                'DELETE FROM command_aliases WHERE guild_id = $1 AND alias_type = $2 AND alias_name = $3',
                guild_id, alias_type.lower(), alias_name
            )
        self.bot.command_aliases[guild_id][alias_type.lower()].pop(alias_name, None)
        await self.handler.send_message(ctx, content=f'✅ Deleted alias `{alias_name}` from `{alias_type}`.')
        
    @commands.hybrid_command(name='aliases', help='List all the aliases in the current guild, or filter by channel.')
    @commands.check(is_owner_developer_coordinator_moderator)
    async def list_aliases(self, ctx, channel: discord.TextChannel = None) -> None:
        guild_id = ctx.guild.id
        aliases = self.bot.command_aliases.get(guild_id, {})
        if not aliases:
            await self.handler.send_message(ctx, content='No aliases defined in this guild.')
            return
        pages = []
        if channel:
            channel_id = channel.id
            found_aliases = False
            for kind in ('mute', 'unmute', 'ban', 'unban', 'flag'):
                entries = aliases.get(kind, {})
                if not entries:
                    continue
                channel_entries = {name: cid for name, cid in entries.items() if cid == channel_id}
                if channel_entries:
                    found_aliases = True
                    embed = discord.Embed(
                        title=f'{kind.capitalize()} Aliases for {channel.mention}',
                        description='\n'.join(f'`{name}` → <#{cid}>' for name, cid in channel_entries.items()),
                        color=discord.Color.blue()
                    )
                    pages.append(embed)
            if not found_aliases:
                await self.handler.send_message(
                    ctx,
                    content=f'No aliases found for {channel.mention}.'
                )
                return
        else:
            for kind in ('mute', 'unmute', 'ban', 'unban', 'flag'):
                entries = aliases.get(kind, {})
                if not entries:
                    continue
                embed = discord.Embed(
                    title=f'{kind.capitalize()} Aliases for {ctx.guild.name}',
                    description='\n'.join(f'`{name}` → <#{cid}>' for name, cid in entries.items()),
                    color=discord.Color.blue()
                )
                pages.append(embed)
        if not pages:
            await self.handler.send_message(ctx, content='No aliases found.')
            return
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()
    
    def create_unmute_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(name=command_name, help='Unmutes a member in a specific VC.')
        @commands.check(is_owner_developer_coordinator_moderator)
        async def unmute_alias(
            ctx,
            member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Include a reason for the unmute.')
        ) -> None:
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member.isdigit():
                member_id = int(member)
            elif member.startswith('<@') and member.endswith('>'):
                try:
                    member_id = int(member.strip('<@!>'))
                except ValueError:
                    pass
            if member_id:
                member_object = ctx.guild.get_member(member_id)
            if not member_object:
                return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('unmute', {}).get(command_name)
            async with self.db_pool.acquire() as conn:
                row = await conn.fetchrow('''
                    SELECT source FROM active_mutes
                    WHERE user_id = $1 AND channel_id = $2
                ''', member_object.id, static_channel_id)
                if not row or row["source"] != "bot":
                    return await self.handler.send_message(ctx, content=f"❌ {member_object.mention} was not muted by the bot in <#{static_channel_id}>.")
                if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                    await member_object.edit(mute=False)
                await conn.execute('''
                    DELETE FROM active_mutes
                    WHERE user_id = $1 AND channel_id = $2 AND source = 'bot'
                ''', member_object.id, static_channel_id)
                await conn.execute('''
                    UPDATE users
                    SET mute_channel_ids = array_remove(mute_channel_ids, $2),
                        updated_at = NOW()
                    WHERE user_id = $1
                ''', member_object.id, static_channel_id)
                await conn.execute('''
                    INSERT INTO mute_reasons (guild_id, user_id, channel_id, reason)
                    VALUES ($1, $2, $3, $4)
                    ON CONFLICT (guild_id, user_id, channel_id)
                    DO UPDATE SET reason = EXCLUDED.reason
                ''', guild_id, member_object.id, static_channel_id, reason)
            if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                await member_object.edit(mute=False)
                await self.handler.send_message(ctx, content=f'{member_object.mention} has been unmuted in <#{static_channel_id}>.')
            else:
                await self.handler.send_message(ctx, content=f'{member_object.mention} is no longer marked as muted in <#{static_channel_id}>.')
        return unmute_command
    

    @commands.hybrid_command(
        name='role',
        help='Associate a voice channel with a role used for bans.'
    )
    @commands.check(is_owner_developer_coordinator)
    async def role(
        self,
        ctx: commands.Context,
        channel: discord.VoiceChannel,
        role: discord.Role
    ):
        guild_id = ctx.guild.id
        async with self.db_pool.acquire() as conn:
            await conn.execute(
                '''
                INSERT INTO channel_roles (guild_id, channel_id, role_id)
                VALUES ($1, $2, $3)
                ON CONFLICT (guild_id, channel_id) DO UPDATE SET role_id = EXCLUDED.role_id
                ''',
                guild_id, channel.id, role.id
            )
        await self.handler.send_message(ctx, f"✅ Associated voice channel {channel.mention} with role {role.mention}.")
        
    def create_flag_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(
            name=command_name,
            help='Flag a user in the database for the voice channel mapped to this alias.'
        )
        @commands.check(is_owner_developer_coordinator_moderator)
        async def flag_alias(
            ctx,
            user: str
        ) -> None:
            guild_id = ctx.guild.id
            flag_aliases = self.bot.command_aliases.get(guild_id, {}).get('flag', {})
            channel_id = flag_aliases.get(command_name)
            print("tezt")
            if not channel_id:
                return await self.handler.send_message(ctx, content=f'❌ No flag alias configured for `{command_name}`.')
            if not user:
                return await self.handler.send_message(ctx, content='❌ You must provide a user ID or mention.')
            if re.fullmatch(r'<@!?\d+>', user):
                user_id = int(re.sub(r'\D', '', user))
            elif user.isdigit():
                user_id = int(user)
            else:
                return await self.handler.send_message(ctx, content='❌ Invalid user ID or mention.')
            print("test")
            select_sql = '''
                SELECT 1 FROM users
                WHERE user_id = $1 AND $2 = ANY(flagged_channel_ids)
            '''
            insert_sql = '''
                INSERT INTO users (user_id, flagged_channel_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[])
                ON CONFLICT (user_id)
                DO UPDATE SET flagged_channel_ids = (
                    SELECT ARRAY(
                        SELECT DISTINCT unnest(users.flagged_channel_ids || EXCLUDED.flagged_channel_ids)
                    )
                )
            '''
            try:
                async with self.db_pool.acquire() as conn:
                    already_flagged = await conn.fetchval(select_sql, user_id, channel_id)
                    if already_flagged:
                        return await self.handler.send_message(ctx, content=f'ℹ️ <@{user_id}> is already flagged for <#{channel_id}>.')
                    await conn.execute(insert_sql, user_id, channel_id)
                    await self.handler.send_message(ctx, content=f'✅ Flagged <@{user_id}> for channel <#{channel_id}>.')
            except Exception as e:
                await self.handler.send_message(ctx, content=f'❌ Database error: {e}')
                raise

    #
    #  Help Command: Provides a scope-limited command to investigate available bot commands.
    #
    @commands.hybrid_command(name='help')
    async def help(
        self,
        ctx,
        *,
        command_name: str = commands.parameter(default=None, description='Include a command name')
    ) -> None:
        bot = ctx.bot
        if command_name:
            cmd = bot.get_command(command_name.lower())
            if not cmd:
                await self.handler.send_message(ctx, f'❌ Command `{command_name}` not found.')
                return
            if cmd.hidden:
                await self.handler.send_message(ctx, f'❌ Command `{command_name}` is hidden.')
                return
            if not await cmd.can_run(ctx):
                await self.handler.send_message(ctx, f'❌ You do not have permission to run `{command_name}`.')
                return
            embed = discord.Embed(
                title=f'/{cmd.name}',
                description=cmd.help or 'No description provided.',
                color=discord.Color.blue()
            )
            sig = inspect.signature(cmd.callback)
            parameters = list(sig.parameters.items())[2:]
            if parameters:
                usage_parts = [f"/{cmd.name}"]
                param_details = []
                for name, param in parameters:
                    is_optional = param.kind == inspect.Parameter.KEYWORD_ONLY
                    default = param.default
                    description = (
                        getattr(default, 'description', None)
                        if isinstance(default, commands.Parameter)
                        else None
                    )
                    annotation = (
                        param.annotation.__name__
                        if hasattr(param.annotation, '__name__')
                        else str(param.annotation)
                    )
                    if is_optional:
                        usage_parts.append(f"[{name}]")
                    else:
                        usage_parts.append(f"<{name}>")
                    detail = f"**{name}** ({annotation})"
                    if description:
                        detail += f": {description}"
                    param_details.append(detail)
                embed.add_field(
                    name="Usage",
                    value=f"`{' '.join(usage_parts)}`",
                    inline=False
                )
                if param_details:
                    embed.add_field(
                        name="Parameter Details",
                        value="\n".join(param_details),
                        inline=False
                    )
            await self.handler.send_message(ctx, embed=embed)
            print("text")
            return
        all_commands = await self.get_available_commands(bot, ctx)
        if not all_commands:
            await self.handler.send_message(ctx, '❌ No commands available to you.')
            return
        cog_map: dict[str, list[commands.Command]] = {}
        for command in all_commands:
            if command.hidden:
                continue
            cog_map.setdefault(command.cog_name or 'Uncategorized', []).append(command)
        pages = []
        for cog_name in sorted(cog_map):
            commands_in_cog = sorted(cog_map[cog_name], key=lambda c: c.name)
            embed = discord.Embed(title=f'{cog_name} Commands', color=discord.Color.green())
            embed.description = '\n'.join(
                f'**/{cmd.name}** – {cmd.help or "No description"}' for cmd in commands_in_cog
            )
            pages.append(embed)
        paginator = Paginator(bot, ctx, pages, self.handler)
        await paginator.start()
    
    #
    #  Reason Command: Fetches the reason for an outstanding mute or unmute triggered by the bot.
    #                  Empty if manually muted.
    #
    @commands.hybrid_command(name='reason', help='Get the reason for a mute, unmute, ban, or unban.')
    @commands.check(is_owner_developer_coordinator_moderator)
    async def get_mute_unmute_reason(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.')
    ) -> None:
        guild_id = ctx.guild.id
        member_id = None
        member_object = None
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='❌ Could not resolve a valid guild member from your input.')
        async with self.db_pool.acquire() as conn:
            mute_rows = await conn.fetch('''
                SELECT channel_id, reason, action
                FROM mute_reasons
                WHERE guild_id = $1 AND user_id = $2
            ''', guild_id, member_object.id)
            ban_rows = await conn.fetch('''
                SELECT channel_id, reason, action
                FROM ban_reasons
                WHERE guild_id = $1 AND user_id = $2
            ''', guild_id, member_object.id)
    
        if not mute_rows and not ban_rows:
            return await self.handler.send_message(ctx, content=f'ℹ️ No mute, unmute, ban, or unban history found for {member_object.mention}.')
        lines = []
        def format_row(row, kind):
            channel_id = row['channel_id']
            reason = row['reason'] or '*No reason provided*'
            action = row.get('action', kind)
            return f'• [{action}] <#{channel_id}>: `{reason}`'
        lines += [format_row(row, 'mute') for row in mute_rows]
        lines += [format_row(row, 'ban') for row in ban_rows]
        content = f'📄 Disciplinary reasons for {member_object.mention}:\n' + '\n'.join(lines)
        await self.handler.send_message(ctx, content=content)
        
    @commands.hybrid_command(name='help')
    async def help(
        self,
        ctx,
        *,
        command_name: str = commands.parameter(default=None, description='Include a command name')
    ) -> None:
        bot = ctx.bot
        if command_name:
            cmd = bot.get_command(command_name.lower())
            if not cmd:
                await self.handler.send_message(ctx, f'❌ Command `{command_name}` not found.')
                return
            if cmd.hidden:
                await self.handler.send_message(ctx, f'❌ Command `{command_name}` is hidden.')
                return
            if not await cmd.can_run(ctx):
                await self.handler.send_message(ctx, f'❌ You do not have permission to run `{command_name}`.')
                return
            embed = discord.Embed(
                title=f'/{cmd.name}',
                description=cmd.help or 'No description provided.',
                color=discord.Color.blue()
            )
            sig = inspect.signature(cmd.callback)
            parameters = list(sig.parameters.items())[2:]
            if parameters:
                usage_parts = [f"/{cmd.name}"]
                param_details = []
                for name, param in parameters:
                    is_optional = param.kind == inspect.Parameter.KEYWORD_ONLY
                    default = param.default
                    description = (
                        getattr(default, 'description', None)
                        if isinstance(default, commands.Parameter)
                        else None
                    )
                    annotation = (
                        param.annotation.__name__
                        if hasattr(param.annotation, '__name__')
                        else str(param.annotation)
                    )
                    if is_optional:
                        usage_parts.append(f"[{name}]")
                    else:
                        usage_parts.append(f"<{name}>")
                    detail = f"**{name}** ({annotation})"
                    if description:
                        detail += f": {description}"
                    param_details.append(detail)
                embed.add_field(
                    name="Usage",
                    value=f"`{' '.join(usage_parts)}`",
                    inline=False
                )
                if param_details:
                    embed.add_field(
                        name="Parameter Details",
                        value="\n".join(param_details),
                        inline=False
                    )
            await self.handler.send_message(ctx, embed=embed)
            return
        all_commands = await self.get_available_commands(bot, ctx)
        if not all_commands:
            await self.handler.send_message(ctx, '❌ No commands available to you.')
            return
        permission_groups = await self.group_commands_by_permission(bot, ctx, all_commands)
        pages = []
        permission_order = [
            ('Owner', '`Owner` inherits `developer`'),
            ('Developer', '`Developer` inherits `coordinator`.'),
            ('Coordinator', '`Coordinator inherits `moderator`.'),
            ('Moderator', 'Moderators can use these commands.'),
            ('Everyone', 'No commands available.')
        ]
        for perm_level, description in permission_order:
            if perm_level not in permission_groups:
                continue
            commands_in_level = sorted(permission_groups[perm_level], key=lambda c: c.name)
            if not commands_in_level:
                continue
            embed = discord.Embed(
                title=f'{perm_level} Commands',
                description=description,
                color=self.get_permission_color(perm_level)
            )
            cog_map = {}
            for command in commands_in_level:
                if command.hidden:
                    continue
                cog_name = command.cog_name or 'Uncategorized'
                cog_map.setdefault(cog_name, []).append(command)
            for cog_name in sorted(cog_map):
                commands_in_cog = sorted(cog_map[cog_name], key=lambda c: c.name)
                command_list = '\n'.join(
                    f'**/{cmd.name}** – {cmd.help or "No description"}'
                    for cmd in commands_in_cog
                )
                if len(command_list) > 1024:
                    chunks = self.split_command_list(commands_in_cog)
                    for i, chunk in enumerate(chunks):
                        field_name = f'{cog_name}' if i == 0 else f'{cog_name} (cont.)'
                        embed.add_field(
                            name=field_name,
                            value=chunk,
                            inline=False
                        )
                else:
                    embed.add_field(
                        name=cog_name,
                        value=command_list,
                        inline=False
                    )
            pages.append(embed)
        if not pages:
            await self.handler.send_message(ctx, '❌ No commands available to you.')
            return
        paginator = Paginator(bot, ctx, pages)
        await paginator.start()
    
    async def group_commands_by_permission(self, bot, ctx, commands_list):
        permission_groups = {
            'Owner': [],
            'Developer': [],
            'Coordinator': [],
            'Moderator': [],
            'Everyone': []
        }
        for command in commands_list:
            perm_level = await self.get_command_permission_level(bot, ctx, command)
            if perm_level in permission_groups:
                permission_groups[perm_level].append(command)
        return permission_groups
    
    async def get_command_permission_level(self, bot, ctx, command):
        if not hasattr(command, 'checks') or not command.checks:
            return 'Everyone'
        permission_levels = {
            'Owner': 4,
            'Developer': 3,
            'Coordinator': 2,
            'Moderator': 1,
            'Everyone': 0
        }
        highest_level = 'Everyone'
        highest_value = 0
        for check in command.checks:
            check_names = []
            if hasattr(check, '__name__'):
                check_names.append(check.__name__)
            if hasattr(check, '__wrapped__'):
                wrapped = check.__wrapped__
                if hasattr(wrapped, '__name__'):
                    check_names.append(wrapped.__name__)
            if hasattr(check, 'predicate') and hasattr(check.predicate, '__name__'):
                check_names.append(check.predicate.__name__)
            check_str = str(check)
            if 'function' in check_str:
                try:
                    func_name = check_str.split('function ')[1].split(' ')[0]
                    check_names.append(func_name)
                except:
                    pass
            for check_name in check_names:
                level = None
                if check_name == 'is_owner_developer_coordinator_moderator':
                    level = 'Moderator'
                elif check_name == 'is_owner_developer_coordinator':
                    level = 'Coordinator'
                elif check_name == 'is_owner_developer':
                    level = 'Developer'
                elif check_name in ['is_owner', 'is_guild_owner', 'is_system_owner']:
                    level = 'Owner'
                elif check_name == 'is_developer':
                    level = 'Developer'
                elif check_name == 'is_coordinator':
                    level = 'Coordinator'
                elif check_name in ['is_moderator', 'is_channel_moderator']:
                    level = 'Moderator'
                if level and permission_levels[level] > highest_value:
                    highest_level = level
                    highest_value = permission_levels[level]
        return highest_level
    
    def get_permission_color(self, perm_level):
        colors = {
            'Owner': discord.Color.red(),
            'Developer': discord.Color.purple(),
            'Coordinator': discord.Color.orange(),
            'Moderator': discord.Color.blue(),
            'Everyone': discord.Color.green()
        }
        return colors.get(perm_level, discord.Color.greyple())
    
    def split_command_list(self, commands_list, max_length=1024):
        chunks = []
        current_chunk = []
        current_length = 0
        for cmd in commands_list:
            cmd_line = f'**/{cmd.name}** – {cmd.help or "No description"}\n'
            cmd_length = len(cmd_line)
            if current_length + cmd_length > max_length and current_chunk:
                chunks.append('\n'.join(current_chunk))
                current_chunk = [cmd_line.rstrip()]
                current_length = cmd_length
            else:
                current_chunk.append(cmd_line.rstrip())
                current_length += cmd_length
        if current_chunk:
            chunks.append('\n'.join(current_chunk))
        return chunks
    
async def setup(bot: commands.Bot):
    cog = Hybrid(bot)
    await bot.add_cog(cog)
