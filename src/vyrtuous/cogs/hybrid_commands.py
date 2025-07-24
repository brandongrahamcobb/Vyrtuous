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
from collections import defaultdict
from discord import app_commands
from discord.ext import commands
from discord.ext.commands import Command
from vyrtuous.service.discord_message_service import DiscordMessageService, Paginator
from vyrtuous.service.check_service import *
from vyrtuous.inc.helpers import *
from types import MethodType
from typing import List

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
                    self.bot.add_command(cmd)
                elif alias_type == 'unmute':
                    cmd = self.create_unmute_alias(alias_name)
                    self.bot.add_command(cmd)
                elif alias_type == 'ban':
                    cmd = self.create_ban_alias(alias_name)
                    self.bot.add_command(cmd)
                elif alias_type == 'unban':
                    cmd = self.create_unban_alias(alias_name)
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
                print(f"❌ Exception while checking command '{command}': {e}")
        return available_commands
    
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
        member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ) -> None:
        guild_id = ctx.guild.id
        member_id = None
        member_object = None
        if member_input.isdigit():
            member_id = int(member_input)
        elif member_input.startswith('<@') and member_input.endswith('>'):
            try:
                member_id = int(member_input.strip('<@!>'))
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
        member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ) -> None:
        member_id = None
        member_object = None
        if member_input.isdigit():
            member_id = int(member_input)
        elif member_input.startswith('<@') and member_input.endswith('>'):
            try:
                member_id = int(member_input.strip('<@!>'))
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
        member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel_input: str = commands.parameter(description='Tag a channel or include its snowflake ID.')
    ) -> None:
        member_id = None
        member_object = None
        if member_input.isdigit():
            member_id = int(member_input)
        elif member_input.startswith('<@') and member_input.endswith('>'):
            try:
                member_id = int(member_input.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        resolved_channel = None
        if channel_input.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel_input))
        elif channel_input.startswith('<#') and channel_input.endswith('>'):
            try:
                channel_id = int(channel_input.strip('<#>'))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel_input.lower():
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
        member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel_input: str = commands.parameter(description='Tag a VC or include its snowflake ID.')
    ) -> None:
        member_id = None
        member_object = None
        if member_input.isdigit():
            member_id = int(member_input)
        elif member_input.startswith('<@') and member_input.endswith('>'):
            try:
                member_id = int(member_input.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        resolved_channel = None
        if channel_input.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel_input))
        elif channel_input.startswith('<#') and channel_input.endswith('>'):
            try:
                channel_id = int(channel_input.strip('<#>'))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel_input.lower():
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

    @commands.hybrid_command(name='mods', help='Lists VC moderators.')
    @commands.check(is_owner_developer_coordinator_moderator)
    async def list_moderators(self, ctx) -> None:
        guild = ctx.guild
        pages = []
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT user_id, moderator_ids
                FROM users
                WHERE cardinality(moderator_ids) > 0
            ''')
        for row in rows:
            user_id = row['user_id']
            moderator_ids = row['moderator_ids'] or []
            valid_channels = [
                guild.get_channel(cid)
                for cid in moderator_ids
                if (guild.get_channel(cid) and isinstance(guild.get_channel(cid), discord.VoiceChannel))
            ]
            if not valid_channels:
                continue
            user = guild.get_member(user_id)
            display_name = user.display_name if user else f'User ID {user_id}'
            embed = discord.Embed(
                title=f'VC Mod: {display_name}',
                description='\n'.join(f'<#{channel.id}> — {channel.name}' for channel in valid_channels),
                color=discord.Color.blue()
            )
            embed.set_footer(text=f'User ID: {user_id}')
            pages.append(embed)
        if not pages:
            await self.handler.send_message(ctx, content='No VC moderators are configured in this guild.')
            return
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()
        
    #
    # Alias Commands: creation
    #                 deletion
    #                 listing
    #
    @commands.hybrid_command(
        name='alias',
        help='Set an alias for a mute, unmute, ban, unban, or role action.'
    )
    @commands.check(is_owner_developer_coordinator)  # Devs can bypass all checks
    async def create_alias(
        self,
        ctx,
        alias_type: str = commands.parameter(description='One of: `mute`, `unmute`, `ban`, `role`, `xrole`, `unban`'),
        alias_name: str = commands.parameter(description='Alias name'),
        target: str = commands.parameter(description='VC/Role ID, mention, or name')
    ) -> None:
        alias_type = alias_type.lower()
        guild_id = ctx.guild.id
        valid_types = {'mute', 'unmute', 'ban', 'role', 'xrole', 'unban'}
        if alias_type not in valid_types:
            return await ctx.send(f'❌ Invalid alias type. Must be one of: {", ".join(valid_types)}', ephemeral=True)
        if not alias_name.strip():
            return await ctx.send('❌ Alias name cannot be empty.', ephemeral=True)
        resolved_id = None
        resolved_type = None
        if target.isdigit():
            resolved_id = int(target)
        elif target.startswith('<#') and target.endswith('>'):
            try:
                resolved_id = int(target.strip('<#>'))
            except ValueError:
                pass
        elif target.startswith('<@&') and target.endswith('>'):
            try:
                resolved_id = int(target.strip('<@&>'))
            except ValueError:
                pass
        if alias_type in {'mute', 'unmute', 'ban', 'unban'}:
            if not resolved_id:
                for vc in ctx.guild.voice_channels:
                    if vc.name.lower() == target.lower():
                        resolved_id = vc.id
                        break
            channel = ctx.guild.get_channel(resolved_id)
            if not resolved_id or not channel or channel.type != discord.ChannelType.voice:
                return await ctx.send('❌ Could not resolve a valid voice channel.', ephemeral=True)
            resolved_type = 'channel_id'
        else:
            if not resolved_id:
                for role in ctx.guild.roles:
                    if role.name.lower() == target.lower():
                        resolved_id = role.id
                        break
            role = ctx.guild.get_role(resolved_id)
            if not resolved_id or not role:
                return await ctx.send('❌ Could not resolve a valid role.', ephemeral=True)
            resolved_type = 'role_id'
            if not await is_owner_or_coordinator(ctx):
                return await ctx.send('❌ Only developers can manage role-related aliases.', ephemeral=True)
        self.bot.command_aliases.setdefault(guild_id, {}).setdefault(alias_type, {})[alias_name] = resolved_id
        async with self.db_pool.acquire() as conn:
            await conn.execute(
                '''
                INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_id, alias_type, alias_name)
                DO UPDATE SET channel_id = EXCLUDED.channel_id
                ''',
                guild_id, alias_type, alias_name, resolved_id
            )
        target_object = channel if resolved_type == 'channel_id' else role
        await self.handler.send_message(
            ctx,
            content=f'✅ Alias `{alias_name}` ({alias_type}) set to {target_object.mention}.'
        )
        if alias_type == 'mute':
            cmd = self.create_mute_alias(alias_name)
            self.bot.add_command(cmd)
        elif alias_type == 'unmute':
            cmd = self.create_unmute_alias(alias_name)
            self.bot.add_command(cmd)
        elif alias_type == 'ban':
            cmd = self.create_ban_alias(alias_name)
            self.bot.add_command(cmd)
        elif alias_type == 'unban':
            cmd = self.create_unban_alias(alias_name)
            self.bot.add_command(cmd)
    
    def create_ban_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(
            name=command_name,
            help='Ban a user from a voice channel. Use 0 hours for permanent ban (coordinators only).'
        )
        @commands.check(is_owner_developer_coordinator_moderator)
        async def ban_alias(
            ctx,
            member_input: str = commands.parameter(description='Mention or user ID of the member to ban.'),
            duration_hours: Optional[int] = commands.parameter(default=6, description='Duration of ban in hours (0 = permanent).'),
            *,
            reason: str = commands.parameter(default='', description='Reason for ban (required for permanent).')
        ) -> None:
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member_input.isdigit():
                member_id = int(member_input)
            elif member_input.startswith('<@') and member_input.endswith('>'):
                try:
                    member_id = int(member_input.strip('<@!>'))
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
            role_id = self.bot.command_aliases.get(guild_id, {}).get('banrole', {}).get(command_name)
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'❌ No channel alias mapping found for `{command_name}`.')
            if not role_id:
                return await self.handler.send_message(ctx, content=f'❌ No ban role alias mapping found for `{command_name}`.')
            role = ctx.guild.get_role(role_id)
            if not role:
                return await self.handler.send_message(ctx, content=f'❌ Could not resolve role ID `{role_id}`.')
            try:
                await member_object.add_roles(role, reason=f"Banned from <#{static_channel_id}>: {reason or 'No reason provided'}")
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content='❌ Missing permissions to assign ban role.')
            async with self.db_pool.acquire() as conn:
                await conn.execute('''
                    INSERT INTO active_bans (user_id, channel_id, expires_at)
                    VALUES ($1, $2, $3)
                    ON CONFLICT (user_id, channel_id) DO UPDATE SET expires_at = EXCLUDED.expires_at
                ''', member_object.id, static_channel_id,
                    None if duration_hours == 0 else discord.utils.utcnow() + timedelta(hours=duration_hours)
                )
                await conn.execute('''
                    UPDATE users
                    SET ban_channel_ids = array_append(ban_channel_ids, $2),
                        updated_at = NOW()
                    WHERE user_id = $1
                ''', member_object.id, static_channel_id)
    
                await conn.execute('''
                    INSERT INTO ban_reasons (guild_id, user_id, channel_id, reason)
                    VALUES ($1, $2, $3, $4)
                    ON CONFLICT (guild_id, user_id, channel_id)
                    DO UPDATE SET reason = EXCLUDED.reason
                ''', guild_id, member_object.id, static_channel_id, reason or 'No reason provided')
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
                        user.id, channel.id, expires_at
                    )
        return ban_alias
    
    def create_mute_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(name=command_name, help='Mutes a member in a specific VC.')
        @commands.check(is_owner_developer_coordinator_moderator)
        async def mute_alias(
            ctx,
            member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Optionally include a reason for the mute.')
        ) -> None:
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member_input.isdigit():
                member_id = int(member_input)
            elif member_input.startswith('<@') and member_input.endswith('>'):
                try:
                    member_id = int(member_input.strip('<@!>'))
                except ValueError:
                    pass
            if member_id:
                member_object = ctx.guild.get_member(member_id)
            if not member_object:
                return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('mute', {}).get(command_name)
            async with self.db_pool.acquire() as conn:
                await conn.execute('''
                    INSERT INTO active_mutes (user_id, channel_id, source)
                    VALUES ($1, $2, 'bot')
                    ON CONFLICT DO NOTHING
                ''', member_object.id, static_channel_id)
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
                await self.handler.send_message(ctx, content=f'{member_object.mention} has been muted in <#{static_channel_id}> with reason {reason}.')
            else:
                await self.handler.send_message(ctx, content=f'{member_object.mention} has been muted in <#{static_channel_id}> with reason {reason}.')
        return mute_alias
    
    def create_unban_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(
            name=command_name,
            help='Unban a user from a voice channel.'
        )
        @commands.check(is_owner_developer_coordinator_moderator)
        async def unban_alias(
            ctx,
            member_input: str = commands.parameter(description='Mention or user ID of the member to unban.'),
            *,
            reason: str = commands.parameter(default='', description='Optional reason for unbanning.')
        ) -> None:
            guild_id = ctx.guild.id
            member_id = None
            if member_input.isdigit():
                member_id = int(member_input)
            elif member_input.startswith('<@') and member_input.endswith('>'):
                try:
                    member_id = int(member_input.strip('<@!>'))
                except ValueError:
                    pass
            member_object = ctx.guild.get_member(member_id) if member_id else None
            if not member_object:
                return await self.handler.send_message(ctx, content='❌ Could not resolve a valid member.')
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('ban', {}).get(command_name)
            role_id = self.bot.command_aliases.get(guild_id, {}).get('banrole', {}).get(command_name)
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'❌ No channel alias mapping found for `{command_name}`.')
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
            member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Include a reason for the unmute.')
        ) -> None:
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member_input.isdigit():
                member_id = int(member_input)
            elif member_input.startswith('<@') and member_input.endswith('>'):
                try:
                    member_id = int(member_input.strip('<@!>'))
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
            await ctx.send('❌ `alias_name` cannot be empty.', ephemeral=True)
            return
        guild_id = ctx.guild.id
        alias_type = None
        for candidate in ('mute', 'unmute', 'ban', 'unban', 'role', 'xrole'):
            if alias_name in self.bot.command_aliases.get(guild_id, {}).get(candidate, {}):
                alias_type = candidate
                break
        if alias_type.lower() not in {'mute', 'unmute', 'ban', 'unban', 'role', 'xrole'}:
            await ctx.send('❌ `alias_type` must be either `mute` or `unmute`.', ephemeral=True)
            return
        if not alias_type:
            await ctx.send(f'❌ Alias `{alias_name}` not found.', ephemeral=True)
            return
        alias_map = self.bot.command_aliases.get(guild_id, {}).get(alias_type.lower(), {})
        if alias_name not in alias_map:
            await ctx.send(f'❌ Alias `{alias_name}` not found in `{alias_type}` for guild `{guild_id}`.', ephemeral=True)
            return
        async with self.db_pool.acquire() as conn:
            await conn.execute(
                'DELETE FROM command_aliases WHERE guild_id = $1 AND alias_type = $2 AND alias_name = $3',
                guild_id, alias_type.lower(), alias_name
            )
        self.bot.command_aliases[guild_id][alias_type.lower()].pop(alias_name, None)
        await self.handler.send_message(ctx, content=f'✅ Deleted alias `{alias_name}` from `{alias_type}`.')
        
    @commands.hybrid_command(name='aliases', help='List all the aliases in the current guild.')
    @commands.check(is_owner_developer_coordinator_moderator)
    async def list_aliases(self, ctx) -> None:
        guild_id = ctx.guild.id
        aliases = self.bot.command_aliases.get(guild_id, {})
        pages = []
        for kind in ('mute', 'unmute', 'ban', 'unban', 'role', 'xrole'):
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
            await self.handler.send_message(ctx, content='No aliases defined in this guild.')
            return
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()
    
    def create_unmute_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(name=command_name, help='Unmutes a member in a specific VC.')
        @commands.check(is_owner_developer_coordinator_moderator)
        async def unmute_alias(
            ctx,
            member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Include a reason for the unmute.')
        ) -> None:
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member_input.isdigit():
                member_id = int(member_input)
            elif member_input.startswith('<@') and member_input.endswith('>'):
                try:
                    member_id = int(member_input.strip('<@!>'))
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
                await ctx.send(f'❌ Command `{command_name}` not found.')
                return
            if cmd.hidden:
                await ctx.send(f'❌ Command `{command_name}` is hidden.')
                return
            if not await cmd.can_run(ctx):
                await ctx.send(f'❌ You do not have permission to run `{command_name}`.')
                return
            embed = discord.Embed(
                title=f'/{cmd.name}',
                description=cmd.help or 'No description provided.',
                color=discord.Color.blue()
            )
            sig = inspect.signature(cmd.callback)
            parameters = list(sig.parameters.items())[2:]
            for name, param in parameters:
                is_optional = param.kind == inspect.Parameter.KEYWORD_ONLY
                default = param.default
                description = (
                    getattr(default, 'description', None)
                    if isinstance(default, commands.Parameter)
                    else None
                )
                default_value = (
                    getattr(default, 'default', None)
                    if isinstance(default, commands.Parameter)
                    else default
                )
                annotation = (
                    param.annotation.__name__
                    if hasattr(param.annotation, '__name__')
                    else str(param.annotation)
                )
                label = 'Optional' if is_optional else 'Required'
                display = f'Type: `{annotation}`\n{label}'
                if description:
                    display += f'\n{description}'
                embed.add_field(name=f'`{name}`', value=display, inline=False)
            await ctx.send(embed=embed)
            return
        all_commands = await self.get_available_commands(bot, ctx)
        if not all_commands:
            await ctx.send('❌ No commands available to you.')
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
        paginator = Paginator(bot, ctx, pages)
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
        member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.')
    ) -> None:
        guild_id = ctx.guild.id
        member_id = None
        member_object = None
        if member_input.isdigit():
            member_id = int(member_input)
        elif member_input.startswith('<@') and member_input.endswith('>'):
            try:
                member_id = int(member_input.strip('<@!>'))
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
    
async def setup(bot: commands.Bot):
    cog = Hybrid(bot)
    await bot.add_cog(cog)
