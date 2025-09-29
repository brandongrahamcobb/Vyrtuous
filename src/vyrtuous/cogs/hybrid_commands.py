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
from typing import Any, Optional, Union
import discord

from discord.ext.commands import Command

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.check_service import *
from vyrtuous.service.discord_message_service import DiscordMessageService, ChannelPaginator, Paginator
from vyrtuous.utils.setup_logging import logger

VEGAN_EMOJIS = [
    '\U0001F436',
    '\U0001F431',
    '\U0001F42D',
    '\U0001F439',
    '\U0001F430',
    '\U0001F98A',
    '\U0001F43B',
    '\U0001F43C',
    '\U0001F428',
    '\U0001F42F',
    '\U0001F981',
    '\U0001F42E',
    '\U0001F437',
    '\U0001F43D',
    '\U0001F438',
    '\U0001F435',
    '\U0001F412',
    '\U0001F98D',
    '\U0001F9A7',
    '\U0001F414',
    '\U0001F427',
    '\U0001F426',
    '\U0001F424',
    '\U0001F423',
    '\U0001F425',
    '\U0001F986',
    '\U0001F9A2',
    '\U0001F989',
    '\U0001F99A',
    '\U0001F99C',
    '\U0001F43A',
    '\U0001F99D',
    '\U0001F9A8',
    '\U0001F9A1',
    '\U0001F417',
    '\U0001F434',
    '\U0001F984',
    '\U0001F41D',
    '\U0001F41B',
    '\U0001F98B',
    '\U0001F40C',
    '\U0001F41E',
    '\U0001F40C',
    '\U0001FAB2',
    '\U0001F997',
    '\U0001F577',
    '\U0001F982',
    '\U0001F422',
    '\U0001F40D',
    '\U0001F98E',
    '\U0001F996',
    '\U0001F995',
    '\U0001F419',
    '\U0001F991',
    '\U0001F990',
    '\U0001F99E',
    '\U0001F980',
    '\U0001F421',
    '\U0001F420',
    '\U0001F41F',
    '\U0001F42C',
    '\U0001F988',
    '\U0001F433',
    '\U0001F40B',
    '\U0001F9AD',
    '\U0001F40A',
    '\U0001F406',
    '\U0001F405',
    '\U0001F403',
    '\U0001F402',
    '\U0001F42B',
    '\U0001F42A',
    '\U0001F999',
    '\U0001F992',
    '\U0001F98F',
    '\U0001F99B',
    '\U0001F418',
    '\U0001F998',
    '\U0001F9A5',
    '\U0001F9A6',
    '\U0001F9A8',
    '\U0001F9A9',
    '\U0001F54A'
]
PERMISSION_ORDER = ['Owner', 'Developer', 'Coordinator', 'Moderator', 'Everyone']

class Hybrid(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.bot.loop.create_task(self.load_server_muters())
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
        self.log_channels: dict[int, list[dict]] = {}
        self.server_muters: dict[int, set[int]] = defaultdict(set)
#        self.super = False

    async def cog_load(self) -> None:
        if not hasattr(self, '_loaded_aliases'):
            self._loaded_aliases = set()
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT guild_id, alias_type, alias_name, channel_id, role_id FROM command_aliases'
            )
            for row in rows:
                guild_id = row['guild_id']
                alias_type = row['alias_type']
                alias_name = row['alias_name']
                if alias_type in ('role', 'unrole'):
                    role_id = row['role_id']
                    channel_id = row['channel_id']
                    if role_id is None:
                        continue
                    self.bot.command_aliases \
                        .setdefault(guild_id, {}) \
                        .setdefault('role_aliases', {}) \
                        .setdefault(alias_type, {})[alias_name] = {
                            'role_id': int(role_id),
                            'channel_id': int(channel_id) if channel_id else None
                        }
                else:
                    channel_id = row['channel_id']
                    if channel_id is None:
                        continue
                    self.bot.command_aliases \
                        .setdefault(guild_id, {}) \
                        .setdefault('channel_aliases', {}) \
                        .setdefault(alias_type, {})[alias_name] = int(channel_id)
                if alias_name in self._loaded_aliases:
                    continue
                cmd = None
                if alias_type == 'mute':
                    cmd = self.create_voice_mute_alias(alias_name)
                elif alias_type == 'unmute':
                    cmd = self.create_unmute_alias(alias_name)
                elif alias_type == 'ban':
                    cmd = self.create_ban_alias(alias_name)
                elif alias_type == 'unban':
                    cmd = self.create_unban_alias(alias_name)
                elif alias_type == 'cow':
                    cmd = self.create_cow_alias(alias_name)
                elif alias_type == 'uncow':
                    cmd = self.create_uncow_alias(alias_name)
                elif alias_type == 'flag':
                    cmd = self.create_flag_alias(alias_name)
                elif alias_type == 'unflag':
                    cmd = self.create_unflag_alias(alias_name)
                elif alias_type == 'tmute':
                    cmd = self.create_text_mute_alias(alias_name)
                elif alias_type == 'untmute':
                    cmd = self.create_untextmute_alias(alias_name)
                elif alias_type == 'role':
                    cmd = self.create_role_alias(alias_name)
                elif alias_type == 'unrole':
                    cmd = self.create_unrole_alias(alias_name)
                if cmd and not self.bot.get_command(alias_name):
                    self.bot.add_command(cmd)
                    self._loaded_aliases.add(alias_name)
            await self.load_log_channels()
    
    async def load_log_channels(self):
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT guild_id, channel_id, type, snowflakes, enabled FROM log_channels;')
            log_channels: dict[int, list[dict]] = {}
            for r in rows:
                log_channels.setdefault(r['guild_id'], []).append({
                    "channel_id": r['channel_id'],
                    "type": r['type'] or "general",
                    "snowflakes": r['snowflakes'] or [],
                    "enabled": r['enabled']
                })
            self.log_channels = log_channels

#    @commands.command(name='toggle', hidden=True)
#    @is_owner_predicator()
#    async def toggle_feature(self, ctx: commands.Context):
#        self.super = not self.super
#        state = f'ON {self.get_random_emoji()}' if self.super else f'OFF \U0001F6AB'
#        await ctx.send(f'{self.get_random_emoji()} Feature switched {state}.')
#        if not self.super:
#            for channel in ctx.guild.channels:
#                if 'vegan' in channel.name.lower():
#                    await self.unrestrict(ctx.guild, channel)
#    
#    async def unrestrict(self, guild, vegan_channel):
#        async with self.bot.db_pool.acquire() as conn:
#            rows = await conn.fetch('''
#                SELECT discord_snowflake FROM users
#                WHERE $1 = ANY(coordinator_channel_ids) OR $1 = ANY(coordinator_ids)
#            ''', vegan_channel.id)
#        for uid in [row['discord_snowflake'] for row in rows]:
#            member = guild.get_member(uid)
#            if not member:
#                continue
#            async with self.bot.db_pool.acquire() as conn:
#                ban_rows = await conn.fetch('SELECT channel_id FROM active_bans WHERE discord_snowflake=$1', uid)
#                mute_rows = await conn.fetch('SELECT channel_id FROM active_voice_mutes WHERE discord_snowflake=$1', uid)
#                text_rows = await conn.fetch('SELECT channel_id FROM active_text_mutes WHERE discord_snowflake=$1', uid)
#            for r in ban_rows:
#                try: await guild.unban(discord.Object(id=uid), reason='Toggle OFF')
#                except: pass
#            for r in mute_rows:
#                ch = guild.get_channel(r['channel_id'])
#                if ch and member.voice and member.voice.mute: await member.edit(mute=False)
#            for r in text_rows:
#                ch = guild.get_channel(r['channel_id'])
#                text_mute_role = discord.utils.get(guild.roles, name='TextMuted')
#                if ch and text_mute_role and text_mute_role in member.roles: await member.remove_roles(text_mute_role)
#            async with self.bot.db_pool.acquire() as conn:
#                await conn.execute('DELETE FROM active_bans WHERE discord_snowflake=$1', uid)
#                await conn.execute('DELETE FROM active_voice_mutes WHERE discord_snowflake=$1', uid)
#                await conn.execute('DELETE FROM active_text_mutes WHERE discord_snowflake=$1', uid)

    @commands.command(name='backup', help='Creates a backup of the database and uploads it')
    @is_owner_developer_predicator()
    async def backup(self, ctx: commands.Context):
        try:
            backup_dir = self.setup_backup_directory('/app/backups')
            backup_file = self.perform_backup(
                db_user=os.getenv('POSTGRES_USER'),
                db_name=os.getenv('POSTGRES_DATABASE'),
                db_host=os.getenv('POSTGRES_HOST'),
                db_password=os.getenv('POSTGRES_PASSWORD'),
                backup_dir=backup_dir
            )
            logger.info(f'Backup completed successfully: {backup_file}')
            await ctx.send(file=discord.File(backup_file))
        except Exception as e:
            logger.error(f'Error during database backup: {e}')
            await ctx.send(f'\U0001F6AB Failed to create backup: {e}')
    
    @commands.command(name='cap', help='Set a duration limit for bans, mutes and text mutes.')
    @is_owner_developer_predicator()
    async def cap(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        moderation_type: Optional[str] = commands.parameter(default=None, description='One of: `mute`, `ban`, `tmute`'),
        *,
        duration: Optional[str] = commands.parameter(default='24', description='Options: (+|-)duration(m|h|d) \n 0 - permanent / 24h - default')
    ):
        channel_obj = await self.resolve_channel(ctx, channel)
        if not channel_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        valid_types = {'mute', 'ban', 'tmute'}
        if moderation_type not in valid_types:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid moderation type. Must be one of: {', '.join(valid_types)}')
        expires_at, duration_str = self.parse_duration(duration)
        original_duration = await self.get_cap(channel_obj.id, ctx.guild.id, moderation_type)
        await self.set_cap(channel_obj.id, ctx.guild.id, moderation_type, duration)
        if original_duration:
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Cap changed on {channel_obj.mention} for {moderation_type} from {duration} to {duration_str}.')
        else:
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Cap set on {channel_obj.mention} for {moderation_type} {duration_str}.')
        
    @commands.command(name='admin', help='Grants server mute privileges to a member for the entire guild.')
    @is_owner_predicator()
    async def create_administrator(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(description='Tag a member or include their snowflake ID')
    ):
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (discord_snowflake, server_muter_guild_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[]) ON CONFLICT (discord_snowflake) DO
                UPDATE
                    SET server_muter_guild_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(COALESCE (u.server_muter_guild_ids, '{}') || ARRAY[$2])
                        )
                        FROM users u WHERE u.discord_snowflake = EXCLUDED.discord_snowflake
                     ),
                    updated_at = NOW()
            ''', member_obj.id, ctx.guild.id)
        self.server_muters.setdefault(ctx.guild.id, set()).add(member_obj.id)
        return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} has been granted server mute permissions.', allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='alias', help='Set an alias for a cow, uncow, mute, unmute, ban, unban, flag, unflag, tmute, untmute, role, or unrole action.')
    @is_owner_developer_coordinator_predicator(None)
    async def create_alias(
        self,
        ctx: commands.Context,
        alias_type: Optional[str] = commands.parameter(default=None, description='One of: `cow`, `uncow`, `mute`, `unmute`, `ban`, `unban`, `flag`, `unflag`, `tmute`, `untmute`, `role`, `unrole`'),
        alias_name: Optional[str] = commands.parameter(default=None, description='Alias/Pseudonym'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        *,
        role: Optional[str] = commands.parameter(default=None, description='Role ID (only for role/unrole)')
    ) -> None:
        cmd = None
        if alias_type:
            alias_type = alias_type.lower()
        valid_types = {'cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute', 'role', 'unrole'}
        if alias_type not in valid_types:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid alias type. Must be one of: `{"`, `".join(valid_types)}`')
        if not alias_name.strip():
            return await self.handler.send_message(ctx, content='\U0001F6AB Alias name cannot be empty.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if not channel_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        is_owner_or_dev, _ = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev:
            async with ctx.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow('SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1', ctx.author.id)
            if not row or channel_obj.id not in (row.get('coordinator_channel_ids') or []):
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`alias`) in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            existing_alias = await conn.fetchrow('''
                SELECT guild_id, channel_id, role_id
                FROM command_aliases
                WHERE alias_type = $1
                    AND alias_name = $2
            ''', alias_type, alias_name)
            if existing_alias:
                existing_channel = ctx.guild.get_channel(existing_alias['channel_id']) if existing_alias['channel_id'] else None
                existing_role = ctx.guild.get_role(existing_alias['role_id']) if existing_alias['role_id'] else None
                mention = existing_channel.mention if existing_channel else (existing_role.mention if existing_role else 'unknown')
                return await ctx.send(f'\U0001F6AB Alias `{alias_name}` ({alias_type}) already exists and is set to {mention}.',     allowed_mentions=discord.AllowedMentions.none())
            if self.bot.get_command(alias_name):
                return await self.handler.send_message(ctx, content=f'\U0001F6AB A command named `{alias_name}` already exists.')
            if alias_type in ('role', 'unrole'):
                if not role:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Role ID is required for role/unrole aliases.')
                try:
                    role_id = int(role.replace('<@&', '').replace('>', ''))
                except ValueError:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid role ID: {role_id}')
                await conn.execute('''
                    INSERT INTO command_aliases (guild_id, alias_type, alias_name, role_id, channel_id)
                    VALUES ($1, $2, $3, $4, $5)
                ''', ctx.guild.id, alias_type, alias_name, role_id, channel_obj.id)
            else:
                await conn.execute('''
                    INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id)
                    VALUES ($1, $2, $3, $4)
                ''', ctx.guild.id, alias_type, alias_name, channel_obj.id)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1, $2, $3, $4, $5, $6)
            ''', 'create_alias', None, ctx.author.id, ctx.guild.id, channel_obj.id, f'Created an alias: {alias_name}')
        if alias_type in ('role', 'unrole'):
            if is_owner_or_dev:
                self.bot.command_aliases.setdefault(ctx.guild.id, {}).setdefault('role_aliases', {}).setdefault(alias_type, {})[alias_name] = {
                    'channel_id': int(channel_obj.id),
                    'role_id': int(role_id)
                }
        else:
            self.bot.command_aliases.setdefault(ctx.guild.id, {}).setdefault('channel_aliases', {}).setdefault(alias_type, {})[alias_name] = channel_obj.id
        if alias_type == 'ban':
            cmd = self.create_ban_alias(alias_name)
        elif alias_type == 'flag':
            cmd = self.create_flag_alias(alias_name)
        elif alias_type == 'unflag':
            cmd = self.create_unflag_alias(alias_name)
        elif alias_type == 'mute':
            cmd = self.create_voice_mute_alias(alias_name)
        elif alias_type == 'tmute':
            cmd = self.create_text_mute_alias(alias_name)
        elif alias_type == 'unban':
            cmd = self.create_unban_alias(alias_name)
        elif alias_type == 'cow':
            if channel_obj.id != '1222056499959042108':
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {channel_obj.mention} has no vegan association.')
            cmd = self.create_cow_alias(alias_name)
        elif alias_type == 'uncow':
            if channel_obj.id != '1222056499959042108':
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {channel_obj.mention} has no vegan association.')
            cmd = self.create_uncow_alias(alias_name)
        elif alias_type == 'unmute':
            cmd = self.create_unmute_alias(alias_name)
        elif alias_type == 'untmute':
            cmd = self.create_untextmute_alias(alias_name)
        elif alias_type == 'role':
            cmd = self.create_role_alias(alias_name)
        elif alias_type == 'unrole':
            cmd = self.create_unrole_alias(alias_name)
        self.bot.add_command(cmd)
        if alias_type in ('role', 'unrole'):
            role_obj = ctx.guild.get_role(int(role_id))
            mention = role_obj.mention if role else f'<@&{role_id}>'
        else:
            mention = channel_obj.mention
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Alias `{alias_name}` ({alias_type}) set to {mention}.')

    def create_ban_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Ban a user from a voice channel.')
        @is_owner_developer_coordinator_moderator_predicator('ban')
        async def ban_alias(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            duration: Optional[str] = commands.parameter(default='24', description='Options: (+|-)duration(m|h|d) \n 0 - permanent / 24h - default'),
            *,
            reason: Optional[str] = commands.parameter(default='', description='Optional reason (required for 7 days or more)')
        ) -> None:
            cmd = ctx.invoked_with
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            static_channel_id = int(
                self.bot.command_aliases
                    .get(ctx.guild.id, {})
                    .get('channel_aliases', {})
                    .get('ban', {})
                    .get(command_name)
            )
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No channel alias mapping found for `{cmd}`.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from the alias.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot ban the bot.')
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permissions to this command (`{command_name}`) in {channel_obj.mention}.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to ban this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                existing_ban = await conn.fetchrow('''
                    SELECT expires_at, reason
                    FROM active_bans
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                ''', ctx.guild.id, member_obj.id, channel_obj.id)
                is_modification = existing_ban is not None
                base_time = existing_ban['expires_at'] if existing_ban else None
                stripped = duration.strip() if duration else ''
                if stripped in ('+', '-', '='):
                    expires_at = base_time
                    duration_display = self.fmt_duration(base_time)
                else:
                    expires_at, duration_display = self.parse_duration(duration, base=base_time)
                is_coordinator = await is_owner_developer_coordinator_via_alias(ctx, 'ban')
                is_relative_duration = stripped.startswith('+') and (len(stripped) > 1 and stripped[1].isdigit())
                is_reason_append = stripped == '+' or (stripped.startswith('+') and not stripped[1].isdigit())
                is_reason_set = stripped.startswith('=')
                is_reason_delete = stripped == '-'
                updated_reason = existing_ban['reason'] if existing_ban else None
                if existing_ban and (is_reason_append or is_reason_set or is_reason_delete):
                    if is_reason_append:
                        new_text = reason.strip() if reason else ''
                        if not new_text:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason to append.')
                        updated_reason = f"{updated_reason}\n{new_text}" if updated_reason else new_text
                    elif is_reason_set:
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset ban reasons.')
                        updated_reason = reason.strip() if reason else ''
                        if not updated_reason:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason after "=" to set.')
                    elif is_reason_delete:
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset ban reasons.')
                        updated_reason = None
                if is_modification and not is_coordinator and not is_relative_duration and not (is_reason_append or is_reason_set or is_reason_delete):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can overwrite an existing ban with an absolute duration.')
                caps = await self.get_caps_for_channel(ctx.guild.id, channel_obj.id)
                active_cap = next((c for c in caps if c[0] == 'ban'), None)
                now = datetime.now(timezone.utc)
                if active_cap:
                    cap_expires_at, _ = self.parse_duration(active_cap[1])
                    if cap_expires_at is None or (expires_at and expires_at > cap_expires_at):
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content=f'\U0001F6AB Only coordinators can create bans longer than the channel cap ({active_cap[1]}).')
                        if not reason.strip() and not is_reason_set:
                            return await self.handler.send_message(ctx, content=f'\U0001F6AB A reason is required for bans longer than the channel cap ({active_cap[1]}).')
                elif expires_at is None or (expires_at - now) > timedelta(days=7):
                    if not is_coordinator:
                        return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can ban permanently or longer than 7 days.')
                    if not reason.strip() and not is_reason_set:
                        return await self.handler.send_message(ctx, content='\U0001F6AB A reason is required for permanent bans or those longer than 7 days.')
                    if existing_ban and (is_relative_duration or is_coordinator) and not (is_reason_append or is_reason_set or is_reason_delete):
                        duration_str = stripped.lower()
                        valid_relative = duration_str in ('0','0h','0d','0m') or (duration_str.startswith('+') and duration_str[1].isdigit()) or duration_str.startswith('-')
                        if not valid_relative:
                            if existing_ban['expires_at'] is None:
                                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already permanently banned in {channel_obj.mention}.')
                            else:
                                remaining = existing_ban['expires_at'] - discord.utils.utcnow()
                                if remaining.total_seconds() > 0:
                                    hours_left = round(remaining.total_seconds() / 3600, 1)
                                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already banned from {channel_obj.mention} for another {hours_left}h.')
                elif existing_ban:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already banned in {channel_obj.mention}.')
            try:
                await channel_obj.set_permissions(
                    member_obj,
                    view_channel=False,
                    reason=f'{reason or "No reason provided"}'
                )
            except discord.Forbidden:
                logger.warning('\U0001F6AB Missing permissions to deny channel access.')
            is_in_channel = False
            if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                is_in_channel = True
                try:
                    await member_obj.move_to(None, reason=f'{reason or "No reason provided"}')
                except discord.Forbidden:
                    await ctx.send(f'\U0001F6AB Could not disconnect {member_obj.mention} from {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                except Exception as e:
                    logger.exception(f'Unexpected error while disconnecting user: {e}')
                    raise
            try:
                async with self.bot.db_pool.acquire() as conn:
                    await conn.execute('''
                        INSERT INTO active_bans (guild_id, discord_snowflake, channel_id, expires_at, reason)
                        VALUES ($1, $2, $3, $4, $5)
                        ON CONFLICT (guild_id, discord_snowflake, channel_id) DO UPDATE
                        SET expires_at = EXCLUDED.expires_at,
                            reason = EXCLUDED.reason
                    ''', ctx.guild.id, member_obj.id, channel_obj.id, expires_at, reason or 'No reason provided')
                    await conn.execute(
                    '''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'ban', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, reason or 'No reason provided')
            except Exception as e:
                logger.warning(f'Database error occurred: {e}')
                raise
            embed = discord.Embed(
                title=f"{self.get_random_emoji()} {member_obj.display_name} has been banned",
                description=f"**Channel:** {channel_obj.mention}\n**Duration:** {duration_display}\n**Reason:** {updated_reason or 'No reason provided'}",
                color=discord.Color.orange()
            )
            await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.none())
            highest_role = await self.get_highest_role(ctx, ctx.author, channel_obj)
            await self.send_log(ctx, 'ban', member_obj, channel_obj, duration_display, reason or 'No reason provided', ctx.author, expires_at, command_name, is_in_channel, is_modification, highest_role)
        return ban_alias
        
    @commands.command(name='coord', help='Grants coordinator access for a specific voice channel.')
    @is_owner_developer_predicator()
    async def create_coordinator(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
    ) -> None:
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj or not member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if not channel_obj or not channel:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        is_owner_or_dev, _ = await check_owner_dev_coord(ctx, channel_obj)
        if not is_owner_or_dev:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permissions to use this command (`{command_name}`) in {ctx.guild.name}')
        if member_obj.bot and not is_owner_or_dev:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot a coordinator.')
        highest_role, success = await check_block(ctx, member_obj, channel_obj)
        if not success:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to make this `{highest_role}` a coordinator because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (discord_snowflake, coordinator_channel_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[])
                ON CONFLICT (discord_snowflake) DO UPDATE
                SET 
                    coordinator_channel_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(
                                COALESCE(users.coordinator_channel_ids, ARRAY[]::BIGINT[]) || 
                                ARRAY[$2]::BIGINT[]
                            )
                        )
                    ),
                    updated_at = NOW()
            ''', member_obj.id, channel_obj.id)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1, $2, $3, $4, $5, $6)
            ''', 'create_coordinator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Created a coordinator')
        return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} has been granted coordinator rights in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
        
    def create_cow_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Label a user as going vegan for tracking purposes.')
        @is_owner_developer_coordinator_moderator_predicator('cow')
        async def going_vegan_alias(
                ctx: commands.Context,
                member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            channel_id = (
                self.bot.command_aliases
                    .get(ctx.guild.id, {})
                    .get('channel_aliases', {})
                    .get('cow', {})
                    .get(command_name)
            )
            channel_obj = await self.resolve_channel(ctx, channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot cow the bot.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to cow this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            select_sql = '''
                SELECT 1
                FROM active_cows
                WHERE discord_snowflake = $1
                AND channel_id = $2
            '''
            insert_sql = '''
                INSERT INTO active_cows (guild_id, discord_snowflake, channel_id)
                VALUES ($1, $2, $3)
                ON CONFLICT (guild_id, discord_snowflake, channel_id) DO NOTHING
            '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    already_cowed = await conn.fetchval(select_sql, member_obj.id, channel_obj.id)
                    if already_cowed:
                        return await ctx.send(f'\U0001F6AB {member_obj.mention} is already going vegan.', allowed_mentions=discord.AllowedMentions.none())
                    await conn.execute(insert_sql, ctx.guild.id, member_obj.id, channel_obj.id)
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'cow', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Cowed a user')
                    await ctx.send(f'\U0001F525 {member_obj.mention} is going vegan!!! \U0001F525', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                return await self.handler.send_message(ctx, content=f'Database error: {e}')
                raise
        return going_vegan_alias
        
    @commands.command(name='dev', help='Elevates a user\'s permissions to a bot developer.')
    @is_owner_predicator()
    async def create_developer(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
    ) -> None:
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj or not member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
             return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot a developer.')
        highest_role, success = await check_block(ctx, member_obj, None)
        if not success:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to make this `{highest_role}` a developer because they are a higher/or equivalent role than you in {ctx.guild.name}.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (discord_snowflake, developer_guild_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[])
                ON CONFLICT (discord_snowflake) DO UPDATE
                SET developer_guild_ids = (
                    SELECT ARRAY(
                        SELECT DISTINCT unnest(u.developer_guild_ids || EXCLUDED.developer_guild_ids)
                    )
                    FROM users u WHERE u.discord_snowflake = EXCLUDED.discord_snowflake
                ),
                updated_at = NOW()
            ''', member_obj.id, ctx.guild.id)
        return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} has been granted developer rights in {ctx.guild.name}.', allowed_mentions=discord.AllowedMentions.none())

    def create_flag_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Flag a user in the database for the channel mapped to this alias.')
        @is_owner_developer_coordinator_moderator_predicator('flag')
        async def flag_alias(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            *,
            reason: Optional[str] = commands.parameter(default='', description='Optional reason')
        ) -> None:
            channel_id = (
                self.bot.command_aliases
                    .get(ctx.guild.id, {})
                    .get('channel_aliases', {})
                    .get('flag', {})
                    .get(command_name)
            )
            channel_obj = await self.resolve_channel(ctx, channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot flag the bot.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to flag this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            select_sql = '''
                SELECT reason
                FROM active_flags
                WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
            '''
            insert_sql = '''
                INSERT INTO active_flags (guild_id, discord_snowflake, channel_id, reason)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_id, discord_snowflake, channel_id) DO NOTHING
            '''
            update_sql = '''
                UPDATE active_flags
                SET reason = $4
                WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
            '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    existing_flag = await conn.fetchrow(select_sql, ctx.guild.id, member_obj.id, channel_obj.id)
                    is_modification = existing_flag is not None
                    updated_reason = existing_flag['reason'] if existing_flag else None
                    stripped = reason.strip() if reason else ''
                    is_reason_append = stripped == '+' or (stripped.startswith('+') and not stripped[1].isdigit())
                    is_reason_set = stripped.startswith('=')
                    is_reason_delete = stripped == '-'
                    if is_modification and (is_reason_append or is_reason_set or is_reason_delete):
                        if is_reason_append:
                            new_text = reason.strip() if reason else ''
                            if not new_text:
                                return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason to append.')
                            updated_reason = f"{updated_reason}\n{new_text}" if updated_reason else new_text
                        elif is_reason_set:
                            is_coordinator = await is_owner_developer_coordinator_via_alias(ctx, 'flag')
                            if not is_coordinator:
                                return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset flag reasons.')
                            updated_reason = reason.strip() if reason else ''
                            if not updated_reason:
                                return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason after "=" to set.')
                        elif is_reason_delete:
                            is_coordinator = await is_owner_developer_coordinator_via_alias(ctx, 'flag')
                            if not is_coordinator:
                                return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can delete flag reasons.')
                            updated_reason = None
                        await conn.execute(update_sql, ctx.guild.id, member_obj.id, channel_obj.id, updated_reason)
                    elif not is_modification:
                        await conn.execute(insert_sql, ctx.guild.id, member_obj.id, channel_obj.id, reason if reason else 'No reason provided')
                        updated_reason = reason if reason else 'No reason provided'
                    else:
                        return await ctx.send(f'\U0001F6AB {member_obj.mention} is already flagged in {channel_obj.mention} for {existing_flag["reason"]}.', allowed_mentions=discord.AllowedMentions.none())
                    embed = discord.Embed(color=discord.Color.orange())
                    embed.set_author(name=f"{member_obj.mention} is flagged", icon_url=member_obj.display_avatar.url)
                    embed.add_field(name="Channel", value=channel_obj.mention, inline=True)
                    embed.add_field(name="Reason", value=updated_reason or 'No reason provided', inline=False)
                    await ctx.send(embed=embed)
            except Exception as e:
                logger.exception(f'Database error in flag_alias: {e}')
                raise
    
        return flag_alias
        
    @commands.command(name='mod', help='Elevates a user\'s permission to VC moderator for a specific channel.')
    @is_owner_developer_coordinator_predicator(None)
    async def create_moderator(
            self,
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ) -> None:
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj or not member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if not channel_obj or not channel:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        is_owner_or_dev, is_coord = await check_owner_dev_coord(ctx, channel_obj)
        if not is_owner_or_dev and not is_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
        if member_obj.bot and not is_owner_or_dev:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot a moderator.')
        highest_role, success = await check_block(ctx, member_obj, channel_obj)
        if not success:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to make this `{highest_role}` a moderator because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        if not is_owner_or_dev:
            async with self.bot.db_pool.acquire() as conn:
                coordinator_row = await conn.fetchrow('''
                    SELECT 1
                    FROM users
                    WHERE discord_snowflake = $1
                      AND $2 = ANY (coordinator_channel_ids)
                ''', ctx.author.id, channel_obj.id)
                if not coordinator_row:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not a coordinator in {channel_obj.mention} and cannot assign moderators there.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (discord_snowflake, moderator_channel_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[])
                ON CONFLICT (discord_snowflake) DO UPDATE
                SET 
                    moderator_channel_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(
                                COALESCE(users.moderator_channel_ids, ARRAY[]::BIGINT[]) ||
                                ARRAY[$2]::BIGINT[]
                            )
                        )
                    ),
                    updated_at = NOW()
            ''', member_obj.id, channel_obj.id)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1, $2, $3, $4, $5, $6)
            ''', 'create_moderator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Created a moderator')
        return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} has been granted moderator rights in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())

    def create_role_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help=f'Gives a specific role to a user.')
        @is_owner_developer_coordinator_predicator('role')
        async def role_alias(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            alias_data = (
                self.bot.command_aliases
                    .get(ctx.guild.id, {})
                    .get('role_aliases', {})
                    .get('role', {})
                    .get(command_name)
            )
            if not alias_data:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No role alias configured for `{command_name}`.')
            static_role_id = int(alias_data.get('role_id'))
            target_channel_id = int(alias_data.get('channel_id')) if alias_data.get('channel_id') else None
            is_owner_or_dev, is_coord = await check_owner_dev_coord(ctx, target_channel_id)
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            channel_obj = await self.resolve_channel(ctx, target_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            if not is_owner_or_dev and not is_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot give the bot a role.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to make this `{highest_role}` have a role because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            role_obj = ctx.guild.get_role(static_role_id)
            if not role_obj:
                return await ctx.send(f'⚠️ Could not resolve role with ID `{static_role_id}`.')
            if role_obj in member_obj.roles:
                return await ctx.send(f'{member_obj.mention} already has {role_obj.mention}.')
            await member_obj.add_roles(role_obj, reason=f'Added role')
            return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} was given {role_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
        return role_alias
    
    def create_text_mute_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Text mutes a user in a specific text channel.')
        @is_owner_developer_coordinator_moderator_predicator('tmute')
        async def text_mute_alias(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            duration: Optional[str] = commands.parameter(default='8', description='Options: (+|-)duration(m|h|d) \n 0 - permanent / 8h - default'),
            *,
            reason: Optional[str] = commands.parameter(default='', description='Optional reason (required for 7 days or more)')
        ) -> None:
            static_channel_id = int(
                self.bot.command_aliases
                    .get(ctx.guild.id, {})
                    .get('channel_aliases', {})
                    .get('tmute', {})
                    .get(command_name)
            )
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No text mute alias configured for `{command_name}`.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot text mute the bot.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to text-mute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                existing_text_mute = await conn.fetchrow('''
                    SELECT expires_at, reason
                    FROM active_text_mutes
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                ''', ctx.guild.id, member_obj.id, channel_obj.id)
                is_modification = existing_text_mute is not None
                base_time = existing_text_mute['expires_at'] if existing_text_mute else None
                stripped = duration.strip() if duration else ''
                if stripped in ('+', '-', '='):
                    expires_at = base_time
                    duration_display = self.fmt_duration(base_time)
                else:
                    expires_at, duration_display = self.parse_duration(duration, base=base_time)
                is_coordinator = await is_owner_developer_coordinator_via_alias(ctx, 'tmute')
                is_relative_duration = stripped.startswith('+') and (len(stripped) > 1 and stripped[1].isdigit())
                is_reason_append = stripped == '+' or (stripped.startswith('+') and not stripped[1].isdigit())
                is_reason_set = stripped.startswith('=')
                is_reason_delete = stripped == '-'
                updated_reason = existing_text_mute['reason'] if existing_text_mute else None
                if existing_text_mute and (is_reason_append or is_reason_set or is_reason_delete):
                    if is_reason_append:
                        new_text = reason.strip() if reason else ''
                        if not new_text:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason to append.')
                        updated_reason = f"{updated_reason}\n{new_text}" if updated_reason else new_text
                    elif is_reason_set:
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset text mute reasons.')
                        updated_reason = reason.strip() if reason else ''
                        if not updated_reason:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason after "=" to set.')
                    elif is_reason_delete:
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset text mute reasons.')
                        updated_reason = None
                if is_modification and not is_coordinator and not is_relative_duration and not (is_reason_append or is_reason_set or is_reason_delete):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can overwrite an existing text mute with an absolute duration.')
                caps = await self.get_caps_for_channel(ctx.guild.id, channel_obj.id)
                active_cap = next((c for c in caps if c[0] == 'tmute'), None)
                now = datetime.now(timezone.utc)
                if active_cap:
                    cap_expires_at, _ = self.parse_duration(active_cap[1])
                    if cap_expires_at is None or (expires_at and expires_at > cap_expires_at):
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content=f'\U0001F6AB Only coordinators can create text mutes longer than the channel cap ({active_cap[1]}).')
                        if not reason.strip() and not is_reason_set:
                            return await self.handler.send_message(ctx, content=f'\U0001F6AB A reason is required for text mutes longer than the channel cap ({active_cap[1]}).')
                elif expires_at is None or (expires_at - now) > timedelta(days=7):
                    if not is_coordinator:
                        return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can text mute permanently or longer than 7 days.')
                    if not reason.strip() and not is_reason_set:
                        return await self.handler.send_message(ctx, content='\U0001F6AB A reason is required for permanent text mutes or those longer than 7 days.')
                    if existing_text_mute and (is_relative_duration or is_coordinator) and not (is_reason_append or is_reason_set or is_reason_delete):
                        duration_str = stripped.lower()
                        valid_relative = duration_str in ('0','0h','0d','0m') or (duration_str.startswith('+') and duration_str[1].isdigit()) or duration_str.startswith('-')
                        if not valid_relative:
                            if existing_text_mute['expires_at'] is None:
                                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already permanently text muted in {channel_obj.mention}.')
                            else:
                                remaining = existing_text_mute['expires_at'] - discord.utils.utcnow()
                                if remaining.total_seconds() > 0:
                                    hours_left = round(remaining.total_seconds() / 3600, 1)
                                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already text muted from {channel_obj.mention} for another {hours_left}h.')
                elif existing_text_mute:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already text muted in {channel_obj.mention}.')
                try:
                    await channel_obj.set_permissions(member_obj, send_messages=False, add_reactions=False)
                except discord.Forbidden:
                    logger.warning('\U0001F6AB The user\'s channel permissions were unable to be updated.')
                try:
                    await conn.execute('''
                        INSERT INTO active_text_mutes (guild_id, discord_snowflake, channel_id, reason, expires_at)
                        VALUES ($1, $2, $3, $4, $5)
                        ON CONFLICT (guild_id, discord_snowflake, channel_id) DO UPDATE
                        SET reason = $4,
                            expires_at = $5
                    ''', ctx.guild.id, member_obj.id, channel_obj.id, updated_reason or 'No reason provided', expires_at)
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'textmute', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, f'Textmuted a user ({updated_reason or "No reason provided"})')
                except Exception as e:
                    logger.warning(f'DB insert failed: {e}')
                    return await self.handler.send_message(ctx, content=str(e))
                embed = discord.Embed(
                    title=f"{self.get_random_emoji()} {member_obj.display_name} is text-muted",
                    description=f"**Channel:** {channel_obj.mention}\n**Duration:** {duration_display}\n**Reason:** {updated_reason or 'No reason provided'}",
                    color=discord.Color.orange()
                )
                await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.none())
                highest_role = await self.get_highest_role(ctx, ctx.author, channel_obj)
                await self.send_log(ctx, 'tmute', member_obj, channel_obj, duration_display, updated_reason or 'No reason provided', ctx.author, expires_at, command_name, True, is_modification, highest_role)
        return text_mute_alias

    def create_voice_mute_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Mutes a member in a specific VC.')
        @is_owner_developer_coordinator_moderator_predicator('mute')
        async def voice_mute_alias(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            duration: Optional[str] = commands.parameter(default='8', description='Options: (+|-)duration(m|h|d) \n 0 - permanent / 8h - default'),
            *,
            reason: Optional[str] = commands.parameter(default='', description='Optional reason (required for 7 days or more)')
        ) -> None:
            static_channel_id = int(
                self.bot.command_aliases
                    .get(ctx.guild.id, {})
                    .get('channel_aliases', {})
                    .get('mute', {})
                    .get(command_name)
            )
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No text unmute alias configured for `{command_name}`.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}.')
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot voice mute the bot.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to mute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                existing_mute = await conn.fetchrow('''
                    SELECT expires_at, reason
                    FROM active_voice_mutes
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                ''', ctx.guild.id, member_obj.id, channel_obj.id)
                is_modification = existing_mute is not None
                base_time = existing_mute['expires_at'] if existing_mute else None
                stripped = duration.strip() if duration else ''
                if stripped in ('+', '-', '='):
                    expires_at = base_time
                    duration_display = self.fmt_duration(base_time)
                else:
                    expires_at, duration_display = self.parse_duration(duration, base=base_time)
                is_coordinator = await is_owner_developer_coordinator_via_alias(ctx, 'mute')
                is_relative_duration = stripped.startswith('+') and (len(stripped) > 1 and stripped[1].isdigit())
                is_reason_append = stripped == '+' or (stripped.startswith('+') and not stripped[1].isdigit())
                is_reason_set = stripped.startswith('=')
                is_reason_delete = stripped == '-'
                updated_reason = existing_mute['reason'] if existing_mute else None
                if existing_mute and (is_reason_append or is_reason_set or is_reason_delete):
                    if is_reason_append:
                        new_text = reason.strip() if reason else ''
                        if not new_text:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason to append.')
                        updated_reason = f"{updated_reason}\n{new_text}" if updated_reason else new_text
                    elif is_reason_set:
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset mute reasons.')
                        updated_reason = reason.strip() if reason else ''
                        if not updated_reason:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason after "=" to set.')
                    elif is_reason_delete:
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset voice mute reasons.')
                        updated_reason = None
                if is_modification and not is_coordinator and not is_relative_duration and not (is_reason_append or is_reason_set or is_reason_delete):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can overwrite an existing voice mute with an absolute duration.')
                caps = await self.get_caps_for_channel(ctx.guild.id, channel_obj.id)
                active_cap = next((c for c in caps if c[0] == 'mute'), None)
                now = datetime.now(timezone.utc)
                if active_cap:
                    cap_expires_at, _ = self.parse_duration(active_cap[1])
                    if cap_expires_at is None or (expires_at and expires_at > cap_expires_at):
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content=f'\U0001F6AB Only coordinators can create voice mutes longer than the channel cap ({active_cap[1]}).')
                        if not reason.strip() and not is_reason_set:
                            return await self.handler.send_message(ctx, content=f'\U0001F6AB A reason is required for voice mutes longer than the channel cap ({active_cap[1]}).')
                elif expires_at is None or (expires_at - now) > timedelta(days=7):
                    if not is_coordinator:
                        return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can ban permanently or longer than 7 days.')
                    if not reason.strip() and not is_reason_set:
                        return await self.handler.send_message(ctx, content='\U0001F6AB A reason is required for permanent voice mutes or those longer than 7 days.')
                    if existing_mute and (is_relative_duration or is_coordinator) and not (is_reason_append or is_reason_set or is_reason_delete):
                        duration_str = stripped.lower()
                        valid_relative = duration_str in ('0','0h','0d','0m') or (duration_str.startswith('+') and duration_str[1].isdigit()) or duration_str.startswith('-')
                        if not valid_relative:
                            if existing_mute['expires_at'] is None:
                                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already permanently banned in {channel_obj.mention}.')
                            else:
                                remaining = existing_mute['expires_at'] - discord.utils.utcnow()
                                if remaining.total_seconds() > 0:
                                    hours_left = round(remaining.total_seconds() / 3600, 1)
                                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already voice muted in {channel_obj.mention} for another {hours_left}h.')
                elif existing_mute:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already voice muted in {channel_obj.mention}.')
            try:
                async with self.bot.db_pool.acquire() as conn:
                    await conn.execute('''
                        INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, reason)
                        VALUES ($1, $2, $3, $4, $5)
                        ON CONFLICT (guild_id, discord_snowflake, channel_id) DO UPDATE
                        SET expires_at = EXCLUDED.expires_at,
                            reason = EXCLUDED.reason
                    ''', ctx.guild.id, member_obj.id, channel_obj.id, expires_at, reason or 'No reason provided')
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'voice_mute', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Voice muted a member')
            except Exception as e:
                logger.warning(f'DB insert failed: {e}')
                return await self.handler.send_message(ctx, content=str(e))
            is_in_channel = False
            if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                is_in_channel = True
                await member_obj.edit(mute=True)
            embed = discord.Embed(
                    title=f"{self.get_random_emoji()} {member_obj.display_name} is voice muted",
                    description=f"**Channel:** {channel_obj.mention}\n**Duration:** {duration_display}\n**Reason:** {updated_reason or 'No reason provided'}",
                    color=discord.Color.orange()
                )
            await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.none())
            highest_role = await self.get_highest_role(ctx, ctx.author, channel_obj)
            await self.send_log(ctx, 'vmute', member_obj, channel_obj, duration_display, reason or 'No reason provided', ctx.author, expires_at, command_name, is_in_channel, is_modification, highest_role)
        return voice_mute_alias


    def create_unban_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Unban a user from a voice channel.')
        @is_owner_developer_coordinator_moderator_predicator('unban')
        async def unban_alias(
            ctx,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            static_channel_id = int(self.bot.command_aliases.get(ctx.guild.id, {}).get('channel_aliases', {}).get('unban', {}).get(command_name))
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No unban alias configured for `{command_name}`.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel for `{command_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot unban the bot.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to unban this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow('''
                    SELECT expires_at FROM active_bans
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                ''', ctx.guild.id, member_obj.id, channel_obj.id)
                if row and row['expires_at'] is None and not await is_owner_developer_coordinator_via_alias(ctx, 'unban'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Coordinator-only for undoing permanent bans.')
            try:
                await channel_obj.set_permissions(member_obj, overwrite=None)
            except discord.Forbidden:
                logger.warning('\U0001F6AB Missing permissions to update channel permissions.')
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute('DELETE FROM active_bans WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3', ctx.guild.id, member_obj.id, channel_obj.id)
                await conn.execute('INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1, $2, $3, $4, $5, $6)', 'unban', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Unbanned a user')
            return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} has been unbanned from {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
        return unban_alias

    def create_uncow_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Unlabel a user for tracking purposes.')
        @is_owner_developer_coordinator_moderator_predicator('uncow')
        async def no_longer_going_vegan_alias(
                ctx: commands.Context,
                member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            static_channel_id = int(
                self.bot.command_aliases
                    .get(ctx.guild.id, {})
                    .get('channel_aliases', {})
                    .get('uncow', {})
                    .get(command_name)
            )
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No uncow alias configured for `{command_name}`.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel for `{command_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot uncow the bot.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to uncow this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            select_sql = '''
                SELECT 1
                FROM active_cows
                WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
            '''
            update_sql = '''
                DELETE FROM active_cows
                WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
            '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    is_flagged = await conn.fetchval(select_sql, ctx.guild.id, member_obj.id, channel_obj.id)
                    if not is_flagged:
                        return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} has no active record in {channel_obj.mention}.')
                    await conn.execute(update_sql, ctx.guild.id, member_obj.id, channel_obj.id)
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'uncow', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Uncowed a user')
                    await ctx.send(f'👎<@{member_obj.id}> is no longer going vegan.👎', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                await self.handler.send_message(ctx, content=f'Database error: {e}')
                raise
        return no_longer_going_vegan_alias
        
    def create_unflag_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Unflag a user in the database for the voice channel mapped to this alias.')
        @is_owner_developer_coordinator_moderator_predicator('unflag')
        async def unflag_alias(
                ctx: commands.Context,
                member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            static_channel_id = int(
                self.bot.command_aliases
                    .get(ctx.guild.id, {})
                    .get('channel_aliases', {})
                    .get('unflag', {})
                    .get(command_name)
            )
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No unflag alias configured for `{command_name}`.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel for `{command_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot unflag the bot.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to unflag this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            select_sql = '''
                SELECT 1
                FROM active_flags
                WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
            '''
            update_sql = '''
                DELETE FROM active_flags
                WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
            '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    is_flagged = await conn.fetchval(select_sql, ctx.guild.id, member_obj.id, channel_obj.id)
                    if not is_flagged:
                        return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not flagged for {channel_obj.mention}.')
                    await conn.execute(update_sql, ctx.guild.id, member_obj.id, channel_obj.id)
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'unflag', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Unflagged a user')
                    await ctx.send(f'{self.get_random_emoji()} Unflagged {member_obj.mention} for channel {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                await self.handler.send_message(ctx, content=f'Database error: {e}')
                raise
        return unflag_alias


    def create_unmute_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Unmutes a member in a specific VC.')
        @is_owner_developer_coordinator_moderator_predicator('unmute')
        async def unmute_alias(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            static_channel_id = int(self.bot.command_aliases.get(ctx.guild.id, {}).get('channel_aliases', {}).get('unmute', {}).get(command_name))
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No voice unmute alias configured for `{command_name}`.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel for `{command_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot unmute the bot.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to unmute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            try:
                async with self.bot.db_pool.acquire() as conn:
                    row = await conn.fetchrow('''
                        SELECT expires_at FROM active_voice_mutes
                        WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                    ''', ctx.guild.id, member_obj.id, channel_obj.id)
                    if not row:
                        return await ctx.send(f'\U0001F6AB {member_obj.mention} is not muted in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                    if row['expires_at'] is None and not await is_owner_developer_coordinator_via_alias(ctx, 'unmute'):
                        return await self.handler.send_message(ctx, content='\U0001F6AB Coordinator-only for undoing permanent voice mutes.')
                    if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                        await member_obj.edit(mute=False)
                    await conn.execute('DELETE FROM active_voice_mutes WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3', ctx.guild.id, member_obj.id,  channel_obj.id)
                    await conn.execute('INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1,  $2, $3, $4, $5, $6)', 'unmute', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Unmuted a member')
            except Exception as e:
                logger.warning(f'\U0001F6AB Database error: {e}')
                raise
            if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} has been unmuted in {channel_obj.mention}.',  allowed_mentions=discord.AllowedMentions.none())
            return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} is no longer marked as muted in {channel_obj.mention}.',  allowed_mentions=discord.AllowedMentions.none())
        return unmute_alias

    def create_unrole_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Removes a specific role from a user.')
        @is_owner_developer_coordinator_predicator('unrole')
        async def unrole_alias(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            alias_data = (
                self.bot.command_aliases
                    .get(ctx.guild.id, {})
                    .get('role_aliases', {})
                    .get('unrole', {})
                    .get(command_name)
            )
            if not alias_data:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No unrole alias configured for `{command_name}`.')
            static_channel_id = int(alias_data.get('channel_id')) if alias_data.get('channel_id') else None
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel for `{command_name}`.')
            is_owner_or_dev, is_coord = await check_owner_dev_coord(ctx, channel_obj)
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot unrole the bot.')
            if not is_coord and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to derole this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            static_role_id = int(alias_data.get('role_id'))
            role_obj = ctx.guild.get_role(static_role_id)
            if not role_obj:
                return await ctx.send(f'⚠️ Could not resolve role with ID `{static_role_id}`.')
            if role_obj not in member_obj.roles:
                return await ctx.send(f'{member_obj.mention} does not have {role_obj.mention}.')
            await member_obj.remove_roles(role_obj)
            return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} had {role_obj.mention} removed.', allowed_mentions=discord.AllowedMentions.none())
        return unrole_alias
        
    def create_untextmute_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Removes a text mute from a user in a specific text channel.')
        @is_owner_developer_coordinator_moderator_predicator('untmute')
        async def untext_mute_alias(
            ctx,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            static_channel_id = int(self.bot.command_aliases.get(ctx.guild.id, {}).get('channel_aliases', {}).get('untmute', {}).get(command_name))
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No text unmute alias configured for `{command_name}`.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                 return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel for `{command_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if member_obj.bot and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot undo a textmute on the bot.')
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to un textmute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow('''
                    SELECT expires_at FROM active_text_mutes
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                ''', ctx.guild.id, member_obj.id, static_channel_id)
                if row and row['expires_at'] is None and not await is_owner_developer_coordinator_via_alias(ctx, 'untmute'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Coordinator-only for undoing permanent text mutes.')
                try:
                    await channel_obj.set_permissions(member_obj, send_messages=None)
                except discord.Forbidden:
                    logger.warning('\U0001F6AB Discord forbidden: Cannot change the user\'s channel permissions.')
                    raise
                await conn.execute('DELETE FROM active_text_mutes WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3', ctx.guild.id, member_obj.id, static_channel_id)
                await conn.execute('INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1, $2, $3, $4, $5, $6)', 'untfmute', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Untextmuted a user')
            return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention}\'s text muted in {channel_obj.mention} has been removed.', allowed_mentions=discord.AllowedMentions.none())
        return untext_mute_alias

    @commands.command(name='xalias', help='Deletes an alias.')
    @is_owner_developer_coordinator_predicator(None)
    async def delete_alias(
        self,
        ctx: commands.Context,
        alias_name: Optional[str] = commands.parameter(default=None, description='Include an alias name')
    ) -> None:
        if not alias_name.strip():
            return await self.handler.send_message(ctx, content='\U0001F6AB `alias_name` cannot be empty.')
        guild_aliases = self.bot.command_aliases.get(ctx.guild.id, {})
        alias_type = None
        alias_dict = None
        for candidate in ('cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute'):
            if alias_name in guild_aliases.get('channel_aliases', {}).get(candidate, {}):
                alias_type = candidate
                alias_dict = 'channel_aliases'
                break
        if not alias_type:
            for candidate in ('role', 'unrole'):
                if alias_name in guild_aliases.get('role_aliases', {}).get(candidate, {}):
                    alias_type = candidate
                    alias_dict = 'role_aliases'
                    break
        if not alias_type:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Alias `{alias_name}` not found.')
        if alias_dict == 'channel_aliases':
            channel_id = guild_aliases['channel_aliases'][alias_type][alias_name]
            channel_obj = ctx.guild.get_channel(channel_id)
            if channel_obj:
                is_owner_or_dev, _ = await check_owner_dev_coord_mod(ctx, channel_obj)
                if not is_owner_or_dev:
                    async with ctx.bot.db_pool.acquire() as conn:
                        row = await conn.fetchrow('SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1', ctx.author.id)
                    if not row or channel_obj.id not in (row.get('coordinator_channel_ids') or []):
                        return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`xalias`) in {channel_obj.mention}.')
        else:
            role_id = guild_aliases['role_aliases'][alias_type][alias_name]
            channel = None
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM command_aliases WHERE guild_id = $1 AND alias_type = $2 AND alias_name = $3', ctx.guild.id, alias_type, alias_name)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1, $2, $3, $4, $5, $6)
            ''', 'delete_alias', None, ctx.author.id, ctx.guild.id, channel_obj.id if channel else None, f'Deleted alias {alias_name}')
        if self.bot.get_command(alias_name):
            self.bot.remove_command(alias_name)
        if alias_dict == 'channel_aliases':
            guild_aliases['channel_aliases'][alias_type].pop(alias_name, None)
        else:
            guild_aliases['role_aliases'][alias_type].pop(alias_name, None)
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Deleted alias `{alias_name}` from `{alias_type}`.')
            
    @commands.command(name='xcoord', help='Revokes coordinator access from a user in a specific voice channel.')
    @is_owner_developer_predicator()
    async def delete_coordinator(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ) -> None:
        channel_obj = await self.resolve_channel(ctx, channel)
        member_obj = await self.resolve_member(ctx, member)
        if not channel_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
            if not row:
                return await ctx.send(f'\U0001F6AB {member_obj.mention} is not found in the coordinator database.', allowed_mentions=discord.AllowedMentions.none())
            current_channel_ids = row.get('coordinator_channel_ids', []) or []
            if channel_obj.id not in current_channel_ids:
                return await ctx.send(f'\U0001F6AB {member_obj.mention} is not a coordinator in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            await conn.execute('''
                UPDATE users
                SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $2),
                    updated_at = NOW()
                WHERE discord_snowflake = $1
            ''', member_obj.id, channel_obj.id)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1, $2, $3, $4, $5, $6)
            ''', 'remove_coordinator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Removed a coordinator from a voice channel')
            updated_row = await conn.fetchrow('SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
            remaining_channels = updated_row.get('coordinator_channel_ids', []) if updated_row else []
            guild_voice_channel_ids = [vc.id for vc in ctx.guild.voice_channels]
            has_remaining_guild_channels = any(ch_id in guild_voice_channel_ids for ch_id in remaining_channels)
            if not has_remaining_guild_channels:
                return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention}\'s coordinator access has been completely revoked from {channel_obj.mention} and in {ctx.guild.name} (no remaining channels).', allowed_mentions=discord.AllowedMentions.none())
            else:
                return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention}\'s coordinator access has been revoked from {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='xdev', help='Removes a developer.')
    @is_owner_predicator()
    async def delete_developer(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
    ) -> None:
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET developer_guild_ids = array_remove(developer_guild_ids, $2),
                    updated_at = NOW()
                WHERE discord_snowflake = $1
            ''', member_obj.id, ctx.guild.id)
        return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention}\'s developer access has been revoked in {ctx.guild.name}.', allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='xmod', help='Revokes a member\'s VC moderator role for a given channel.')
    @is_owner_developer_coordinator_predicator(None)
    async def delete_moderator(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ) -> None:
        channel_obj = await self.resolve_channel(ctx, channel)
        member_obj = await self.resolve_member(ctx, member)
        if not channel_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`xmod`) in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT moderator_channel_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
            if not row:
                return await ctx.send(f'\U0001F6AB {member_obj.mention} is not found in the moderator database.', allowed_mentions=discord.AllowedMentions.none())
            current_channel_ids = row.get('moderator_channel_ids', []) or []
            if channel_obj.id not in current_channel_ids:
                return await ctx.send(f'\U0001F6AB {member_obj.mention} is not a moderator in {channel_obj.name}.', allowed_mentions=discord.AllowedMentions.none())
            await conn.execute('''
                UPDATE users
                SET moderator_channel_ids = array_remove(moderator_channel_ids, $2),
                    updated_at = NOW()
                WHERE discord_snowflake = $1
            ''', member_obj.id, channel_obj.id)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1, $2, $3, $4, $5, $6)
            ''', 'remove_moderator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Removed a moderator from the channel')
            updated_row = await conn.fetchrow('SELECT moderator_channel_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
        return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} has been revoked moderator access in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='bans', help='Lists ban statistics.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_bans(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        if not channel_obj and not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel or member from input: {target}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use command (`bans`) in {channel_obj.mention}.')
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners or developers can list all bans across the server.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT discord_snowflake, channel_id, expires_at, reason
                    FROM active_bans
                    WHERE guild_id = $1
                    ORDER BY channel_id, expires_at NULLS LAST
                ''', ctx.guild.id)
            if not rows:
                return await self.handler.send_message(ctx, content='\U0001F6AB No active bans found in {ctx.guild.name}.')
            grouped = defaultdict(list)
            for row in rows:
                grouped[row['channel_id']].append(row)
            embeds = []
            for ch_id, records in grouped.items():
                ch = ctx.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                embed = discord.Embed(title=f'⛔ Bans in {ch_name}', color=discord.Color.red())
                for record in records:
                    user = ctx.guild.get_member(record['discord_snowflake'])
                    reason = record['reason'] or 'No reason provided'
                    if record['expires_at'] is None:
                        duration_str = 'Permanent'
                    else:
                        now =  discord.utils.utcnow()
                        delta = record['expires_at'] - now
                        if delta.total_seconds() <= 0:
                            duration_str = 'Expired'
                        else:
                            days, seconds = delta.days, delta.seconds
                            hours = seconds // 3600
                            minutes = (seconds % 3600) // 60
                            duration_str = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                    mention = user.mention if user else f'`{record['discord_snowflake']}`'
                    embed.add_field(name='User', value=f'{mention}\nReason: {reason}\nDuration: {duration_str}', inline=False)
                embeds.append(embed)
            paginator = Paginator(self.bot, ctx, embeds)
            return await paginator.start()
        _, is_coord = await check_owner_dev_coord(ctx, channel_obj)
        if (is_owner_or_dev or is_coord) and member_obj:
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch('''
                    SELECT channel_id, expires_at, reason
                    FROM active_bans
                    WHERE guild_id = $1 AND discord_snowflake = $2
                ''', ctx.guild.id, member_obj.id)
            bans = [b for b in bans if ctx.guild.get_channel(b['channel_id'])]
            if not bans:
                return await ctx.send(f'\U0001F6AB {member_obj.mention} is not banned in any channels.', allowed_mentions=discord.AllowedMentions.none())
            embed = discord.Embed(title=f'Ban Records', description=f'For {member_obj.mention}', color=discord.Color.red())
            for record in bans:
                channel_obj = ctx.guild.get_channel(record['channel_id'])
                channel_mention = channel_obj.mention if channel_obj else f'Channel ID `{record['channel_id']}`'
                reason = record['reason'] or 'No reason provided'
                if record['expires_at'] is None:
                    duration_str = 'Permanent'
                else:
                    now =  discord.utils.utcnow()
                    delta = record['expires_at'] - now
                    if delta.total_seconds() <= 0:
                        duration_str = 'Expired'
                    else:
                        days, seconds = delta.days, delta.seconds
                        hours = seconds // 3600
                        minutes = (seconds % 3600) // 60
                        duration_str = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                embed.add_field(name=channel_mention, value=f'Reason: {reason}\nDuration: {duration_str}', inline=False)
            return await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.none())
        elif (is_mod_or_coord or is_owner_or_dev) and channel_obj:
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch('''
                    SELECT discord_snowflake, expires_at, reason
                    FROM active_bans
                    WHERE guild_id = $1 AND channel_id = $2
                    ORDER BY expires_at NULLS LAST
                ''', ctx.guild.id, channel_obj.id)
            if not bans:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No active bans found for {channel_obj.mention}.')
            lines = []
            for record in bans:
                uid = record['discord_snowflake']
                member_obj = ctx.guild.get_member(uid)
                if not member_obj:
                    continue
                name = member_obj.display_name
                if record['expires_at'] is None:
                    time_left = 'Permanent'
                else:
                    now =  discord.utils.utcnow()
                    delta = record['expires_at'] - now
                    if delta.total_seconds() <= 0:
                        time_left = 'Expired'
                    else:
                        days, seconds = delta.days, delta.seconds
                        hours = seconds // 3600
                        minutes = (seconds % 3600) // 60
                        time_left = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                lines.append(f'• {name} — {time_left} — <@{uid}>')
            if not lines:
                return await ctx.send(f'\U0001F6AB No active bans for users currently in {ctx.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(title=f'⛔ Active Bans in {channel_obj.mention}', description='\n'.join(chunk), color=discord.Color.red())
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F6AB You must specify a member, a text channel or use "all".')
    
    @commands.command(name='caps', help='List active caps for a channel or all channels if "all" is provided.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_caps(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        if not channel_obj and not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel or member from input: {target}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`caps`) {channel_obj.mention}.')
        lines, found_caps = [], False
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners or developers can list all caps.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch("SELECT channel_id, moderation_type, duration FROM active_caps WHERE guild_id=$1", ctx.guild.id)
            for row in rows:
                ch = ctx.guild.get_channel(row["channel_id"])
                ch_name = ch.mention if ch else f'Channel ID `{row["channel_id"]}`'
                lines.append(f'**{row["moderation_type"]} in {ch_name}** → `{row["duration"]}`')
                found_caps = True
            if not found_caps:
                return await self.handler.send_message(ctx, content='\U0001F6AB No caps found server-wide.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(title='All Active Caps in Server', description='\n'.join(chunk), color=discord.Color.red())
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        else:
            caps = await self.get_caps_for_channel(ctx.guild.id, channel_obj.id)
            for moderation_type, duration in caps:
                lines.append(f'**{moderation_type} in {channel_obj.mention}** → `{duration}`')
                found_caps = True
        if not found_caps:
            return await self.handler.send_message(ctx, content='\U0001F6AB No caps found for the specified channel or server-wide.')
        embed_title = 'All Active Caps in Server' if target and target.lower() == 'all' else f'Active Caps for {channel_obj.mention}'
        embed = discord.Embed(title=embed_title, description='\n'.join(lines), color=discord.Color.red())
        await self.handler.send_message(ctx, embed=embed)
        
    @commands.command(name='cmds', help='List command aliases routed to a specific channel or all channels if "all" is provided.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_commands(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        aliases = self.bot.command_aliases.get(ctx.guild.id, {})
        if not aliases:
            return await self.handler.send_message(ctx, content=f'No aliases defined in {ctx.guild.name}.')
        channel_obj = await self.resolve_channel(ctx, target)
        if not channel_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {target}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`cmds`) in {channel_obj.mention}.')
        lines = []
        found_aliases = False
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners or developers can list all aliases across the server.')
            for kind, type_map in aliases.get('channel_aliases', {}).items():
                grouped_by_channel = defaultdict(list)
                for name, cid in type_map.items():
                    grouped_by_channel[cid].append(name)
                for ch_id, names in grouped_by_channel.items():
                    ch = ctx.guild.get_channel(ch_id)
                    ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                    lines.append(f'**{kind.capitalize()} in {ch_name}**')
                    lines.extend(f'`{name}` → <#{ch_id}>' for name in names)
                    found_aliases = True
            for kind, type_map in aliases.get('role_aliases', {}).items():
                grouped_by_channel = defaultdict(list)
                for name, data in type_map.items():
                    if isinstance(data, dict) and 'channel_id' in data:
                        grouped_by_channel[data['channel_id']].append((name, data.get('role_id')))
                for ch_id, entries in grouped_by_channel.items():
                    ch = ctx.guild.get_channel(ch_id)
                    ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                    lines.append(f'**{kind.capitalize()} Role Aliases in {ch_name}**')
                    for name, rid in entries:
                        role = ctx.guild.get_role(rid)
                        mention = role.mention if role else f'<@&{rid}>'
                        lines.append(f'`{name}` → {mention}')
                    found_aliases = True
        else:
            channel_id = channel_obj.id
            for kind, type_map in aliases.get('channel_aliases', {}).items():
                channel_entries = {name: cid for name, cid in type_map.items() if cid == channel_id}
                if channel_entries:
                    found_aliases = True
                    lines.append(f'**{kind.capitalize()}**')
                    lines.extend(f'`{name}` → <#{cid}>' for name, cid in channel_entries.items())
            for kind, type_map in aliases.get('role_aliases', {}).items():
                role_entries = {name: data for name, data in type_map.items() if isinstance(data, dict) and data.get('channel_id') == channel_id}
                if role_entries:
                    found_aliases = True
                    lines.append(f'**{kind.capitalize()}**')
                    for name, data in role_entries.items():
                        rid = data.get('role_id')
                        role = ctx.guild.get_role(rid)
                        mention = role.mention if role else f'<@&{rid}>'
                        lines.append(f'`{name}` → {mention}')
        if not found_aliases:
            return await self.handler.send_message(ctx, content='\U0001F6AB No aliases found for the specified channel or server-wide.')
        embed_title = 'All Aliases in Server' if target and target.lower() == 'all' else f'Aliases for {channel_obj.mention}'
        embed = discord.Embed(
            title=embed_title,
            description='\n'.join(lines),
            color=discord.Color.blue()
        )
        await self.handler.send_message(ctx, embed=embed)

    @commands.command(name='coords', help='Lists coordinators for a specific voice channel, all, or a member.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_coordinators(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Voice channel name, mention, ID, "all", or member ID')
    ) -> None:
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        if not channel_obj and not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel or member from input: {target}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`coords`) in {channel_obj.mention}.')
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all coordinators.')
            query = '''
                SELECT unnest(coordinator_channel_ids) AS channel_id, discord_snowflake
                FROM users
                WHERE coordinator_channel_ids IS NOT NULL
            '''
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(query)
            if not rows:
                return await self.handler.send_message(ctx, content='\U0001F6AB No coordinators found in any voice channels.')
            channel_map = defaultdict(list)
            for row in rows:
                channel_map[row['channel_id']].append(row['discord_snowflake'])
            pages = []
            for ch_id, user_ids in sorted(channel_map.items()):
                vc = ctx.guild.get_channel(ch_id)
                vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(title=f'🧭 Coordinators for {vc_name}', color=discord.Color.gold())
                for uid in user_ids:
                    m = ctx.guild.get_member(uid)
                    name = m.display_name if m else f'User ID {uid}'
                    embed.add_field(name='\u200b', value=f'• {name} (<@{uid}>)', inline=False)
                pages.append(embed)
            if len(pages) == 1:
                return await ctx.send(embed=pages[0], allowed_mentions=discord.AllowedMentions.none())
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        if is_owner_or_dev and member_obj:
            if member_obj.id != ctx.author.id:
                if not is_owner_or_dev:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use (`coords`) for {member_obj.mention}.')
            query = '''
                SELECT coordinator_channel_ids
                FROM users
                WHERE discord_snowflake = $1
            '''
            async with self.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(query, member_obj.id)
            if not row or not row['coordinator_channel_ids']:
                return await ctx.send(f'\U0001F6AB {member_obj.display_name} is not a coordinator in any channels.')
            channel_mentions = []
            for ch_id in row['coordinator_channel_ids']:
                vc = ctx.guild.get_channel(ch_id)
                channel_mentions.append(vc.mention if vc else f'Unknown Channel ({ch_id})')
            chunk_size = 18
            embeds = []
            for i in range(0, len(channel_mentions), chunk_size):
                chunk = channel_mentions[i:i+chunk_size]
                embed = discord.Embed(
                    title=f'🧭 {member_obj.display_name} is a coordinator in:',
                    description='\n'.join(f'• {ch}' for ch in chunk),
                    color=discord.Color.gold()
                )
                embeds.append(embed)
            if len(embeds) == 1:
                return await ctx.send(embed=embeds[0], allowed_mentions=discord.AllowedMentions.none())
            paginator = Paginator(self.bot, ctx, embeds)
            return await paginator.start()
        elif (is_mod_or_coord or is_owner_or_dev) and channel_obj:
            query = '''
                SELECT discord_snowflake
                FROM users
                WHERE $1 = ANY (coordinator_channel_ids)
            '''
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(query, channel_obj.id)
            if not rows:
                return await ctx.send(
                    f'\U0001F6AB No coordinators found for {channel_obj.mention}.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            lines = []
            for row in rows:
                uid = row['discord_snowflake']
                m = ctx.guild.get_member(uid)
                if m:
                    lines.append(f'• {m.display_name} — <@{uid}>')
            if not lines:
                return await ctx.send(f'\U0001F6AB No coordinators currently in {ctx.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'🧭 Coordinators for {channel_obj.name}',
                    description='\n'.join(chunk),
                    color=discord.Color.gold()
                )
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        else:
            return await self.handler.send_message(ctx, content='\U0001F6AB You must specify a member, a voice channel, or use "all, or are not a moderator or above".')
    
    @commands.command(name='devs', hidden=True, help='Lists developers.')
    @is_owner_developer_predicator()
    async def list_developers(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Guild ID, "all", or user mention/ID')
    ) -> None:
        async with self.bot.db_pool.acquire() as conn:
            guild, pages = ctx.guild, []
            member_obj = await self.resolve_member(ctx, target)
            if target is None:
                rows = await conn.fetch('SELECT discord_snowflake FROM users WHERE $1 = ANY(developer_guild_ids)', guild.id)
                if not rows: return await self.handler.send_message(ctx, content=f'\U0001F6AB No developers are configured in {ctx.guild.name}.')
                for row in rows:
                    user = guild.get_member(row['discord_snowflake'])
                    name = user.display_name if user else f'User ID {row["discord_snowflake"]}'
                    embed = discord.Embed(title=f'Developer: {name}', color=discord.Color.blurple())
                    pages.append(embed)
            elif target.lower() == 'all':
                rows = await conn.fetch('SELECT discord_snowflake, developer_guild_ids FROM users WHERE array_length(developer_guild_ids, 1) > 0')
                if not rows: return await self.handler.send_message(ctx, content='\U0001F6AB No developers are configured.')
                for row in rows:
                    user = self.bot.get_user(row['discord_snowflake'])
                    name = user.name if user else f'User ID {row["discord_snowflake"]}'
                    guilds = [self.bot.get_guild(gid).name for gid in row['developer_guild_ids'] if self.bot.get_guild(gid)]
                    embed = discord.Embed(title=f'Developer: {name}', description=', '.join(guilds) if guilds else 'No known guilds', color=discord.Color.blurple())
                    pages.append(embed)
            else:
                if not member_obj:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {target}.')
                if member_obj.id != ctx.author.id:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`devs`) for {member_obj.mention}.')
                row = await conn.fetchrow('SELECT developer_guild_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
                if not row or not row['developer_guild_ids']: return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not a developer in any guilds.')
                guilds = [self.bot.get_guild(gid).name for gid in row['developer_guild_ids'] if self.bot.get_guild(gid)]
                embed = discord.Embed(title=f'Developer guilds for {member_obj.display_name}', description=', '.join(guilds) if guilds else 'No known guilds', color=discord.Color.blurple())
                pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        return await paginator.start()

    @commands.command(name='flags', help='List flag statistics.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_flags(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        if not channel_obj and not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel or member from input: {target}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`flags`) in {channel_obj.mention}.')
        guild = ctx.guild
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners or developers can list flags across all channels.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT channel_id, discord_snowflake
                    FROM active_flags
                    WHERE guild_id = $1
                ''', ctx.guild.id)
            if not rows:
                return await self.handler.send_message(ctx, content='\U0001F6AB No users are flagged in any voice channels.')
            channel_map = defaultdict(list)
            for row in rows:
                channel_map[row['channel_id']].append(row['discord_snowflake'])
            pages = []
            for ch_id, user_ids in sorted(channel_map.items()):
                ch = guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(title=f'🚩 Flagged Users in {ch_name}', color=discord.Color.yellow())
                for uid in user_ids:
                    m = guild.get_member(uid)
                    mention = m.mention if m else f'<@{uid}>'
                    embed.add_field(name='\u200b', value=f'• {mention}', inline=False)
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        _, is_coord = await check_owner_dev_coord(ctx, channel_obj)
        if (is_owner_or_dev or is_coord) and member_obj:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT channel_id, reason
                    FROM active_flags
                    WHERE guild_id = $1 AND discord_snowflake = $2
                ''', ctx.guild.id, member_obj.id)
                
                rows = [r for r in rows if guild.get_channel(r['channel_id'])]
                if not rows:
                    return await ctx.send(
                        f'\U0001F6AB {member_obj.mention} is not flagged in any voice channels.',
                        allowed_mentions=discord.AllowedMentions.none()
                    )
                lines = []
                for r in rows:
                    ch = guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                    reason = r['reason'] or "No reason given"
                    lines.append(f'• {ch_name} — {reason}')
                
                embed = discord.Embed(
                    title=f'🚩 Channels Where {member_obj.display_name} is Flagged',
                    description='\n'.join(lines),
                    color=discord.Color.orange()
                )
                return await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.all())
        elif (is_mod_or_coord or is_owner_or_dev) and channel_obj:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT discord_snowflake
                    FROM active_flags
                    WHERE guild_id = $1 AND channel_id = $2
                ''', ctx.guild.id, channel_obj.id)
            if not rows:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are flagged for {channel_obj.mention}.')
            pages = []
            chunk_size = 18
            for i in range(0, len(rows), chunk_size):
                chunk = rows[i:i + chunk_size]
                formatted_lines = []
                for record in chunk:
                    uid = record['discord_snowflake']
                    member_obj = guild.get_member(uid)
                    if not member_obj:
                        continue
                    formatted_lines.append(f'• {member_obj.display_name} — <@{uid}>')
                if formatted_lines:
                    embed = discord.Embed(
                        title=f'🚩 Flagged Users in {channel_obj.mention}',
                        color=discord.Color.red()
                    )
                    embed.add_field(name='\u200b', value='\n'.join(formatted_lines), inline=False)
                    pages.append(embed)
            if not pages:
                return await ctx.send(f'\U0001F6AB No flagged users currently in {guild.name}.')
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F6AB You must specify a member, a voice channel, or use "all".')

    @commands.command(name='logs', help='Lists log channels for this guild or a specified guild ID.')
    @is_owner_developer_predicator()
    async def list_logs(
        self,
        ctx: commands.Context,
        guild_id: Optional[int] = commands.parameter(default=None, description='Guild ID')
    ):
        guild_id = guild_id or (ctx.guild.id if ctx.guild else None)
        if not guild_id:
            await self.handler.send_message(ctx, content='\U0001F6AB No guild context or ID provided.')
            return
        entries = self.log_channels.get(guild_id, [])
        if not entries:
            await self.handler.send_message(ctx, content=f'\U0001F6AB No log channels configured in {ctx.guild.name}.')
            return
        embed = discord.Embed(
            title=f'{self.get_random_emoji()} Log Channels',
            color=discord.Color.blue()
        )
        embed.set_footer(text=f'Guild ID: {guild_id}')
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT * FROM log_channels WHERE guild_id=$1;', guild_id)
            for row in rows:
                channel_obj = self.bot.get_channel(row['channel_id'])
                mention = channel_obj.mention if channel_obj else f'`{row["channel_id"]}`'
                enabled = row.get('enabled', False)
                log_type = row.get('type') or 'general'
                snowflakes = row.get('snowflakes') or []
                if log_type == 'general':
                    detail = f"Logs all events in {ctx.guild.name}"
                elif log_type == 'channel':
                    detail = f"Logs only events from channels: {', '.join(f'<#{s}>' for s in snowflakes)}"
                elif log_type == 'member':
                    detail = f"Logs only events for members: {', '.join(f'<@{s}>' for s in snowflakes)}"
                else:
                    detail = "Unknown filter"
                embed.add_field(
                    name=f"{mention} {'✅' if enabled else '❌'}",
                    value=f"Type: **{log_type}**\n{detail}",
                    inline=False
                )
        await self.handler.send_message(ctx, embed=embed)

    @commands.command(name='mods', help='Lists moderator statistics.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_moderators(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Voice channel name/mention/ID, "all", or member mention/ID')
    ) -> None:
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        if not channel_obj and not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel or member from input: {target}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`mods`) in {channel_obj.mention}.')
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all moderators.')
            query = '''
                SELECT unnest(moderator_channel_ids) AS channel_id, discord_snowflake
                FROM users
                WHERE moderator_channel_ids IS NOT NULL
            '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    rows = await conn.fetch(query)
                if not rows:
                    return await self.handler.send_message(ctx, content='\U0001F6AB No moderators found in any voice channels.')
                channel_map = defaultdict(list)
                for row in rows:
                    channel_map[row['channel_id']].append(row['discord_snowflake'])
                pages = []
                for ch_id, user_ids in sorted(channel_map.items()):
                    vc = ctx.guild.get_channel(ch_id)
                    vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
                    embed = discord.Embed(title=f'🛡️ Moderators for {vc_name}', color=discord.Color.magenta())
                    for uid in user_ids:
                        m = ctx.guild.get_member(uid)
                        name = m.display_name if m else f'User ID {uid}'
                        embed.add_field(name='\u200b', value=f'• {name} (<@{uid}>)', inline=False)
                    pages.append(embed)
                if len(pages) == 1:
                    return await ctx.send(embed=pages[0], allowed_mentions=discord.AllowedMentions.none())
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            except Exception as e:
                await self.handler.send_message(ctx, content=f'Database error: {e}')
                raise
        _, is_coord = await check_owner_dev_coord(ctx, channel_obj)
        if is_owner_or_dev and member_obj:
            if member_obj.id != ctx.author.id:
                if not is_coord:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`mods`) for {member_obj.mention}.')
            query = '''
                SELECT moderator_channel_ids
                FROM users
                WHERE discord_snowflake = $1
            '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    row = await conn.fetchrow(query, member_obj.id)
                if not row or not row['moderator_channel_ids']:
                    return await ctx.send(f'\U0001F6AB {member_obj.display_name} is not a moderator in any channels.')
                channel_mentions = []
                for ch_id in row['moderator_channel_ids']:
                    vc = ctx.guild.get_channel(ch_id)
                    channel_mentions.append(vc.mention if vc else f'Unknown Channel ({ch_id})')
                chunk_size = 18
                embeds = []
                for i in range(0, len(channel_mentions), chunk_size):
                    chunk = channel_mentions[i:i+chunk_size]
                    embed = discord.Embed(
                        title=f'🛡️ {member_obj.display_name} moderates:',
                        description='\n'.join(f'• {ch}' for ch in chunk),
                        color=discord.Color.magenta()
                    )
                    embeds.append(embed)
                if len(embeds) == 1:
                    return await ctx.send(embed=embeds[0], allowed_mentions=discord.AllowedMentions.none())
                paginator = Paginator(self.bot, ctx, embeds)
                return await paginator.start()
            except Exception as e:
                logger.warning(f'\U0001F6AB Database error: {e}')
                raise
        elif (is_mod_or_coord or is_owner_or_dev) and channel_obj:
            query = '''
                SELECT discord_snowflake
                FROM users
                WHERE $1 = ANY (moderator_channel_ids)
            '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    rows = await conn.fetch(query, channel_obj.id)
                if not rows:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No moderators found for {channel_obj.mention}.')
                lines = []
                for row in rows:
                    uid = row['discord_snowflake']
                    m = ctx.guild.get_member(uid)
                    if not m:
                        continue
                    lines.append(f'• {m.display_name} — <@{uid}>')
                if not lines:
                    return await ctx.send(f'\U0001F6AB No moderators currently in {ctx.guild.name}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i + chunk_size]
                    embed = discord.Embed(
                        title=f'🛡️ Moderators for {channel_obj.name}',
                        description='\n'.join(chunk),
                        color=discord.Color.magenta()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            except Exception as e:
                logger.warning(f'\U0001F6AB Database error: {e}')
                raise
        return await self.handler.send_message(ctx, content='\U0001F6AB You must specify a member, a voice channel, or use "all".')

    @commands.command(name='mutes', help='Lists mute statistics.')
    @is_owner_developer_coordinator_predicator(None)
    async def list_mutes(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        if not channel_obj and not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel or member from input: {target}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`mutes`) in {channel_obj.mention}.')
        if target and target.lower() == 'all':
            if is_owner_or_dev:
                async with self.bot.db_pool.acquire() as conn:
                    records = await conn.fetch('''
                        SELECT discord_snowflake, channel_id, expires_at, COALESCE(reason, 'No reason provided') AS reason
                        FROM active_voice_mutes
                        WHERE guild_id = $1
                        ORDER BY channel_id, discord_snowflake
                    ''', ctx.guild.id)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No muted users currently in {ctx.guild.name}.')
                grouped = defaultdict(list)
                for record in records:
                    grouped[record['channel_id']].append(record)
                pages = []
                for channel_id, user_entries in sorted(grouped.items()):
                    channel = ctx.guild.get_channel(channel_id)
                    channel_name = channel_obj.mention if channel else f'Unknown Channel ({channel_id})'
                    chunk_size = 18
                    for i in range(0, len(user_entries), chunk_size):
                        embed = discord.Embed(title=f'🔇 Active Mutes in {channel_name}', color=discord.Color.orange())
                        for record in user_entries[i:i + chunk_size]:
                            user_id = record['discord_snowflake']
                            member = ctx.guild.get_member(user_id)
                            name = member_obj.display_name if member else f'User ID {user_id}'
                            mention = member_obj.mention if member else f'`{user_id}`'
                            reason = record['reason'] or 'No reason provided'
                            duration_str = self.fmt_duration(record['expires_at'])
                            embed.add_field(name=name, value=f'{mention}\nReason: {reason}\nDuration: {duration_str}', inline=False)
                        pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
        _, is_coord = await check_owner_dev_coord(ctx, channel_obj)
        if (is_owner_or_dev or is_coord) and member_obj:
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''
                    SELECT guild_id, channel_id, expires_at, reason
                    FROM active_voice_mutes
                    WHERE discord_snowflake = $1 AND guild_id = $2
                ''', member_obj.id, ctx.guild.id)
            records = [r for r in records if ctx.guild.get_channel(r['channel_id'])]
            if not records:
                return await ctx.send(f'\U0001F6AB {member_obj.mention} is not muted in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
            description_lines = []
            for record in records:
                channel_obj = ctx.guild.get_channel(record['channel_id'])
                channel_mention = channel_obj.mention if channel_obj else f'`{record['channel_id']}`'
                reason = record['reason']
                duration_str = self.fmt_duration(record['expires_at'])
                description_lines.append(f'• {channel_mention} — {reason} — {duration_str}')
            embed = discord.Embed(title=f'Mute Records for {member_obj.mention}', description='\n'.join(description_lines), color=discord.Color.orange())
            return await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.none())
        elif (is_mod_or_coord or is_owner_or_dev) and channel_obj:
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''
                    SELECT discord_snowflake, expires_at, COALESCE(reason, 'No reason provided') AS reason
                    FROM active_voice_mutes
                    WHERE channel_id = $1 AND guild_id = $2
                ''', channel_obj.id, ctx.guild.id)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are currently muted in {channel_obj.mention}.')
                description_lines = []
                for record in records:
                    uid = record['discord_snowflake']
                    member_obj = ctx.guild.get_member(uid)
                    if not member_obj:
                        continue
                    name = member_obj.display_name
                    duration_str = self.fmt_duration(record['expires_at'])
                    description_lines.append(f'• {name} — <@{uid}> — {duration_str}')
                if not description_lines:
                    return await ctx.send(f'\U0001F6AB No muted users currently in {ctx.guild.name}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(description_lines), chunk_size):
                    chunk = description_lines[i:i + chunk_size]
                    embed = discord.Embed(title=f'\U0001F507 Muted Users in {channel_obj.mention}', color=discord.Color.orange())
                    embed.add_field(name='Muted Users', value='\n'.join(chunk), inline=False)
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F6AB You must specify a member, a voice channel or be connected to a voice channel.')

    @commands.command(name='ls', help='List users cowed as going vegan in this guild.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_members(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ) -> None:
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        if not channel_obj and not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel or member from input: {target}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`ls`) in {channel_obj.mention}.')
        guild = ctx.guild
        try:
            async with self.bot.db_pool.acquire() as conn:
                if member_obj:
                    rows = await conn.fetch('''
                        SELECT channel_id
                        FROM active_cows
                        WHERE guild_id = $1 AND discord_snowflake = $2
                    ''', guild.id, member_obj.id)
                    if not rows:
                        return await ctx.send(f'\U0001F6AB {member_obj.mention} is not cowed in any channels.', allowed_mentions=discord.AllowedMentions.none())
                    lines = []
                    for r in rows:
                        ch = guild.get_channel(r['channel_id'])
                        ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                        lines.append(f'• {ch_name}')
                    embed = discord.Embed(
                        title=f'🐮 {member_obj.display_name}',
                        description='\n'.join(lines),
                        color=discord.Color.green()
                    )
                    return await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.all())
                elif channel_obj:
                    rows = await conn.fetch('''
                        SELECT discord_snowflake
                        FROM active_cows
                        WHERE guild_id = $1 AND channel_id = $2
                    ''', guild.id, channel_obj.id)
                    if not rows:
                        return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are cowed in {channel_obj.mention}.')
                    lines = []
                    for row in rows:
                        uid = row['discord_snowflake']
                        m = guild.get_member(uid)
                        if not m:
                            continue
                        lines.append(f'• {m.display_name} — <@{uid}>')
                    if not lines:
                        return await ctx.send(f'\U0001F6AB No new vegans currently in {guild.name}.')
                    chunk_size = 18
                    pages = []
                    for i in range(0, len(lines), chunk_size):
                        chunk = lines[i:i + chunk_size]
                        embed = discord.Embed(
                            title=f'🐮 New Vegans in {channel_obj.mention}',
                            description='\n'.join(chunk),
                            color=discord.Color.green()
                        )
                        pages.append(embed)
                    paginator = Paginator(self.bot, ctx, pages)
                    return await paginator.start()
        except Exception as e:
            await logger.warning(ctx, content=f'Database error: {e}')
            raise
            
    @commands.hybrid_command(name='rms', help='Lists all members with server mute privileges in this guild.')
    @is_owner_predicator()
    async def list_server_muters(
        self,
        ctx: commands.Context
    ) -> None:
        if not await is_owner(ctx, ctx.author.id):
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission because you are not an Owner.')
        async with self.bot.db_pool.acquire() as conn:
            records = await conn.fetch('''
                SELECT discord_snowflake
                FROM users
                WHERE $1 = ANY(server_muter_guild_ids)
                ORDER BY discord_snowflake
            ''', ctx.guild.id)
            if not records:
                return await ctx.send(f'\U0001F6AB No admins found in {ctx.guild.name}.', allowed_mentions=discord.AllowedMentions.none())
            description_lines = []
            for record in records:
                uid = record['discord_snowflake']
                member_obj = await self.resolve_member(ctx, uid)
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
            
    @commands.command(name='tmutes', help='Lists text-mute statistics.')
    @is_owner_developer_coordinator_predicator(None)
    async def list_text_mutes(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        if not channel_obj and not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel or member from input: {target}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`tmutes`) in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if target and target.lower() == 'all' and is_owner_or_dev:
                records = await conn.fetch('''
                    SELECT discord_snowflake, channel_id, reason, expires_at
                    FROM active_text_mutes
                    WHERE guild_id = $1
                    ORDER BY channel_id, discord_snowflake
                ''', ctx.guild.id)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are currently text-muted in {ctx.guild.name}.')
                grouped = defaultdict(list)
                for r in records: grouped[r['channel_id']].append(r)
                pages, chunk_size = [], 18
                for ch_id, entries in sorted(grouped.items()):
                    ch = ctx.guild.get_channel(ch_id)
                    ch_name = ch.mention if ch else f'Unknown Channel ({ch_id})'
                    for i in range(0, len(entries), chunk_size):
                        embed = discord.Embed(title=f'🔇 Text Mutes in {ch_name}', color=discord.Color.orange())
                        for e in entries[i:i + chunk_size]:
                            user = ctx.guild.get_member(e['discord_snowflake'])
                            mention = user.mention if user else f'`{e["discord_snowflake"]}`'
                            reason = e['reason'] or 'No reason provided'
                            duration_str = self.fmt_duration(e['expires_at'])
                            embed.add_field(name=mention, value=f'Reason: {reason}\nDuration: {duration_str}', inline=False)
                        pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            _, is_coord = await check_owner_dev_coord(ctx, channel_obj)
            if member_obj and (is_owner_or_dev or is_coord):
                records = await conn.fetch('''
                    SELECT channel_id, reason, expires_at
                    FROM active_text_mutes
                    WHERE discord_snowflake = $1 AND guild_id = $2
                ''', member_obj.id, ctx.guild.id)
                if not records:
                    return await ctx.send(f'\U0001F6AB {member_obj.mention} is not text-muted in any channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in records:
                    ch = ctx.guild.get_channel(r['channel_id'])
                    if not ch: continue
                    duration_str = self.fmt_duration(r['expires_at'])
                    lines.append(f'• {ch.mention} — {r["reason"]} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(title=f'Text Mute Records for {member_obj.mention}', description='\n'.join(lines[i:i+chunk_size]), color=discord.Color.orange())
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            if (is_owner_or_dev or is_mod_or_coord) and channel_obj:
                records = await conn.fetch('''
                    SELECT discord_snowflake, reason, expires_at
                    FROM active_text_mutes
                    WHERE channel_id = $1 AND guild_id = $2
                ''', channel_obj.id, ctx.guild.id)
                if not records:
                    return await ctx.send(f'\U0001F6AB No users are currently text-muted in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in records:
                    user = ctx.guild.get_member(r['discord_snowflake'])
                    if not user: continue
                    duration_str = self.fmt_duration(r['expires_at'])
                    lines.append(f'• {user.mention} — {r["reason"]} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(title=f'Text-Muted Users in {channel_obj.mention}', description='\n'.join(lines[i:i+chunk_size]), color=discord.Color.orange())
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F6AB You must specify "all", a member, or a text channel.')
        
    @commands.command(name='mlog', help='Create, modify, or delete a log channel.')
    @is_owner_developer_predicator()
    async def modify_log(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        action: Optional[str] = commands.parameter(default=None, description='create | modify | delete'),
        log_type: Optional[str] = commands.parameter(default=None, description='Type of logs: member, channel, etc'),
        *snowflakes: Optional[int]
    ):
        sf = [int(s) for s in snowflakes] if snowflakes else []
        channel_obj = await self.resolve_channel(ctx, channel)
        if not channel_obj or not channel:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        guild_id = ctx.guild.id
        current_entries = self.log_channels.setdefault(guild_id, [])
        async with self.bot.db_pool.acquire() as conn:
            if action.lower() == 'delete':
                await conn.execute('DELETE FROM log_channels WHERE guild_id=$1 AND channel_id=$2;', guild_id, channel_obj.id)
                current_entries[:] = [e for e in current_entries if e['channel_id'] != channel_obj.id]
                await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Log channel {channel_obj.mention} deleted.')
                return
            existing = await conn.fetchrow('SELECT * FROM log_channels WHERE guild_id=$1 AND channel_id=$2;', guild_id, channel_obj.id)
            if existing:
                await conn.execute('UPDATE log_channels SET type=$1, snowflakes=$2, enabled=TRUE WHERE guild_id=$3 AND channel_id=$4;', log_type, sf if sf else None, guild_id,     channel_obj.id)
                for e in current_entries:
                    if e['channel_id'] == channel_obj.id:
                        e.update({'type': log_type, 'snowflakes': sf if sf else None, 'enabled': True})
                msg = f'{self.get_random_emoji()} Log channel {channel_obj.mention} updated with type `{log_type or "general"}`.'
            else:
                await conn.execute('INSERT INTO log_channels (guild_id, channel_id, type, snowflakes, enabled) VALUES ($1, $2, $3, $4, TRUE);', guild_id, channel_obj.id, log_type, sf     if sf else None)
                current_entries.append({
                    'guild_id': guild_id,
                    'channel_id': channel_obj.id,
                    'type': log_type,
                    'snowflakes': sf if sf else None,
                    'enabled': True
                })
                msg = f'{self.get_random_emoji()} Log channel {channel_obj.mention} created with type `{log_type or "general"}`.'
        self.log_channels[guild_id] = current_entries
        await self.handler.send_message(ctx, content=msg)
        
    @commands.hybrid_command(name='rmute', help='Mutes the whole room (except yourself).')
    @is_owner_predicator()
    async def room_mute(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        duration: Optional[str] = commands.parameter(default='8', description='Duration of mute in hours. Example 0 (permanent), 30m, 2h, 5d'),
        *,
        reason: Optional[str] = commands.parameter(default='', description='Optional reason (required for 7 days or more)')
    ) -> None:
        channel_obj = await self.resolve_channel(ctx, channel)
        if not channel_obj or not channel:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        expires_at, duration_display = self.parse_duration(duration)
        is_owner_or_dev, is_coord = await check_owner_dev_coord(ctx, channel_obj)
        if (expires_at == '0' or expires_at is None) and (not is_coord or not reason.strip()):
            return await self.handler.send_message(ctx, content='\U0001F6AB Reason required and coordinator-only for permanent mutes.')
        skipped_members = []
        muted_members = []
        failed_members = []
        async with self.bot.db_pool.acquire() as conn:
            for member_obj in channel_obj.members:
                if await is_owner(ctx, member_obj.id):
                    skipped_members.append(member_obj)
                    continue
                try:
                    await conn.execute('''
                        INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, reason)
                        VALUES ($1, $2, $3, $4, $5)
                        ON CONFLICT (guild_id, discord_snowflake, channel_id) DO UPDATE
                        SET expires_at = EXCLUDED.expires_at,
                            reason = EXCLUDED.reason
                    ''', ctx.guild.id, member_obj.id, channel_obj.id, expires_at, reason or 'No reason provided')
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'voice_mute', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, reason or 'No reason provided')
                    if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                        await member_obj.edit(mute=True)
                    muted_members.append(member_obg)
                except Exception as e:
                    logger.warning(f'Failed to mute {member_obj.name}: {e}')
                    failed_members.append(member_obj)
        summary = f'{self.get_random_emoji()} Muted {len(muted_members)} member(s) in {channel_obj.mention} for {duration_display}.\nReason: {reason or 'No reason provided'}'
        if skipped_members:
            summary += f'\n Skipped {len(skipped_members)}.'
        if failed_members:
            summary += f'\n\U0001F6AB Failed to mute {len(failed_members)} member(s).'
        return await ctx.send(summary, allowed_mentions=discord.AllowedMentions.none())
    
    @commands.command(name='rmv', help='Move all the members in one room to another.')
    @is_owner_predicator()
    async def room_move_all(
        self,
        ctx: commands.Context,
        source_id: Optional[int] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        target_id: Optional[int] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        source_channel = ctx.guild.get_channel(source_id)
        target_channel = ctx.guild.get_channel(target_id)
        if not source_channel or not target_channel:
            await ctx.send('\U0001F6AB One or both channel IDs are invalid.')
            return
        await self.move_all_members(source_channel, target_channel)
        await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Moved all members from `{source_channel.name}` to `{target_channel.name}`.')
        
    @commands.command(name='runmute', help='Unmutes all members in a specific VC (except yourself).')
    @is_owner_predicator()
    async def room_unmute_all(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        *,
        reason: Optional[str] = commands.parameter(default='', description='Optional reason (required for 7 days or more)')
    ) -> None:
        channel_obj = await self.resolve_channel(ctx, channel)
        if not channel_obj or not channel:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        unmuted_members = []
        skipped_members = []
        failed_members = []
        async with self.bot.db_pool.acquire() as conn:
            for member_obj in channel_obj.members:
                if member_obj.id == ctx.author.id:
                    continue
                try:
                    row = await conn.fetchrow('''
                        SELECT guild_id FROM active_voice_mutes
                        WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                    ''', ctx.guild.id, member_obj.id, channel_obj.id)
                    if not row:
                        skipped_members.append(member_obj)
                        continue
                    if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                        await member_obj.edit(mute=False)
                    await conn.execute('''
                        DELETE FROM active_voice_mutes
                        WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                    ''', ctx.guild.id, member_obj.id, channel_obj.id)
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'unmute', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, f'Unmuted via unmute_all: {reason}')
                    unmuted_members.append(member_obj)
                except Exception as e:
                    logger.warning(f'Unmute failed for {member_obj.mention}: {e}')
                    failed_members.append(member_obj)
        summary = f'{self.get_random_emoji()} Unmuted {len(unmuted_members)} member(s) in {channel_obj.mention}.'
        if skipped_members:
            summary += f'\n\U000026A0\U0000FE0F Skipped {len(skipped_members)}.'
        if failed_members:
            summary += f'\n\U0001F6AB Failed to unmute {len(failed_members)}.'
        return await ctx.send(summary, allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='xadmin', help='Revokes server mute privileges from a user.')
    @is_owner_predicator()
    async def revoke_server_muter(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
    ):
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj or not member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot and ctx.author.id != int(ctx.bot.config['discord_owner_id']):
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot revoke server mute for the bot.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET server_muter_guild_ids = (SELECT ARRAY(
                    SELECT unnest(server_muter_guild_ids)
                     EXCEPT SELECT $2
                ))
                WHERE discord_snowflake = $1
            ''', member_obj.id, ctx.guild.id)
        self.server_muters.get(ctx.guild.id, set()).discard(member_obj.id)
        return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} no longer has server mute privileges.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='smute', help='Mutes a member throughout the entire guild.')
    @commands.check(lambda ctx: ctx.bot.get_cog('Hybrid').can_server_mute(ctx))
    async def server_mute(
            self,
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            *,
            reason: Optional[str] = commands.parameter(default='', description='Optional reason (required for 7 days or more)')
    ) -> None:
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj or not member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot and ctx.author.id != int(ctx.bot.config['discord_owner_id']):
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot server mute the bot.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (discord_snowflake, server_mute_guild_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[]) ON CONFLICT (discord_snowflake) DO
                UPDATE
                    SET server_mute_guild_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(COALESCE (u.server_mute_guild_ids, '{}') || ARRAY[$2])
                        )
                        FROM users u WHERE u.discord_snowflake = EXCLUDED.discord_snowflake
                    ),
                    updated_at = NOW()
            ''', member_obj.id, ctx.guild.id)
            await conn.execute('''
                INSERT INTO active_server_voice_mutes (guild_id, discord_snowflake, reason)
                VALUES ($1, $2, $3)
                ON CONFLICT (guild_id, discord_snowflake) DO UPDATE
                SET reason = EXCLUDED.reason
            ''', ctx.guild.id, member_obj.id, reason or 'No reason provided')
        if member_obj.voice and member_obj.voice.channel:
            await member_obj.edit(mute=True)
            return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} has been server muted in {channel_obj.mention} for reason: {reason}', allowed_mentions=discord.AllowedMentions.none())
        else:
            return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} has been server muted. They are not currently in a voice channel.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='log', help='Toggle logging for a channel on or off.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def toggle_log(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(description='Tag a channel or include its snowflake ID')
    ):
        channel_obj = await self.resolve_channel(ctx, channel)
        if not channel_obj or not channel:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        target_role, allowed = await check_block(ctx, ctx.author, channel_obj)
        if not allowed:
            await self.handler.send_message(ctx, content=f'❌ You must be {target_role} or higher to toggle logging for {channel_obj.mention}.')
            return
        guild_id = ctx.guild.id
        current_channels = self.log_channels.setdefault(guild_id, [])
        async with self.bot.db_pool.acquire() as conn:
            existing = await conn.fetchrow(
                'SELECT * FROM log_channels WHERE guild_id=$1 AND channel_id=$2;',
                guild_id, channel_obj.id
            )
            if existing:
                new_status = not existing['enabled']
                await conn.execute(
                    'UPDATE log_channels SET enabled=$1 WHERE guild_id=$2 AND channel_id=$3;',
                    new_status, guild_id, channel_obj.id
                )
                msg = f'{self.get_random_emoji()} Logging for {channel_obj.mention} toggled {"on" if new_status else "off"}.'
            else:
                await conn.execute(
                    'INSERT INTO log_channels (guild_id, channel_id, enabled) VALUES ($1, $2, TRUE);',
                    guild_id, channel_obj.id
                )
                current_channels.append(channel_obj.id)
                msg = f'{self.get_random_emoji()} Logging for {channel_obj.mention} enabled.'
        self.log_channels[guild_id] = current_channels
        await self.handler.send_message(ctx, content=msg)
        
    @commands.command(name='xsmute', help='Unmutes a member throughout the entire guild.')
    @commands.check(lambda ctx: ctx.bot.get_cog('Hybrid').can_server_mute(ctx))
    async def unsmute(
            self,
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(description='Tag a member or include their snowflake ID')
    ) -> None:
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj or not member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET server_mute_guild_ids = (SELECT ARRAY(
                    SELECT unnest(server_mute_guild_ids)
                    EXCEPT SELECT $2
                ))
                WHERE discord_snowflake = $1
            ''', member_obj.id, ctx.guild.id)
            await conn.execute('''
                DELETE FROM active_server_voice_mutes
                WHERE discord_snowflake = $1 AND guild_id = $2
            ''', member_obj.id, ctx.guild.id)
        if member_obj.voice and member_obj.voice.channel:
            await member_obj.edit(mute=False)
        return await ctx.send(f'{self.get_random_emoji()} {member_obj.mention} has been server unmuted.', allowed_mentions=discord.AllowedMentions.none())


    @staticmethod
    def can_server_mute(ctx):
        cog = ctx.bot.get_cog('Hybrid')
        return cog and ctx.author.id in cog.server_muters.get(ctx.guild.id, set())
        
    def create_ban_log_pages(
        self,
        ctx: commands.Context,
        member: discord.Member,
        channel: discord.VoiceChannel,
        duration_display: Optional[str],
        reason: Optional[str],
        executor: discord.Member,
        expires_at: Optional[datetime],
        command_used: Optional[str],
        was_in_channel: bool = False,
        is_modification: bool = False,
        guild: discord.Guild = None,
        highest_role: Optional[str] = ''
    ):
        # Duration-based styling
        if expires_at is None:
            color = 0xDC143C
            ban_type = '🔒 Permanent'
            duration_emoji = '♾️'
        elif (expires_at - datetime.now(timezone.utc)).days >= 7:
            color = 0xFF6B35
            ban_type = '⏰ Extended'
            duration_emoji = '📅'
        else:
            color = 0xFF8C00
            ban_type = '⏱️ Temporary'
            duration_emoji = '⏰'

        title = '🔄 Ban Modified' if is_modification else '🔨 User Banned'
        guild = guild or ctx.guild
        
        # PRIORITY EMBED 1: User Identity & Images (Highest Priority)
        embed_user = discord.Embed(title=f"{title} - User Identity", color=color, timestamp=datetime.now(timezone.utc))
        embed_user.description = f"**Target:** {member.mention} banned from {channel.mention}"
        
        # High priority user info - Discord IGNs and snowflakes
        user_priority = f"**Display Name:** {member.display_name}\n"
        user_priority += f"**Username:** @{member.name}\n"
        user_priority += f"**User ID:** `{member.id}`\n"
        user_priority += f"**Account Age:** <t:{int(member.created_at.timestamp())}:R>"
        if member.joined_at:
            user_priority += f"\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
        
        embed_user.add_field(name='👤 Target User (Priority Info)', value=user_priority, inline=False)
        
        # Set user avatar as main image (high priority)
        embed_user.set_image(url=member.display_avatar.url)
        embed_user.set_thumbnail(url=executor.display_avatar.url)
        
        # Executor priority info
        exec_priority = f"**Moderator:** {executor.display_name} (@{executor.name})\n"
        exec_priority += f"**Mod ID:** `{executor.id}`\n"
        exec_priority += f"**Top Role:** {highest_role or executor.top_role.mention}"
        embed_user.add_field(name='👮‍♂️ Executed By', value=exec_priority, inline=True)
        
        # Command context reference
        ctx_info = f"**Original Message ID:** `{ctx.message.id}`\n"
        ctx_info += f"**Command Channel:** {ctx.channel.mention}\n"
        ctx_info += f"**Command Used:** `{command_used}`"
        embed_user.add_field(name='📱 Command Context', value=ctx_info, inline=True)
        
        embed_user.set_footer(
            text=f"Ban Ref: {member.id}-{channel.id} | Msg: {ctx.message.id}",
            icon_url=guild.icon.url if guild and guild.icon else None
        )

        # PRIORITY EMBED 2: Duration & Action Details
        embed_duration = discord.Embed(title=f"{title} - Duration Info", color=color, timestamp=datetime.now(timezone.utc))
        
        if expires_at:
            time_left = expires_at - datetime.now(timezone.utc)
            hours_left = round(time_left.total_seconds() / 3600, 1)
            days_left = time_left.days
            duration_info = f'**Type:** {ban_type}\n**Duration:** {duration_display}\n**Expires:** <t:{int(expires_at.timestamp())}:F>\n**Time Left:** '
            duration_info += f'{days_left}d {hours_left % 24:.1f}h' if days_left > 0 else f'{hours_left}h'
        else:
            duration_info = f'**Type:** {ban_type}\n**Duration:** {duration_display}\n**Status:** Permanent Ban'
        
        embed_duration.add_field(name=f'{duration_emoji} Ban Duration', value=duration_info, inline=False)
        
        # Action details
        action_details = f"**Was in Channel:** {'✅ Yes' if was_in_channel else '❌ No'}\n"
        action_details += f"**Action Type:** {'Modification' if is_modification else 'New Ban'}\n"
        action_details += f"**Server:** {guild.name} (`{guild.id}`)"
        embed_duration.add_field(name='⚙️ Action Details', value=action_details, inline=True)
        
        # Minimal channel info (lower priority)
        channel_basic = f"**Channel:** {channel.mention} (`{channel.id}`)\n**Category:** {channel.category.name if channel.category else 'None'}"
        embed_duration.add_field(name='📍 Channel Info', value=channel_basic, inline=True)

        # Start with priority embeds
        embeds = [embed_user, embed_duration]
        
        # REASON EMBEDS: Split long reasons
        reason_chunks = [reason[i:i+1000] for i in range(0, len(reason), 1000)] if reason else ['No reason provided']
        for i, chunk in enumerate(reason_chunks):
            reason_embed = discord.Embed(title=f"{title} - Reason", color=color, timestamp=datetime.now(timezone.utc))
            reason_embed.add_field(name=f'📝 Reason (Part {i+1})', value=f'```{chunk}```', inline=False)
            embeds.append(reason_embed)
        
        return embeds


    def create_text_mute_log_pages(
            self,
            ctx: commands.Context,
            member: discord.Member,
            channel: discord.VoiceChannel,
            duration_display: Optional[str],
            reason: Optional[str],
            executor: discord.Member,
            expires_at: Optional[datetime],
            command_used: Optional[str],
            was_in_channel: bool = False,
            is_modification: bool = False,
            guild: discord.Guild = None,
            highest_role: Optional[str] = ''
        ):
            # Duration-based styling
            if expires_at is None:
                color = 0xDC143C
                mute_type = '🔒 Permanent'
                duration_emoji = '♾️'
            elif (expires_at - datetime.now(timezone.utc)).days >= 7:
                color = 0xFF6B35
                mute_type = '⏰ Extended'
                duration_emoji = '📅'
            else:
                color = 0xFF8C00
                mute_type = '⏱️ Temporary'
                duration_emoji = '⏰'
    
            title = '🔄 Text Mute Modified' if is_modification else '🔨 User Text Muted'
            guild = guild or ctx.guild
    
            # PRIORITY EMBED 1: User Identity & Images
            embed_user = discord.Embed(title=f"{title} - User Identity", color=color, timestamp=datetime.now(timezone.utc))
            embed_user.description = f"**Target:** {member.mention} text muted in {channel.mention}"
    
            user_priority = f"**Display Name:** {member.display_name}\n"
            user_priority += f"**Username:** @{member.name}\n"
            user_priority += f"**User ID:** `{member.id}`\n"
            user_priority += f"**Account Age:** <t:{int(member.created_at.timestamp())}:R>"
            if member.joined_at:
                user_priority += f"\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
    
            embed_user.add_field(name='👤 Target User (Priority Info)', value=user_priority, inline=False)
            embed_user.set_image(url=member.display_avatar.url)
            embed_user.set_thumbnail(url=executor.display_avatar.url)
    
            exec_priority = f"**Moderator:** {executor.display_name} (@{executor.name})\n"
            exec_priority += f"**Mod ID:** `{executor.id}`\n"
            exec_priority += f"**Top Role:** {highest_role or executor.top_role.mention}"
            embed_user.add_field(name='👮‍♂️ Executed By', value=exec_priority, inline=True)
    
            ctx_info = f"**Original Message ID:** `{ctx.message.id}`\n"
            ctx_info += f"**Command Channel:** {ctx.channel.mention}\n"
            ctx_info += f"**Command Used:** `{command_used}`"
            embed_user.add_field(name='📱 Command Context', value=ctx_info, inline=True)
    
            embed_user.set_footer(
                text=f"Text Mute Ref: {member.id}-{channel.id} | Msg: {ctx.message.id}",
                icon_url=guild.icon.url if guild and guild.icon else None
            )
    
            # PRIORITY EMBED 2: Duration & Action Details
            embed_duration = discord.Embed(title=f"{title} - Duration Info", color=color, timestamp=datetime.now(timezone.utc))
    
            if expires_at:
                time_left = expires_at - datetime.now(timezone.utc)
                hours_left = round(time_left.total_seconds() / 3600, 1)
                days_left = time_left.days
                duration_info = f'**Type:** {mute_type}\n**Duration:** {duration_display}\n**Expires:** <t:{int(expires_at.timestamp())}:F>\n**Time Left:** '
                duration_info += f'{days_left}d {hours_left % 24:.1f}h' if days_left > 0 else f'{hours_left}h'
            else:
                duration_info = f'**Type:** {mute_type}\n**Duration:** {duration_display}\n**Status:** Permanent Text Mute'
    
            embed_duration.add_field(name=f'{duration_emoji} Mute Duration', value=duration_info, inline=False)
    
            action_details = f"**Was in Channel:** {'✅ Yes' if was_in_channel else '❌ No'}\n"
            action_details += f"**Action Type:** {'Modification' if is_modification else 'New Text Mute'}\n"
            action_details += f"**Server:** {guild.name} (`{guild.id}`)"
            embed_duration.add_field(name='⚙️ Action Details', value=action_details, inline=True)
    
            channel_basic = f"**Channel:** {channel.mention} (`{channel.id}`)\n**Category:** {channel.category.name if channel.category else 'None'}"
            embed_duration.add_field(name='📍 Channel Info', value=channel_basic, inline=True)
    
            embeds = [embed_user, embed_duration]
    
            # REASON EMBEDS
            reason_chunks = [reason[i:i+1000] for i in range(0, len(reason), 1000)] if reason else ['No reason provided']
            for i, chunk in enumerate(reason_chunks):
                reason_embed = discord.Embed(title=f"{title} - Reason", color=color, timestamp=datetime.now(timezone.utc))
                reason_embed.add_field(name=f'📝 Reason (Part {i+1})', value=f'```{chunk}```', inline=False)
                embeds.append(reason_embed)
    
            return embeds
    
    
    def create_voice_mute_log_pages(
            self,
            ctx: commands.Context,
            member: discord.Member,
            channel: discord.VoiceChannel,
            duration_display: Optional[str],
            reason: Optional[str],
            executor: discord.Member,
            expires_at: Optional[datetime],
            command_used: Optional[str],
            was_in_channel: bool = False,
            is_modification: bool = False,
            guild: discord.Guild = None,
            highest_role: Optional[str] = ''
        ):
            # Duration-based styling
            if expires_at is None:
                color = 0xDC143C
                mute_type = '🔒 Permanent'
                duration_emoji = '♾️'
            elif (expires_at - datetime.now(timezone.utc)).days >= 7:
                color = 0xFF6B35
                mute_type = '⏰ Extended'
                duration_emoji = '📅'
            else:
                color = 0xFF8C00
                mute_type = '⏱️ Temporary'
                duration_emoji = '⏰'
    
            title = '🔄 Voice Mute Modified' if is_modification else '🔨 User Voice Muted'
            guild = guild or ctx.guild
    
            # PRIORITY EMBED 1: User Identity & Images
            embed_user = discord.Embed(title=f"{title} - User Identity", color=color, timestamp=datetime.now(timezone.utc))
            embed_user.description = f"**Target:** {member.mention} voice muted in {channel.mention}"
    
            user_priority = f"**Display Name:** {member.display_name}\n"
            user_priority += f"**Username:** @{member.name}\n"
            user_priority += f"**User ID:** `{member.id}`\n"
            user_priority += f"**Account Age:** <t:{int(member.created_at.timestamp())}:R>"
            if member.joined_at:
                user_priority += f"\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
    
            embed_user.add_field(name='👤 Target User (Priority Info)', value=user_priority, inline=False)
            embed_user.set_image(url=member.display_avatar.url)
            embed_user.set_thumbnail(url=executor.display_avatar.url)
    
            exec_priority = f"**Moderator:** {executor.display_name} (@{executor.name})\n"
            exec_priority += f"**Mod ID:** `{executor.id}`\n"
            exec_priority += f"**Top Role:** {highest_role or executor.top_role.mention}"
            embed_user.add_field(name='👮‍♂️ Executed By', value=exec_priority, inline=True)
    
            ctx_info = f"**Original Message ID:** `{ctx.message.id}`\n"
            ctx_info += f"**Command Channel:** {ctx.channel.mention}\n"
            ctx_info += f"**Command Used:** `{command_used}`"
            embed_user.add_field(name='📱 Command Context', value=ctx_info, inline=True)
    
            embed_user.set_footer(
                text=f"Voice Mute Ref: {member.id}-{channel.id} | Msg: {ctx.message.id}",
                icon_url=guild.icon.url if guild and guild.icon else None
            )
    
            # PRIORITY EMBED 2: Duration & Action Details
            embed_duration = discord.Embed(title=f"{title} - Duration Info", color=color, timestamp=datetime.now(timezone.utc))
    
            if expires_at:
                time_left = expires_at - datetime.now(timezone.utc)
                hours_left = round(time_left.total_seconds() / 3600, 1)
                days_left = time_left.days
                duration_info = f'**Type:** {mute_type}\n**Duration:** {duration_display}\n**Expires:** <t:{int(expires_at.timestamp())}:F>\n**Time Left:** '
                duration_info += f'{days_left}d {hours_left % 24:.1f}h' if days_left > 0 else f'{hours_left}h'
            else:
                duration_info = f'**Type:** {mute_type}\n**Duration:** {duration_display}\n**Status:** Permanent Voice Mute'
    
            embed_duration.add_field(name=f'{duration_emoji} Mute Duration', value=duration_info, inline=False)
    
            action_details = f"**Was in Channel:** {'✅ Yes' if was_in_channel else '❌ No'}\n"
            action_details += f"**Action Type:** {'Modification' if is_modification else 'New Voice Mute'}\n"
            action_details += f"**Server:** {guild.name} (`{guild.id}`)"
            embed_duration.add_field(name='⚙️ Action Details', value=action_details, inline=True)
    
            channel_basic = f"**Channel:** {channel.mention} (`{channel.id}`)\n**Category:** {channel.category.name if channel.category else 'None'}"
            embed_duration.add_field(name='📍 Channel Info', value=channel_basic, inline=True)
    
            embeds = [embed_user, embed_duration]
    
            # REASON EMBEDS
            reason_chunks = [reason[i:i+1000] for i in range(0, len(reason), 1000)] if reason else ['No reason provided']
            for i, chunk in enumerate(reason_chunks):
                reason_embed = discord.Embed(title=f"{title} - Reason", color=color, timestamp=datetime.now(timezone.utc))
                reason_embed.add_field(name=f'📝 Reason (Part {i+1})', value=f'```{chunk}```', inline=False)
                embeds.append(reason_embed)
    
            return embeds
    
    def fmt_duration(self, expires_at: Optional[datetime], base: Optional[datetime] = None) -> str:
        if not expires_at:
            return "permanently"
        now = discord.utils.utcnow()
        delta = expires_at - (base or now)
        total_seconds = int(delta.total_seconds())
        if total_seconds <= 0:
            return "expired"
        days, remainder = divmod(total_seconds, 86400)
        hours, remainder = divmod(remainder, 3600)
        minutes, _ = divmod(remainder, 60)
        if base:
            if expires_at > base:
                prefix = "for "
            else:
                prefix = "reduced by "
        else:
            prefix = "for "
        if days > 0:
            return f"{prefix}{days} day(s)"
        elif hours > 0:
            return f"{prefix}{hours} hour(s)"
        elif minutes > 0:
            return f"{prefix}{minutes} minute(s)"
        else:
            return f"{prefix}less than a minute"
        
    async def get_cap(self, channel_id: int, guild_id: int, moderation_type: Optional[str]) -> Optional[str]:
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3',
                guild_id, channel_id, moderation_type
            )
            return row['duration'] if row else None
            
    async def get_caps_for_channel(self, guild_id: int, channel_id: int) -> list[tuple[str,str]]:
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT moderation_type, duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2',
                guild_id, channel_id
            )
            return [(r['moderation_type'], r['duration']) for r in rows]
    
    async def resolve_member(self, ctx: commands.Context, value: Optional[Union[int, str, discord.Member]]) -> Optional[discord.Member]:
        try:
            if isinstance(value, discord.Member):
                logger.debug(f"Direct member: {value.id}")
                return value
            if isinstance(value, int):
                m = ctx.guild.get_member(value)
                if not m:
                    try: m = await ctx.guild.fetch_member(value)
                    except discord.NotFound: m = None
                if m:
                    logger.debug(f"Resolved member by int ID: {m.id}")
                    return m
            if isinstance(value, str):
                if value.isdigit():
                    mid = int(value)
                    m = ctx.guild.get_member(mid)
                    if not m:
                        try: m = await ctx.guild.fetch_member(mid)
                        except discord.NotFound: m = None
                    if m:
                        logger.debug(f"Resolved member by str ID: {m.id}")
                        return m
                if value.startswith('<@') and value.endswith('>'):
                    mid = int(value[2:-1].replace('!', ''))
                    m = ctx.guild.get_member(mid)
                    if not m:
                        try: m = await ctx.guild.fetch_member(mid)
                        except discord.NotFound: m = None
                    if m:
                        logger.debug(f"Member mention resolved: {m.id}")
                        return m
                    return m
        except Exception as e:
            logger.warning(f"Member resolution error: {e}")
        return None
    
    async def resolve_channel(self, ctx: commands.Context, value: Optional[Union[int, str, discord.TextChannel, discord.VoiceChannel]]) -> Optional[Union[discord.TextChannel,     discord.VoiceChannel]]:
        try:
            if isinstance(value, (discord.TextChannel, discord.VoiceChannel)):
                logger.debug(f"Direct channel: {value.id}")
                return value
            if isinstance(value, int):
                c = ctx.guild.get_channel(value)
                if isinstance(c, (discord.TextChannel, discord.VoiceChannel)):
                    logger.debug(f"Resolved channel by int ID: {c.id}")
                    return c
            if isinstance(value, str):
                if value.isdigit():
                    cid = int(value)
                    c = ctx.guild.get_channel(cid)
                    if isinstance(c, (discord.TextChannel, discord.VoiceChannel)):
                        logger.debug(f"Resolved channel by str ID: {c.id}")
                        return c
                if value.startswith('<#') and value.endswith('>'):
                    cid = int(value[2:-1])
                    c = ctx.guild.get_channel(cid)
                    if c:
                        logger.debug(f"Channel mention resolved: {c.id}")
                        return c
        except Exception as e:
            logger.warning(f"Channel resolution error: {e}")
        return ctx.channel

    async def get_highest_role(self, ctx: commands.Context, member_obj: discord.Member, channel_obj: discord.abc.GuildChannel) -> str:
        bot = ctx.bot
        role_hierarchy = ['Everyone', 'Moderator', 'Coordinator', 'Developer', 'Owner']
        member_roles = []
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT coordinator_channel_ids, moderator_channel_ids, developer_guild_ids FROM users WHERE discord_snowflake = $1',
                member_obj.id
            )
        if row:
            if row.get('moderator_channel_ids') and channel_obj.id in row['moderator_channel_ids']:
                member_roles.append('Moderator')
            if row.get('coordinator_channel_ids') and channel_obj.id in row['coordinator_channel_ids']:
                member_roles.append('Coordinator')
            if row.get('developer_guild_ids') and ctx.guild and ctx.guild.id in row['developer_guild_ids']:
                member_roles.append('Developer')
        if ctx.guild and member_obj.id == ctx.guild.owner_id:
            member_roles.append('Owner')
        if member_obj.id == int(bot.config['discord_owner_id']):
            member_roles.append('Owner')
        return max(member_roles, key=lambda r: role_hierarchy.index(r)) if member_roles else 'Everyone'

    def get_random_emoji(self):
        return random.choice(VEGAN_EMOJIS)
    
    async def load_server_muters(self):
        await self.bot.wait_until_ready()
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT discord_snowflake, server_muter_guild_ids FROM users WHERE server_muter_guild_ids IS NOT NULL')
            for row in rows:
                user_id = row['discord_snowflake']
                for guild_id in row['server_muter_guild_ids']:
                    self.server_muters[guild_id].add(user_id)
                    
    async def move_all_members(self, source_channel: discord.VoiceChannel, target_channel: discord.VoiceChannel):
        if not isinstance(source_channel, discord.VoiceChannel) or not isinstance(target_channel, discord.VoiceChannel):
            raise ValueError('\U0001F6AB Both source and target must be voice channels.')
        for member in source_channel.members:
            try:
                await member.move_to(target_channel)
            except discord.Forbidden:
                logger.warning(f'\U0001F6AB Missing permissions to move {member}.')
            except discord.HTTPException:
                logger.warning(f'\U0001F6AB Failed to move {member} due to a network error.')
                
                
    def parse_duration(self, duration: Optional[str], base: Optional[datetime] = None) -> tuple[Optional[datetime], str]:
        if duration is None:
            delta = timedelta(hours=24)
            return (datetime.now(timezone.utc) + delta), 'for 24 hour(s)'
        duration = duration.strip().lower()
        sign = 1
        is_relative = False
        if duration.startswith('+'):
            is_relative = True
            duration = duration[1:]
        elif duration.startswith('-'):
            is_relative = True
            sign = -1
            duration = duration[1:]
        if duration in ('0', '0h', '0hr', '0hrs', '0hour', '0hours', '0m', '0min', '0mins', '0minute', '0minutes', '0d', '0day', '0days'):
            return None, 'permanently'
        if duration.endswith(('d', 'day', 'days')):
            value = int(duration.rstrip('daysy'))
            delta = timedelta(days=value * sign)
            target = (base if (is_relative and base) else datetime.now(timezone.utc)) + delta
            return target, f'for {value} day(s)' if sign > 0 else f'reduced by {value} day(s)'
        if duration.endswith(('h', 'hr', 'hrs', 'hour', 'hours')):
            value = int(duration.rstrip('hrshours'))
            delta = timedelta(hours=value * sign)
            target = (base if (is_relative and base) else datetime.now(timezone.utc)) + delta
            return target, f'for {value} hour(s)' if sign > 0 else f'reduced by {value} hour(s)'
        if duration.endswith(('m', 'min', 'mins', 'minute', 'minutes')):
            value = int(duration.rstrip('minsmutes'))
            delta = timedelta(minutes=value * sign)
            target = (base if (is_relative and base) else datetime.now(timezone.utc)) + delta
            return target, f'for {value} minute(s)' if sign > 0 else f'reduced by {value} minute(s)'
        value = int(duration)
        delta = timedelta(hours=value * sign)
        target = (base if (is_relative and base) else datetime.now(timezone.utc)) + delta
        return target, f'for {value} hour(s)' if sign > 0 else f'reduced by {value} hour(s)'

    def perform_backup(self, db_user: Optional[str], db_name: Optional[str], db_host: Optional[str], db_password: Optional[str], backup_dir: Optional[str]) -> str:
        timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        backup_file = os.path.join(backup_dir, f'backup_{timestamp}.sql')
        dump_command = [
            'pg_dump',
            '-U', db_user,
            '-h', db_host,
            '-d', db_name,
            '-F', 'p',
            '-f', backup_file,
        ]
        env = os.environ.copy()
        env['PGPASSWORD'] = db_password
        result = subprocess.run(
            dump_command,
            capture_output=True,
            text=True,
            env=env,
        )
        if result.returncode != 0:
            raise RuntimeError(f'Backup failed: {result.stderr}')
        return backup_file
        
    async def send_log(
        self,
        ctx: commands.Context,
        moderation_type: Optional[str],
        member: discord.Member,
        channel: Optional[discord.VoiceChannel],
        duration_display: Optional[str],
        reason: Optional[str],
        executor: discord.Member,
        expires_at: Optional[datetime],
        command_used: Optional[str],
        was_in_channel: bool = False,
        is_modification: bool = False,
        highest_role: Optional[str] = 'Everyone'
    ):
        guild_id = ctx.guild.id
        guild = ctx.guild
        if guild_id not in self.log_channels:
            print(f"No log channels for guild {guild_id}")
            return
        for entry in self.log_channels[guild_id]:
            log_channel = self.bot.get_channel(entry["channel_id"])
            if not log_channel:
                continue
            log_type = entry.get("type")
            snowflakes = entry.get("snowflakes") or []
            if log_type == "member" and member.id not in snowflakes:
                continue
            if log_type == "channel" and channel and channel.id not in snowflakes:
                continue
            if moderation_type == "ban":
                pages = self.create_ban_log_pages(
                    ctx=ctx, member=member, channel=channel, duration_display=duration_display,
                    reason=reason, executor=executor, expires_at=expires_at,
                    command_used=command_used, was_in_channel=was_in_channel,
                    is_modification=is_modification, guild=guild, highest_role=highest_role
                )
            elif moderation_type == "tmute":
                pages = self.create_text_mute_log_pages(
                    ctx=ctx, member=member, channel=channel, duration_display=duration_display,
                    reason=reason, executor=executor, expires_at=expires_at,
                    command_used=command_used, was_in_channel=was_in_channel,
                    is_modification=is_modification, guild=guild, highest_role=highest_role
                )
            elif moderation_type == "vmute":
                pages = self.create_voice_mute_log_pages(
                    ctx=ctx, member=member, channel=channel, duration_display=duration_display,
                    reason=reason, executor=executor, expires_at=expires_at,
                    command_used=command_used, guild=guild, was_in_channel=was_in_channel,
                    is_modification=is_modification, highest_role=highest_role
                )
            else:
                continue
            paginator = ChannelPaginator(self.bot, log_channel, pages)
            await paginator.start()

    async def set_cap(self, channel_id: int, guild_id: int, moderation_type: Optional[str], duration: Optional[str]):
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                '''INSERT INTO active_caps (guild_id, channel_id, moderation_type, duration)
                   VALUES ($1,$2,$3,$4)
                   ON CONFLICT (guild_id, channel_id, moderation_type)
                   DO UPDATE SET duration=EXCLUDED.duration''',
                guild_id, channel_id, moderation_type, duration
            )

    def setup_backup_directory(self, backup_dir: Optional[str]) -> str:
        os.makedirs(backup_dir, exist_ok=True)
        return backup_dir

    async def vegan_db(self, ctx: commands.Context) -> bool:
        channel = ctx.channel
        user_id = ctx.author.id
        channel_id = channel.id
        guild_id = ctx.guild.id
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT coordinator_channel_ids
                FROM users
                WHERE discord_snowflake = $1
            ''', user_id)
        if not row:
            await ctx.send('\U0001F6AB You are not registered in the database.')
            return
        coordinator_channel_ids = row['coordinator_channel_ids'] or []
        is_coordinator = channel_id in coordinator_channel_ids
        is_vegan_channel = 'vegan' in channel.name.lower()
        if is_coordinator and is_vegan_channel:
            return True
        else:
            return False
        
async def setup(bot: DiscordBot):
    cog = Hybrid(bot)
    await bot.add_cog(cog)
