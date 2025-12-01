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
from vyrtuous.service.discord_message_service import DiscordMessageService, ChannelPaginator, Paginator, UserPaginator
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

#print(f'{self.get_random_emoji()}')
        
class Hybrid(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
        self.log_channels: dict[int, list[dict]] = {}
        self.super = {"state": False, "members": set()}
        self.temp_rooms: dict[int, dict[str, TempChannel]] = {}
        self._loaded_aliases = set()

    async def cog_load(self) -> None:
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT guild_id, alias_type, alias_name, channel_id, role_id, room_name FROM command_aliases'
            )
            for row in rows:
                guild_id = row['guild_id']
                alias_type = row['alias_type']
                alias_name = row['alias_name']
                room_name = row['room_name']
                channel_id = row['channel_id']
                role_id = row.get('role_id')
                if room_name != '':
                    self.bot.command_aliases.setdefault(guild_id, self.bot.command_aliases.default_factory()).setdefault('temp_room_aliases', {}).setdefault(alias_type, {})[alias_name] = {
                        'room_name': room_name,
                        'channel_id': int(channel_id) if channel_id else None,
                        'role_id': int(role_id) if role_id else None
                    }
                elif alias_type in ('role', 'unrole'):
                    if role_id:
                        self.bot.command_aliases.setdefault(guild_id, self.bot.command_aliases.default_factory()).setdefault('role_aliases', {}).setdefault(alias_type, {})[alias_name] = {
                            'role_id': int(role_id),
                            'channel_id': int(channel_id) if channel_id else None
                        }
                elif channel_id:
                    self.bot.command_aliases.setdefault(guild_id, self.bot.command_aliases.default_factory()).setdefault('channel_aliases', {}).setdefault(alias_type, {})[alias_name] = int(channel_id)
                if alias_name not in self._loaded_aliases:
                    cmd = None
                    if alias_type == 'mute': cmd = self.create_voice_mute_alias(alias_name)
                    elif alias_type == 'unmute': cmd = self.create_unmute_alias(alias_name)
                    elif alias_type == 'ban': cmd = self.create_ban_alias(alias_name)
                    elif alias_type == 'unban': cmd = self.create_unban_alias(alias_name)
                    elif alias_type == 'cow': cmd = self.create_cow_alias(alias_name)
                    elif alias_type == 'uncow': cmd = self.create_uncow_alias(alias_name)
                    elif alias_type == 'flag': cmd = self.create_flag_alias(alias_name)
                    elif alias_type == 'unflag': cmd = self.create_unflag_alias(alias_name)
                    elif alias_type == 'tmute': cmd = self.create_text_mute_alias(alias_name)
                    elif alias_type == 'untmute': cmd = self.create_untextmute_alias(alias_name)
                    elif alias_type == 'role': cmd = self.create_role_alias(alias_name)
                    elif alias_type == 'unrole': cmd = self.create_unrole_alias(alias_name)
                    if cmd and not self.bot.get_command(alias_name):
                        self.bot.add_command(cmd)
                        self._loaded_aliases.add(alias_name)
        await self.load_log_channels()

    async def load_temp_rooms(self):
        async with self.bot.db_pool.acquire() as conn:
            temp_rows = await conn.fetch(
                'SELECT guild_snowflake, room_name, room_snowflake, owner_snowflake FROM temporary_rooms'
            )
            for row in temp_rows:
                guild_id = row['guild_snowflake']
                room_name = row['room_name']
                room_snowflake = row['room_snowflake']
                guild = self.bot.get_guild(guild_id)
                if not guild:
                    continue
                channel_obj = guild.get_channel(room_snowflake)
                if not channel_obj:
                    continue
                temp_channel = TempChannel(channel_obj, room_name)
                self.temp_rooms.setdefault(guild_id, {})[room_name] = temp_channel
    
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
    
    async def unrestrict(self, guild, member):
        uid = member.id
        async with self.bot.db_pool.acquire() as conn:
            ban_rows = await conn.fetch('SELECT channel_id FROM active_bans WHERE discord_snowflake=$1', uid)
            mute_rows = await conn.fetch('SELECT channel_id FROM active_voice_mutes WHERE discord_snowflake=$1', uid)
            text_rows = await conn.fetch('SELECT channel_id FROM active_text_mutes WHERE discord_snowflake=$1', uid)
        for r in ban_rows:
            try: await guild.unban(discord.Object(id=uid), reason='Toggle OFF')
            except: pass
        for r in mute_rows:
            ch = guild.get_channel(r['channel_id'])
            if ch and member.voice and member.voice.mute: await member.edit(mute=False)
        for r in text_rows:
            ch = guild.get_channel(r['channel_id'])
            text_mute_role = discord.utils.get(guild.roles, name='TextMuted')
            if ch and text_mute_role and text_mute_role in member.roles: await member.remove_roles(text_mute_role)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_bans WHERE discord_snowflake=$1', uid)
            await conn.execute('DELETE FROM active_voice_mutes WHERE discord_snowflake=$1', uid)
            await conn.execute('DELETE FROM active_text_mutes WHERE discord_snowflake=$1', uid)

    def create_ban_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Ban a user from a voice channel.')
        @is_owner_developer_coordinator_moderator_predicator('ban')
        async def ban_alias_text_command(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            duration: Optional[str] = commands.parameter(default='24', description='(+|-)duration(m|h|d) \n 0 - permanent / 24h - default \n `+` to append, `-` to delete, `=` to overwrite reason'),
            *,
            reason: Optional[str] = commands.parameter(default='', description='Optional reason (required for 7 days or more)')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            channel_aliases = self.bot.command_aliases.get(
                ctx.guild.id,
                self.bot.command_aliases.default_factory()
            ).get('channel_aliases', {}).get('ban', {})
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj:
                if member_obj.id in self.super['members']:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB You cannot ban a superhero.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot ban the bot.')
            alias_entry = channel_aliases.get(command_name)
            static_channel_id = None
            room_name = ''
            if alias_entry is not None:
                static_channel_id = int(alias_entry)
            else:
                temp_aliases = self.bot.command_aliases.get(
                    ctx.guild.id,
                    self.bot.command_aliases.default_factory()
                ).get('temp_room_aliases', {}).get('ban', {})
                alias_entry = temp_aliases.get(command_name)
                room_name = alias_entry.get('room_name', '')
                static_channel_id = alias_entry.get('channel_id')
            channel_obj = None
            if static_channel_id is not None:
                channel_obj = await self.resolve_channel(ctx, static_channel_id)
                if not channel_obj:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias')
            if room_name:
                channel_obj = discord.utils.get(ctx.guild.channels, name=room_name)
                if not channel_obj:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from the room name `{room_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                target_name = channel_obj.mention if channel_obj else room_name or ''
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No permission to use `{command_name}` in {target_name}.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                target_name = channel_obj.mention if channel_obj else room_name or ''
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Cannot ban `{highest_role}` in {target_name}.')
            async with self.bot.db_pool.acquire() as conn:
                existing_ban = await conn.fetchrow('''
                    SELECT expires_at, reason
                    FROM active_bans
                    WHERE guild_id = $1
                        AND discord_snowflake = $2
                        AND channel_id = $3
                        AND room_name = $4
                ''', ctx.guild.id, member_obj.id, static_channel_id, room_name)
                stripped = duration.strip() if duration else ''
                is_modification = (
                    (existing_ban is not None and stripped in ('+', '-', '=')) or
                    (duration in ('0', '0h', '0hr', '0hrs', '0hour', '0hours',
                                    '0m', '0min', '0mins', '0minute', '0minutes',
                                    '0d', '0day', '0days'))
                )
                base_time = existing_ban['expires_at'] if existing_ban else None
                if not existing_ban and stripped in ('+', '-', '='):
                    return await self.handler.send_message(ctx, content='\U0001F6AB There is no existing ban to modify.')
                is_coordinator = await is_owner_developer_coordinator_via_alias(ctx, 'ban')
                if stripped in ('+', '-', '='):
                    expires_at = base_time
                    duration_display = self.fmt_duration(base_time)
                else:
                    try:
                        expires_at, duration_display = self.parse_duration(duration, base=base_time)
                    except ValueError:
                        return await self.handler.send_message(ctx, content='\U0001F6AB Invalid duration format.')
                is_relative_duration = stripped.startswith('+') and (len(stripped) > 1 and stripped[1].isdigit())
                is_reason_append = stripped == '+' or (stripped.startswith('+') and not stripped[1].isdigit())
                is_reason_set = stripped.startswith('=')
                is_reason_delete = stripped == '-'
                updated_reason = existing_ban['reason'] if existing_ban and not is_modification else reason
                if existing_ban and (is_reason_append or is_reason_set or is_reason_delete):
                    if is_reason_append:
                        new_text = reason.strip() if reason else ''
                        if not new_text:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason to append.')
                        updated_reason = f'{updated_reason}\n{new_text}' if updated_reason else new_text
                    elif is_reason_set:
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset ban reasons.')
                        updated_reason = reason.strip() if reason else ''
                        if not updated_reason:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason after "=" to set.')
                    elif is_reason_delete:
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can delete ban reasons.')
                        updated_reason = None
                if is_modification and not is_relative_duration and not (is_reason_append or is_reason_set or is_reason_delete):
                    if not is_coordinator:
                        allowed = False
                        if expires_at and existing_ban['expires_at']:
                            caps = await self.get_caps_for_channel(ctx.guild.id, channel_obj.id)
                            active_cap = next((c for c in caps if c[0] == 'ban'), None)
                            cap_expires_at, _ = self.parse_duration(active_cap[1]) if active_cap else (timedelta(days=7) + datetime.now(timezone.utc), None)
                            if expires_at < existing_ban['expires_at'] and expires_at <= cap_expires_at:
                                allowed = True
                        if not allowed:
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
                if expires_at is None or (expires_at - now) > timedelta(days=7):
                    if not is_coordinator:
                        return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can ban permanently or longer than 7 days.')
                    if not reason.strip() and not is_reason_set:
                        return await self.handler.send_message(ctx, content='\U0001F6AB A reason is required for permanent bans or those longer than 7 days.')
                if existing_ban and (is_reason_append or is_reason_set or is_reason_delete):
                    expires_at = existing_ban['expires_at']
                else:
                    if existing_ban:
                        if not is_coordinator and is_relative_duration:
                            pass
                        elif not is_coordinator and expires_at and existing_ban['expires_at']:
                            caps = await self.get_caps_for_channel(ctx.guild.id, channel_obj.id)
                            active_cap = next((c for c in caps if c[0] == 'ban'), None)
                            cap_expires_at, _ = self.parse_duration(active_cap[1]) if active_cap else (timedelta(days=7) + datetime.now(timezone.utc), None)
                            if expires_at < existing_ban['expires_at'] and expires_at <= cap_expires_at:
                                pass
                            else:
                                remaining = existing_ban['expires_at'] - discord.utils.utcnow()
                                hours_left = round(remaining.total_seconds() / 3600, 1)
                                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already banned in {channel_obj.mention} for another {hours_left}h.')
                        elif is_coordinator:
                            pass
                        else:
                            remaining = existing_text_mute['expires_at'] - discord.utils.utcnow()
                            hours_left = round(remaining.total_seconds() / 3600, 1)
                            return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already banned in {channel_obj.mention} for another {hours_left}h.')
                if expires_at and expires_at <= now:
                    return await self.handler.send_message(ctx, content='\U0001F6AB You cannot reduce a ban below the current time.')
            try:
                await channel_obj.set_permissions(member_obj, view_channel=False, reason=f'{reason or "No reason provided"}')
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} was not successfully banned.', allowed_mentions=discord.AllowedMentions.none())
            is_in_channel = False
            if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                is_in_channel = True
                try:
                    await member_obj.move_to(None, reason=f'{reason or "No reason provided"}')
                except discord.Forbidden:
                    await self.handler.send_message(ctx, content=f'\U0001F6AB Could not disconnect {member_obj.mention} from {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                except Exception as e:
                    logger.exception(f'Unexpected error while disconnecting user: {e}')
                    raise
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute('''
                    INSERT INTO active_bans (guild_id, discord_snowflake, channel_id, reason, expires_at, room_name)
                    VALUES ($1,$2,$3,$4,$5,$6)
                    ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name) DO UPDATE
                    SET reason=$4, expires_at=$5
                ''', ctx.guild.id, member_obj.id, channel_obj.id, updated_reason or 'No reason provided', expires_at, room_name or '')
                await conn.execute('''
                    INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                    VALUES ($1, $2, $3, $4, $5, $6)
                ''', 'ban', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, updated_reason or 'No reason provided')
            embed = discord.Embed(
                title=f"{self.get_random_emoji()} {member_obj.display_name} has been banned",
                description=f"**User:** {member_obj.mention}\n**Channel:** {channel_obj.mention}\n**Duration:** {duration_display}\n**Reason:** {updated_reason or 'No reason provided'}",
                color=discord.Color.orange()
            )
            await self.handler.send_message(ctx, embed=embed, allowed_mentions=discord.AllowedMentions.none())
            await self.send_log(ctx, 'ban', member_obj, channel_obj, duration_display, updated_reason or 'No reason provided', ctx.author, expires_at, command_name, is_in_channel, is_modification, highest_role)
        return ban_alias_text_command
            
    def create_cow_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Label a user as going vegan for tracking purposes.')
        @is_owner_developer_coordinator_moderator_predicator('cow')
        async def going_vegan_alias_text_command(
                ctx: commands.Context,
                member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            channel_id = (
                self.bot.command_aliases
                    .get(ctx.guild.id, self.bot.command_aliases.default_factory())
                    .get('channel_aliases', {})
                    .get('cow', {})
                    .get(command_name)
            )
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot cow the bot.')
            channel_obj = await self.resolve_channel(ctx, channel_id)
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to cow this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            select_sql = '''
                SELECT 1
                FROM active_cows
                WHERE discord_snowflake = $1
                AND channel_id = $2
            '''
            insert_cow_sql = '''
                INSERT INTO active_cows (guild_id, discord_snowflake, channel_id, created_at)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_id, discord_snowflake, channel_id) DO NOTHING
            '''
            insert_log_sql = '''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1, $2, $3, $4, $5, $6)
                RETURNING created_at
            '''
            async with self.bot.db_pool.acquire() as conn:
                already_cowed = await conn.fetchval(select_sql, member_obj.id, channel_obj.id)
                if already_cowed:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already going vegan.', allowed_mentions=discord.AllowedMentions.none())
                created_at = await conn.fetchval(insert_log_sql, 'cow', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Cowed a user')
                await conn.execute(insert_cow_sql, ctx.guild.id, member_obj.id, channel_obj.id, created_at)
                await self.handler.send_message(ctx, content=f'\U0001F525 {member_obj.mention} is going vegan!!! \U0001F525', allowed_mentions=discord.AllowedMentions.none())
        return going_vegan_alias_text_command
    
    def create_flag_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Flag a user for the channel mapped to this alias.')
        @is_owner_developer_coordinator_moderator_predicator('flag')
        async def flag_alias_text_command(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            *,
            reason: Optional[str] = commands.parameter(default='', description='An optional reason. For modifications, `+` to append the reason, `=` to overwrite the reason, `-` to delete the reason')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            channel_id = (
                self.bot.command_aliases
                    .get(ctx.guild.id, self.bot.command_aliases.default_factory())
                    .get('channel_aliases', {})
                    .get('flag', {})
                    .get(command_name)
            )
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot flag the bot.')
            if member_obj:
                if member_obj.id in self.super['members']:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB You cannot flag a superhero.')
            channel_obj = await self.resolve_channel(ctx, channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}.')
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
            async with self.bot.db_pool.acquire() as conn:
                existing_flag = await conn.fetchrow(select_sql, ctx.guild.id, member_obj.id, channel_obj.id)
                is_modification = existing_flag is not None
                stripped = reason.strip() if reason else ''
                is_reason_append = stripped.startswith('+') and not stripped[1].isdigit()
                is_reason_set = stripped.startswith('=')
                is_reason_delete = stripped.startswith('-')
                cleaned_reason = stripped[1:].strip() if (is_reason_append or is_reason_set or is_reason_delete) else stripped
                updated_reason = cleaned_reason
                if not is_modification and (is_reason_append or is_reason_set or is_reason_delete):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Cannot modify a flag that does not exist. Use a plain reason to create one.')
                if is_modification and (is_reason_append or is_reason_set or is_reason_delete):
                    if is_reason_append:
                        if not cleaned_reason:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason to append.')
                        updated_reason = f'{existing_flag['reason']}\n{cleaned_reason}' if existing_flag['reason'] else cleaned_reason
                    elif is_reason_set:
                        is_coordinator = await is_owner_developer_coordinator_via_alias(ctx, 'flag')
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset flag reasons.')
                        if not cleaned_reason:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason after "=" to set.')
                        updated_reason = cleaned_reason
                    elif is_reason_delete:
                        is_coordinator = await is_owner_developer_coordinator_via_alias(ctx, 'flag')
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can delete flag reasons.')
                        updated_reason = None
                    await conn.execute(update_sql, ctx.guild.id, member_obj.id, channel_obj.id, updated_reason)
                elif not is_modification:
                    updated_reason = cleaned_reason if cleaned_reason and not is_modification else 'No reason provided'
                    await conn.execute(insert_sql, ctx.guild.id, member_obj.id, channel_obj.id, updated_reason)
                else:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already flagged in {channel_obj.mention} for {existing_flag["reason"]}.', allowed_mentions=discord.AllowedMentions.none())
                embed = discord.Embed(color=discord.Color.orange())
                embed.set_author(name=f'{member_obj.display_name} is flagged', icon_url=member_obj.display_avatar.url)
                embed.add_field(name='User', value=member_obj.mention, inline=True)
                embed.add_field(name='Channel', value=channel_obj.mention, inline=False)
                full_reason = None
                if existing_flag:
                    if is_reason_append: full_reason = f'{existing_flag["reason"]} {cleaned_reason}'.strip()
                    elif is_reason_delete: full_reason = ''
                    else: full_reason = cleaned_reason
                else: full_reason = cleaned_reason
                embed.add_field(name='Reason', value=full_reason or 'No reason provided', inline=False)
                await self.handler.send_message(ctx, embed=embed, allowed_mentions=discord.AllowedMentions.none())
        return flag_alias_text_command
   
    def create_role_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help=f'Gives a specific role to a user.')
        @is_owner_developer_coordinator_moderator_predicator('role')
        async def role_alias_text_command(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            alias_data = (
                self.bot.command_aliases
                    .get(ctx.guild.id, self.bot.command_aliases.default_factory())
                    .get('role_aliases', {})
                    .get('role', {})
                    .get(command_name)
            )
            if not alias_data:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No role alias configured for `{command_name}`.')
            static_role_id = int(alias_data.get('role_id'))
            target_channel_id = int(alias_data.get('channel_id')) if alias_data.get('channel_id') else None
            is_owner_or_dev, is_coord_or_mod = await check_owner_dev_coord_mod(ctx, target_channel_id)
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot give the bot a role.')
            channel_obj = await self.resolve_channel(ctx, target_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            if not is_owner_or_dev and not is_coord_or_mod:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
            role_obj = ctx.guild.get_role(static_role_id)
            if not role_obj:
                return await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Could not resolve role with ID `{static_role_id}`.')
            if role_obj in member_obj.roles:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB{member_obj.mention} already has {role_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            try:
                await member_obj.add_roles(role_obj, reason=f'Added role')
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} was not successfully roled.', allowed_mentions=discord.AllowedMentions.none())
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} was given {role_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
        return role_alias_text_command
            
    def create_text_mute_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name,help='Text mutes a user in a specific text channel or temporary room.')
        @is_owner_developer_coordinator_moderator_predicator('tmute')
        async def text_mute_alias_text_command(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            duration: Optional[str] = commands.parameter(default='8', description='(+|-)duration(m|h|d) \n 0 - permanent / 8h - default \n `+` to append, `-` to delete, `=` to overwrite reason'),
            *,
            reason: Optional[str] = commands.parameter(default='', description='Optional reason (required for 7 days or more)')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            channel_aliases = self.bot.command_aliases.get(
                ctx.guild.id,
                self.bot.command_aliases.default_factory()
            ).get('channel_aliases', {}).get('tmute', {})
            alias_entry = channel_aliases.get(command_name)
            static_channel_id = None
            room_name = ''
            if alias_entry is not None:
                static_channel_id = int(alias_entry)
            else:
                temp_aliases = self.bot.command_aliases.get(
                    ctx.guild.id,
                    self.bot.command_aliases.default_factory()
                ).get('temp_room_aliases', {}).get(alias_type, {})
                alias_entry = temp_aliases.get(command_name)
                if alias_entry is None:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No alias configured for `{command_name}`.')
                room_name = alias_entry.get('room_name', '')
                static_channel_id = alias_entry.get('channel_id')
            channel_obj = None
            if static_channel_id is not None:
                channel_obj = await self.resolve_channel(ctx, static_channel_id)
                if not channel_obj:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            if room_name:
                channel_obj = discord.utils.get(ctx.guild.channels, name=room_name)
                if not channel_obj:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from the room name `{room_name}`.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or member_obj.id in self.super['members']:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid target member: {member}.')
            is_owner_or_dev, is_mod_or_coord = (await check_owner_dev_coord_mod(ctx, channel_obj) if channel_obj else (False, True))
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No permission to use `{command_name}` in {channel_obj.mention}.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj) if channel_obj else (None, True)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to mute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                existing_text_mute = await conn.fetchrow('''
                    SELECT expires_at, reason
                    FROM active_text_mutes
                    WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND room_name=$4
                ''', ctx.guild.id, member_obj.id, static_channel_id, room_name)
                stripped = duration.strip() if duration else ''
                base_time = existing_text_mute['expires_at'] if existing_text_mute else None
                try:
                    expires_at, duration_display = (
                        (base_time, self.fmt_duration(base_time))
                        if stripped in ('+', '-', '=')
                        else self.parse_duration(duration, base=base_time)
                    )
                except ValueError:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Invalid duration format.')
                updated_reason = existing_text_mute['reason'] if existing_text_mute else reason
                if existing_text_mute and stripped in ('+', '-', '='):
                    if stripped == '+':
                        updated_reason = f'{existing_text_mute['reason']}\n{reason}' if existing_text_mute['reason'] else reason
                    elif stripped == '=':
                        if not await is_owner_developer_coordinator_via_alias(ctx,'tmute'):
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset reasons.')
                        updated_reason = reason
                    elif stripped == '-':
                        if not await is_owner_developer_coordinator_via_alias(ctx,'tmute'):
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can delete reasons.')
                        updated_reason = None
                if expires_at and expires_at <= datetime.now(timezone.utc):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Cannot set a mute in the past.')
                if channel_obj:
                    try:
                        await channel_obj.set_permissions(member_obj, send_messages=False, add_reactions=False)
                    except discord.Forbidden:
                        return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} was not successfully text-muted.', allowed_mentions=discord.AllowedMentions.none())
                await conn.execute('''
                    INSERT INTO active_text_mutes (guild_id, discord_snowflake, channel_id, reason, expires_at, room_name)
                    VALUES ($1,$2,$3,$4,$5,$6)
                    ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name) DO UPDATE
                    SET reason=$4, expires_at=$5
                ''', ctx.guild.id, member_obj.id, static_channel_id, updated_reason or 'No reason provided', expires_at, room_name)
                await conn.execute(
                    '''
                    INSERT INTO moderation_logs
                    (action_type,target_discord_snowflake,executor_discord_snowflake,guild_id,channel_id,reason)
                    VALUES ($1,$2,$3,$4,$5,$6)
                    ''', 'textmute', member_obj.id, ctx.author.id, ctx.guild.id, static_channel_id, f'Textmuted ({updated_reason or "No reason"})'
                )
            embed = discord.Embed(
                title=f"{self.get_random_emoji()} {member_obj.display_name} is text-muted",
                description=f"**User:** {member_obj.mention}\n**Channel/Room:** {channel_obj.mention}\n**Duration:** {duration_display}\n**Reason:** {updated_reason or 'No reason provided'}",
                color=discord.Color.orange()
            )
            await self.handler.send_message(ctx, embed=embed, allowed_mentions=discord.AllowedMentions.none())
            return await self.send_log(ctx, 'tmute', member_obj, channel_obj, duration_display, updated_reason,ctx.author, expires_at, command_name, True, stripped in ('+', '-', '='), highest_role)
        return text_mute_alias_text_command
        
    def create_voice_mute_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Mutes a member in a specific VC.')
        @is_owner_developer_coordinator_moderator_predicator('mute')
        async def voice_mute_alias_text_command(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            duration: Optional[str] = commands.parameter(default='8', description='(+|-)duration(m|h|d) \n 0 - permanent / 8h - default \n `+` to append, `-` to delete, `=` to overwrite reason'),
            *,
            reason: Optional[str] = commands.parameter(default='', description='Optional reason (required for 7 days or more)')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            channel_aliases = self.bot.command_aliases.get(
                ctx.guild.id,
                self.bot.command_aliases.default_factory()
            ).get('channel_aliases', {}).get('mute', {})
            alias_entry = channel_aliases.get(command_name)
            static_channel_id = None
            room_name = ''
            if alias_entry is not None:
                static_channel_id = int(alias_entry)
            else:
                temp_aliases = self.bot.command_aliases.get(
                    ctx.guild.id,
                    self.bot.command_aliases.default_factory()
                ).get('temp_room_aliases', {}).get('mute', {})
                alias_entry = temp_aliases.get(command_name)
                if alias_entry is None:
                    return await self.handler.send_message(
                        ctx, content=f'\U0001F6AB No alias configured for `{command_name}`.'
                    )
                room_name = alias_entry.get('room_name', '')
                static_channel_id = alias_entry.get('channel_id')
            channel_obj = None
            if static_channel_id is not None:
                channel_obj = await self.resolve_channel(ctx, static_channel_id)
                if not channel_obj:
                    return await self.handler.send_message(
                        ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.'
                    )
            if room_name:
                channel_obj = discord.utils.get(ctx.guild.channels, name=room_name)
                if not channel_obj:
                    return await self.handler.send_message(
                        ctx, content=f'\U0001F6AB Could not resolve a valid channel from the room name `{room_name}`.'
                    )
            member_obj = await self.resolve_member(ctx, member)
            if member_obj:
                if member_obj.id in self.super['members']:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB You cannot mute a superhero.')
            if not member_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot voice mute the bot.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No permission to use `{command_name}` in {channel_obj.mention}.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to mute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                existing_mute = await conn.fetchrow('''
                    SELECT expires_at, reason
                    FROM active_voice_mutes
                    WHERE guild_id = $1
                      AND discord_snowflake = $2
                      AND channel_id = $3
                      AND room_name = $4
                ''', ctx.guild.id, member_obj.id, static_channel_id, room_name)
                stripped = duration.strip() if duration else ''
                is_modification = (
                    (existing_mute is not None and stripped in ('+', '-', '=')) or
                    (duration in ('0', '0h', '0hr', '0hrs', '0hour', '0hours',
                                    '0m', '0min', '0mins', '0minute', '0minutes',
                                    '0d', '0day', '0days'))
                )
                base_time = existing_mute['expires_at'] if existing_mute else None
                if not existing_mute and stripped in ('+', '-', '='):
                    return await self.handler.send_message(ctx,content='\U0001F6AB There is no existing voice mute to modify.')
                is_coordinator = await is_owner_developer_coordinator_via_alias(ctx, 'mute')
                if stripped in ('+', '-', '='):
                    expires_at = base_time
                    duration_display = self.fmt_duration(base_time)
                else:
                    try:
                        expires_at, duration_display = self.parse_duration(duration, base=base_time)
                    except ValueError:
                        return await self.handler.send_message(ctx, content='\U0001F6AB Invalid duration format.')
                is_relative_duration = stripped.startswith('+') and (len(stripped) > 1 and stripped[1].isdigit())
                is_reason_append = stripped == '+' or (stripped.startswith('+') and not stripped[1].isdigit())
                is_reason_set = stripped == '='
                is_reason_delete = stripped == '-'
                updated_reason = existing_mute['reason'] if existing_mute and not is_modification else reason
                if existing_mute and (is_reason_append or is_reason_set or is_reason_delete):
                    if is_reason_append:
                        new_text = reason.strip() if reason else ''
                        if not new_text:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason to append.')
                        updated_reason = f'{existing_mute["reason"]}\n{new_text}' if updated_reason else new_text
                    elif is_reason_set:
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can reset mute reasons.')
                        updated_reason = reason.strip() if reason else ''
                        if not updated_reason:
                            return await self.handler.send_message(ctx, content='\U0001F6AB You must provide a reason after "=" to set.')
                    elif is_reason_delete:
                        if not is_coordinator:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can delete voice mute reasons.')
                        updated_reason = None
                if is_modification and not is_relative_duration and not (is_reason_append or is_reason_set or is_reason_delete):
                    if not is_coordinator:
                        allowed = False
                        if expires_at and existing_mute['expires_at']:
                            caps = await self.get_caps_for_channel(ctx.guild.id, channel_obj.id)
                            active_cap = next((c for c in caps if c[0] == 'mute'), None)
                            cap_expires_at, _ = self.parse_duration(active_cap[1]) if active_cap else (timedelta(days=7) + datetime.now(timezone.utc), None)
                            if expires_at < existing_mute['expires_at'] and expires_at <= cap_expires_at:
                                allowed = True
                        if not allowed:
                            return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can overwrite an existing mute with an absolute duration.')
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
                if expires_at is None or (expires_at - now) > timedelta(days=7):
                    if not is_coordinator:
                        return await self.handler.send_message(ctx, content='\U0001F6AB Only coordinators can voice mute permanently or longer than 7 days.')
                    if not reason.strip() and not is_reason_set:
                        return await self.handler.send_message(ctx, content='\U0001F6AB A reason is required for permanent voice mutes or those longer than 7 days.')
            if existing_mute and (is_reason_append or is_reason_set or is_reason_delete):
                expires_at = existing_mute['expires_at']
            else:
                if existing_mute:
                    if not is_coordinator and is_relative_duration:
                        pass
                    elif not is_coordinator and expires_at and existing_mute['expires_at']:
                        caps = await self.get_caps_for_channel(ctx.guild.id, channel_obj.id)
                        active_cap = next((c for c in caps if c[0] == 'mute'), None)
                        cap_expires_at, _ = self.parse_duration(active_cap[1]) if active_cap else (timedelta(days=7) + datetime.now(timezone.utc), None)
                        if expires_at < existing_mute['expires_at'] and expires_at <= cap_expires_at:
                            pass
                        else:
                            remaining = existing_mute['expires_at'] - discord.utils.utcnow()
                            hours_left = round(remaining.total_seconds() / 3600, 1)
                            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} is already voice muted in {channel_obj.mention} for another {hours_left}h.')
                    elif is_coordinator:
                        pass
                    else:
                        remaining = existing_mute['expires_at'] - discord.utils.utcnow()
                        hours_left = round(remaining.total_seconds() / 3600, 1)
                        return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is already voice muted in {channel_obj.mention} for another {hours_left}h.')
                if expires_at and expires_at <= now:
                    return await self.handler.send_message(ctx, content='\U0001F6AB You cannot reduce a voice mute below the current time.')
            try:
                async with self.bot.db_pool.acquire() as conn:
                    await conn.execute('''
                        INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, reason, target, room_name)
                        VALUES ($1, $2, $3, $4, $5, $6, $7)
                        ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
                        DO UPDATE SET
                            expires_at = EXCLUDED.expires_at,
                            reason = EXCLUDED.reason
                    ''', ctx.guild.id, member_obj.id, static_channel_id, expires_at, updated_reason or 'No reason provided', 'user', room_name)
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'voice_mute', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Voice muted a member')
            except Exception as e:
                logger.warning(f'DB insert failed: {e}')
                raise
            is_in_channel = False
            if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                is_in_channel = True
                try:
                    await member_obj.edit(mute=True)
                except discord.Forbidden:
                    return await self.handler.send_message(ctx, content=f"\U0001F6AB {member_obj.mention} was not successfully voice muted.", allowed_mentions=discord.AllowedMentions.none())
            embed = discord.Embed(
                    title=f"{self.get_random_emoji()} {member_obj.display_name} is voice muted",
                    description=f"**User:** {member_obj.mention}\n**Channel:** {channel_obj.mention}\n**Duration:** {duration_display}\n**Reason:** {updated_reason or 'No reason provided'}",
                    color=discord.Color.orange()
                )
            await self.handler.send_message(ctx, embed=embed, allowed_mentions=discord.AllowedMentions.none())
            await self.send_log(ctx, 'vmute', member_obj, channel_obj, duration_display, updated_reason or 'No reason provided', ctx.author, expires_at, command_name, is_in_channel, is_modification, highest_role)
        return voice_mute_alias_text_command

    def create_unban_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Unban a user from a voice channel.')
        @is_owner_developer_coordinator_moderator_predicator('unban')
        async def unban_alias_text_command(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
            room_name: Optional[str] = ''
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot unban the bot.')
            channel_aliases = self.bot.command_aliases.get(
                ctx.guild.id,
                self.bot.command_aliases.default_factory()
            ).get('channel_aliases', {}).get('unban', {})
            alias_entry = channel_aliases.get(command_name)
            static_channel_id = None
            room_name = ''
            if alias_entry is not None:
                static_channel_id = int(alias_entry)
            else:
                temp_aliases = self.bot.command_aliases.get(
                    ctx.guild.id,
                    self.bot.command_aliases.default_factory()
                ).get('temp_room_aliases', {}).get('unban', {})
                alias_entry = temp_aliases.get(command_name)
                if alias_entry is None:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No alias configured for `{command_name}`.')
                room_name = alias_entry.get('room_name', '')
                static_channel_id = alias_entry.get('channel_id')
            channel_obj = None
            if static_channel_id is not None:
                channel_obj = await self.resolve_channel(ctx, static_channel_id)
                if not channel_obj:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            if room_name:
                channel_obj = discord.utils.get(ctx.guild.channels, name=room_name)
                if not channel_obj:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from the room name `{room_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj) if channel_obj else (False, False)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to unban in {channel_obj.mention}.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to unban this `{highest_role}`.')
            async with self.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow('''
                    SELECT expires_at
                    FROM active_bans
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3 AND room_name = $4
                ''', ctx.guild.id, member_obj.id, static_channel_id, room_name)
                if row and row['expires_at'] is None and not await is_owner_developer_coordinator_via_alias(ctx, 'unban'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Coordinator-only for undoing permanent bans.')
                try:
                    if channel_obj:
                        await channel_obj.set_permissions(member_obj, overwrite=None)
                except discord.Forbidden:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} was not successfully unbanned.', allowed_mentions=discord.AllowedMentions.none())
                await conn.execute('''
                    DELETE FROM active_bans
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3 AND room_name = $4
                ''', ctx.guild.id, member_obj.id, static_channel_id, room_name)
                await conn.execute('''
                    INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                    VALUES ($1,$2,$3,$4,$5,$6)
                ''', 'unban', member_obj.id, ctx.author.id, ctx.guild.id, static_channel_id if static_channel_id else -1, f'Unbanned a user from room `{room_name}`' if room_name != '' else 'Unbanned a user')
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been unbanned from {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
        return unban_alias_text_command

    def create_uncow_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Unlabel a user for tracking purposes.')
        @is_owner_developer_coordinator_moderator_predicator('uncow')
        async def no_longer_going_vegan_alias_text_command(
                ctx: commands.Context,
                member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            static_channel_id = int(
                self.bot.command_aliases
                    .get(ctx.guild.id, self.bot.command_aliases.default_factory())
                    .get('channel_aliases', {})
                    .get('uncow', {})
                    .get(command_name)
            )
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No uncow alias configured for `{command_name}`.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot uncow the bot.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel for `{command_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
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
                    await self.handler.send_message(ctx, content=f'👎 {member_obj.mention} is no longer going vegan. 👎', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                logger.warning(f'Database error occurred: {e}')
                raise
        return no_longer_going_vegan_alias_text_command
        
    def create_unflag_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Unflag a user in the database for the voice channel mapped to this alias.')
        @is_owner_developer_coordinator_moderator_predicator('unflag')
        async def unflag_alias_text_command(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            static_channel_id = int(
                self.bot.command_aliases
                    .get(ctx.guild.id, self.bot.command_aliases.default_factory())
                    .get('channel_aliases', {})
                    .get('unflag', {})
                    .get(command_name)
            )
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No unflag alias configured for `{command_name}`.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj or not member:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot unflag the bot.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel for `{command_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
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
                        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} is not flagged for {channel_obj.mention}.')
                    await conn.execute(update_sql, ctx.guild.id, member_obj.id, channel_obj.id)
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'unflag', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Unflagged a user')
                    await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Unflagged {member_obj.mention} for channel {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                logger.warning(f'Database error occurred: {e}')
                raise
        return unflag_alias_text_command

    def create_unmute_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Unmutes a member in a specific VC.')
        @is_owner_developer_coordinator_moderator_predicator('unmute')
        async def unmute_alias_text_command(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot unmute the bot.')
            channel_aliases = self.bot.command_aliases.get(
                ctx.guild.id,
                self.bot.command_aliases.default_factory()
            ).get('channel_aliases', {}).get('unmute', {})
            alias_entry = channel_aliases.get(command_name)
            static_channel_id = None
            room_name = ''
            if alias_entry is not None:
                static_channel_id = int(alias_entry)
            else:
                temp_aliases = self.bot.command_aliases.get(
                    ctx.guild.id,
                    self.bot.command_aliases.default_factory()
                ).get('temp_room_aliases', {}).get('unmute', {})
                alias_entry = temp_aliases.get(command_name)
                if alias_entry is None:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No alias configured for `{command_name}`.')
                room_name = alias_entry.get('room_name', '')
                static_channel_id = alias_entry.get('channel_id')
            channel_obj = None
            if static_channel_id is not None:
                channel_obj = await self.resolve_channel(ctx, static_channel_id)
                if not channel_obj:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            if room_name:
                channel_obj = discord.utils.get(ctx.guild.channels, name=room_name)
                if not channel_obj:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from the room name `{room_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to unmute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow('''
                    SELECT expires_at
                    FROM active_voice_mutes
                    WHERE guild_id = $1
                      AND discord_snowflake = $2
                      AND channel_id = $3
                      AND room_name = $4
                      AND target = 'user'
                ''', ctx.guild.id, member_obj.id, channel_obj.id, room_name)
                if not row:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not muted in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                if row['expires_at'] is None and not await is_owner_developer_coordinator_via_alias(ctx, 'unmute'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Coordinator-only for undoing permanent voice mutes.')
                await conn.execute('''
                    DELETE FROM active_voice_mutes
                    WHERE guild_id = $1
                      AND discord_snowflake = $2
                      AND channel_id = $3
                      AND room_name = $4
                      AND target = $5
                ''', ctx.guild.id, member_obj.id, channel_obj.id, room_name, 'user')
                if member_obj.voice and member_obj.voice.channel:
                    try:
                        await member_obj.edit(mute=False)
                    except discord.Forbidden:
                        return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} was not successfully unmuted.', allowed_mentions=discord.AllowedMentions.none())
                await conn.execute('INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1,  $2, $3, $4, $5, $6)', 'unmute', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Unmuted a member')
            if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been unmuted in {channel_obj.mention}.',  allowed_mentions=discord.AllowedMentions.none())
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} is no longer marked as muted in {channel_obj.mention}.',  allowed_mentions=discord.AllowedMentions.none())
        return unmute_alias_text_command

    def create_unrole_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Removes a specific role from a user.')
        @is_owner_developer_coordinator_predicator('unrole')
        async def unrole_alias_text_command(
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            alias_data = (
                self.bot.command_aliases
                    .get(ctx.guild.id, self.bot.command_aliases.default_factory())
                    .get('role_aliases', {})
                    .get('unrole', {})
                    .get(command_name)
            )
            if not alias_data:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No unrole alias configured for `{command_name}`.')
            static_channel_id = int(alias_data.get('channel_id'))
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot unrole the bot.')
            channel_obj = await self.resolve_channel(ctx, static_channel_id)
            if not channel_obj:
                return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel for `{command_name}`.')
            is_owner_or_dev, is_coord = await check_owner_dev_coord(ctx, channel_obj)
            if not is_coord and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to derole this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            static_role_id = int(alias_data.get('role_id'))
            role_obj = ctx.guild.get_role(static_role_id)
            if not role_obj:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve role with ID `{static_role_id}`.')
            if role_obj not in member_obj.roles:
                return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} does not have {role_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            try:
                await member_obj.remove_roles(role_obj)
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} was not successfully unroled.', allowed_mentions=discord.AllowedMentions.none())
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} had {role_obj.mention} removed.', allowed_mentions=discord.AllowedMentions.none())
        return unrole_alias_text_command
        
    def create_untextmute_alias(self, command_name: Optional[str]) -> Command:
        @commands.command(name=command_name, help='Removes a text mute from a user in a specific text channel.')
        @is_owner_developer_coordinator_moderator_predicator('untmute')
        async def untext_mute_alias_text_command(
            ctx,
            member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
        ) -> None:
            if not ctx.guild:
                return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
            member_obj = await self.resolve_member(ctx, member)
            if not member_obj:
                 return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
            if member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot undo a textmute on the bot.')
            channel_aliases = self.bot.command_aliases.get(
                ctx.guild.id,
                self.bot.command_aliases.default_factory()
            ).get('channel_aliases', {}).get('untmute', {})
            alias_entry = channel_aliases.get(command_name)
            static_channel_id = None
            room_name = ''
            if alias_entry is not None:
                static_channel_id = int(alias_entry)
            else:
                temp_aliases = self.bot.command_aliases.get(
                    ctx.guild.id,
                    self.bot.command_aliases.default_factory()
                ).get('temp_room_aliases', {}).get('untmute', {})
                alias_entry = temp_aliases.get(command_name)
                if alias_entry is None:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No alias configured for `{command_name}`.')
                room_name = alias_entry.get('room_name', '')
                static_channel_id = alias_entry.get('channel_id')
            channel_obj = None
            if static_channel_id is not None:
                channel_obj = await self.resolve_channel(ctx, static_channel_id)
                if not channel_obj:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve a valid channel from the alias.')
            if room_name:
                channel_obj = discord.utils.get(ctx.guild.channels, name=room_name)
                if not channel_obj:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from the room name `{room_name}`.')
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_mod_or_coord:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`{command_name}`) in {channel_obj.mention}.')
            highest_role, success = await check_block(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to un textmute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow('''
                    SELECT expires_at
                    FROM active_text_mutes
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3 AND room_name = $4
                ''', ctx.guild.id, member_obj.id, static_channel_id, room_name)
                if row and row['expires_at'] is None and not await is_owner_developer_coordinator_via_alias(ctx, 'untmute'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Coordinator-only for undoing permanent text mutes.')
                try:
                    await channel_obj.set_permissions(member_obj, send_messages=None)
                except discord.Forbidden:
                    return await self.handler.send_message(ctx, content=f"\U0001F6AB {member_obj.mention} was not successfully untext-muted.", allowed_mentions=discord.AllowedMentions.none())
                await conn.execute('''
                    DELETE FROM active_text_mutes
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3 AND room_name = $4
                ''', ctx.guild.id, member_obj.id, static_channel_id, room_name)
                await conn.execute('INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1, $2, $3, $4, $5, $6)', 'untmute', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Untextmuted a user')
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention}\'s text muted in {channel_obj.mention} has been removed.', allowed_mentions=discord.AllowedMentions.none())
        return untext_mute_alias_text_command
            
    @app_commands.command(name='admin', description='Grants server mute privileges to a member for the entire guild.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_app_predicator()
    async def create_administrator_app_command(
        self,
        interaction: discord.Interaction,
        member: str
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in a server.')
        member_obj = await self.resolve_member_app(interaction, member)
        if member_obj.bot:
            return await self.send(interaction, content='\U0001F6AB You cannot make the bot an administrator.')
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
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
            ''', member_obj.id, interaction.guild.id)
        return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has been granted server mute permissions.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='admin', help='Grants server mute privileges to a member for the entire guild.')
    @is_owner_predicator()
    async def create_administrator_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(description='Tag a member or include their snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot an administrator.')
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
        await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been granted server mute permissions.', allowed_mentions=discord.AllowedMentions.none())
            
    @app_commands.command(name='admins', description='Lists all members with server mute privileges in this guild.')
    @is_owner_app_predicator()
    async def list_server_muters_app_command(
        self,
        interaction: discord.Interaction
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in a server.')
        async with self.bot.db_pool.acquire() as conn:
            records = await conn.fetch('''SELECT discord_snowflake FROM users WHERE $1=ANY(server_muter_guild_ids) ORDER BY discord_snowflake''', interaction.guild.id)
            if not records:
                return await self.send(interaction, content=f'\U0001F6AB No admins found in {guild.name}.', allowed_mentions=discord.AllowedMentions.none())
            description_lines = []
            for record in records:
                uid = record['discord_snowflake']
                member_obj = await self.resolve_member_app(interaction, uid)
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
            
    @commands.command(name='admins', help='Lists all members with server mute privileges in this guild.')
    @is_owner_predicator()
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
            
    @app_commands.command(name='alias', description='Set an alias for a Vyrtuous action.')
    @is_owner_developer_coordinator_app_predicator(None)
    @app_commands.describe(
        alias_type='One of: cow, uncow, mute, unmute, ban, unban, flag, unflag, tmute, untmute, role, unrole',
        alias_name='Alias/Pseudonym',
        channel='Tag a channel or include its snowflake ID',
        role='Role ID (only for role/unrole)'
    )
    async def create_alias_app_command(
        self,
        interaction: discord.Interaction,
        alias_type: Optional[str] = None,
        alias_name: Optional[str] = None,
        channel: Optional[str] = None,
        role: Optional[str] = None
    ):
        valid_types = {'cow','uncow','mute','unmute','ban','unban','flag','unflag','tmute','untmute','role','unrole'}
        if not alias_type or alias_type.lower() not in valid_types:
            return await self.send(interaction, content=f'\U0001F6AB Invalid alias type. Must be one of: `{"`, `".join(valid_types)}`')
        alias_type = alias_type.lower()
        if not alias_name or not alias_name.strip():
            return await self.send(interaction, content='\U0001F6AB Alias name cannot be empty.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        room_name = ''
        temp_room_obj = None
        for guild_temp_rooms in self.temp_rooms.values():
            for temp_channel in guild_temp_rooms.values():
                if temp_channel.channel.id == channel_obj.id:
                    temp_room_obj = temp_channel
                    room_name = temp_channel.room_name
                    break
            if temp_room_obj:
                break
        if not channel_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        is_owner_or_dev, _ = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        if not is_owner_or_dev and is_temp_room:
            async with self.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(
                    'SELECT coordinator_room_names, coordinator_channel_ids FROM users WHERE discord_snowflake=$1',
                    interaction.user.id
                )
            allowed_rooms = row.get('coordinator_room_names') or []
            allowed_channels = row.get('coordinator_channel_ids') or []
            if channel_obj.name not in allowed_rooms and channel_obj.id not in allowed_channels:
                return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use this command in `{channel_obj.name}`.')
        async with self.bot.db_pool.acquire() as conn:
            existing_alias = await conn.fetchrow('''
                SELECT guild_id, channel_id, role_id, room_name
                FROM command_aliases
                WHERE guild_id=$1 AND alias_type=$2 AND alias_name=$3 AND room_name=$4
            ''', interaction.guild.id, alias_type, alias_name, room_name)
            if existing_alias:
                existing_channel = interaction.guild.get_channel(existing_alias['channel_id']) if existing_alias['channel_id'] else None
                existing_role = interaction.guild.get_role(existing_alias['role_id']) if existing_alias['role_id'] else None
                existing_room = existing_alias.get('room_name')
                mention = existing_channel.mention if existing_channel else (existing_role.mention if existing_role else (existing_room or 'unknown'))
                return await self.send(interaction, content=f'\U0001F6AB Alias `{alias_name}` ({alias_type}) already exists and is set to {mention}.')
        if self.bot.get_command(alias_name):
            return await self.send(interaction, content=f'\U0001F6AB A command named `{alias_name}` already exists.')
        if alias_type in ('role','unrole'):
            if not role:
                return await self.send(interaction, content='\U0001F6AB Role ID is required for role/unrole aliases.')
            try:
                role_id = int(role.replace('<@&','').replace('>',''))
            except ValueError:
                return await self.send(interaction, content=f'\U0001F6AB Invalid role ID: {role}')
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute('''
                    INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id, role_id, room_name)
                    VALUES ($1,$2,$3,$4,$5,$6)
                ''', interaction.guild.id, alias_type, alias_name, channel_obj.id, role_id, room_name)
        else:
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute('''
                    INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id, role_id, room_name)
                    VALUES ($1,$2,$3,$4,$5,$6)
                ''', interaction.guild.id, alias_type, alias_name, channel_obj.id, None, room_name)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'create_alias', None, interaction.user.id, interaction.guild.id, channel_obj.id, f'Created an alias: {alias_name}')
        if alias_type in ('role','unrole') and is_owner_or_dev:
            self.bot.command_aliases.setdefault(interaction.guild.id, self.bot.command_aliases.default_factory()).setdefault('role_aliases', {}).setdefault(alias_type, {})[alias_name] = {'channel_id': channel_obj.id, 'role_id': role_id}
        elif room_name != '':
            self.bot.command_aliases.setdefault(interaction.guild.id, self.bot.command_aliases.default_factory()).setdefault('temp_room_aliases', {}).setdefault(alias_type, {})[alias_name] = {'room_name': room_name, 'channel_id': channel_obj.id}
        else:
            self.bot.command_aliases.setdefault(interaction.guild.id, self.bot.command_aliases.default_factory()).setdefault('channel_aliases', {}).setdefault(alias_type, {})[alias_name] = channel_obj.id
        cmd = None
        if alias_type == 'ban': cmd = self.create_ban_alias(alias_name)
        elif alias_type == 'flag': cmd = self.create_flag_alias(alias_name)
        elif alias_type == 'unflag': cmd = self.create_unflag_alias(alias_name)
        elif alias_type == 'mute': cmd = self.create_voice_mute_alias(alias_name)
        elif alias_type == 'tmute': cmd = self.create_text_mute_alias(alias_name)
        elif alias_type == 'unban': cmd = self.create_unban_alias(alias_name)
        elif alias_type == 'cow':
            if channel_obj.id != 1222056499959042108: return await self.send(interaction, content=f'\U0001F6AB {channel_obj.mention} has no vegan association.')
            cmd = self.create_cow_alias(alias_name)
        elif alias_type == 'uncow':
            if channel_obj.id != 1222056499959042108: return await self.send(interaction, content=f'\U0001F6AB {channel_obj.mention} has no vegan association.')
            cmd = self.create_uncow_alias(alias_name)
        elif alias_type == 'unmute': cmd = self.create_unmute_alias(alias_name)
        elif alias_type == 'untmute': cmd = self.create_untextmute_alias(alias_name)
        elif alias_type == 'role': cmd = self.create_role_alias(alias_name)
        elif alias_type == 'unrole': cmd = self.create_unrole_alias(alias_name)
        if cmd: self.bot.add_command(cmd)
        mention = interaction.guild.get_role(role_id).mention if alias_type in ('role','unrole') else channel_obj.mention
        await self.send(interaction, content=f'{self.get_random_emoji()} Alias `{alias_name}` ({alias_type}) set to {mention}.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(
        name='alias',
        help='Set an alias for a cow, uncow, mute, unmute, ban, unban, flag, unflag, tmute, untmute, role, or unrole action.'
    )
    @is_owner_developer_coordinator_predicator(None)
    async def create_alias_text_command(
        self,
        ctx: commands.Context,
        alias_type: Optional[str] = commands.parameter(default=None, description='One of: `cow`, `uncow`, `mute`, `unmute`, `ban`, `unban`, `flag`, `unflag`, `tmute`, `untmute`, `role`, `unrole`'),
        alias_name: Optional[str] = commands.parameter(default=None, description='Alias/Pseudonym'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        *,
        role: Optional[str] = commands.parameter(default=None, description='Role ID (only for role/unrole)')
    ) -> None:
        valid_types = {'cow','uncow','mute','unmute','ban','unban','flag','unflag','cow','uncow','tmute','untmute','role','unrole'}
        if not alias_type or alias_type.lower() not in valid_types:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid alias type. Must be one of: `{"`, `".join(valid_types)}`')
        alias_type = alias_type.lower()
        if not alias_name or not alias_name.strip():
            return await self.handler.send_message(ctx, content='\U0001F6AB Alias name cannot be empty.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        room_name = ''
        temp_room_obj = None
        for guild_temp_rooms in self.temp_rooms.values():
            for temp_channel in guild_temp_rooms.values():
                if temp_channel.channel.id == channel_obj.id:
                    temp_room_obj = temp_channel
                    room_name = temp_channel.room_name
                    break
            if temp_room_obj:
                break
        if not channel_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid channel from input: {channel}.')
        is_owner_or_dev, _ = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev:
            async with ctx.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(
                    'SELECT coordinator_room_names, coordinator_channel_ids FROM users WHERE discord_snowflake=$1',
                    ctx.author.id
                )
            allowed_rooms = row.get('coordinator_room_names') or []
            allowed_channels = row.get('coordinator_channel_ids') or []
            if channel_obj.name not in allowed_rooms and channel_obj.id not in allowed_channels:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command in `{channel_obj.name}`.')
        async with ctx.bot.db_pool.acquire() as conn:
            existing_alias = await conn.fetchrow('''
                SELECT guild_id, channel_id, role_id, room_name
                FROM command_aliases
                WHERE guild_id=$1 AND alias_type=$2 AND alias_name=$3 AND room_name = $4
            ''', ctx.guild.id, alias_type, alias_name, room_name)
            if existing_alias:
                existing_channel = ctx.guild.get_channel(existing_alias['channel_id']) if existing_alias['channel_id'] else None
                existing_role = ctx.guild.get_role(existing_alias['role_id']) if existing_alias['role_id'] else None
                existing_room = existing_alias.get('room_name')
                mention = existing_channel.mention if existing_channel else (existing_role.mention if existing_role else (existing_room or 'unknown'))
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Alias `{alias_name}` ({alias_type}) already exists and is set to {mention}.')
        if self.bot.get_command(alias_name):
            return await self.handler.send_message(ctx, content=f'\U0001F6AB A command named `{alias_name}` already exists.')
        if alias_type in ('role','unrole'):
            if not role:
                return await self.handler.send_message(ctx, content='\U0001F6AB Role ID is required for role/unrole aliases.')
            try:
                role_id = int(role.replace('<@&','').replace('>',''))
            except ValueError:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid role ID: {role}')
            async with ctx.bot.db_pool.acquire() as conn:
                await conn.execute('''
                    INSERT INTO command_aliases (
                        guild_id, alias_type, alias_name, channel_id, role_id, room_name
                    )
                    VALUES ($1,$2,$3,$4,$5,$6)
                ''', ctx.guild.id, alias_type, alias_name, channel_obj.id, role_id, room_name)
        else:
            async with ctx.bot.db_pool.acquire() as conn:
                await conn.execute('''
                    INSERT INTO command_aliases (
                        guild_id, alias_type, alias_name, channel_id, role_id, room_name
                    )
                    VALUES ($1,$2,$3,$4,$5,$6)
                ''', ctx.guild.id, alias_type, alias_name, channel_obj.id, None, room_name)
        async with ctx.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'create_alias', None, ctx.author.id, ctx.guild.id, channel_obj.id, f'Created an alias: {alias_name}')
        if alias_type in ('role','unrole') and is_owner_or_dev:
            self.bot.command_aliases.setdefault(ctx.guild.id, self.bot.command_aliases.default_factory()).setdefault('role_aliases', {}).setdefault(alias_type, {})[alias_name] = {'channel_id': int(channel_obj.id), 'role_id': int(role_id)}
        elif room_name != '':
            self.bot.command_aliases.setdefault(ctx.guild.id, self.bot.command_aliases.default_factory()).setdefault('temp_room_aliases', {}).setdefault(alias_type, {})[alias_name] = {'room_name': room_name, 'channel_id': int(channel_obj.id)}
        else:
            self.bot.command_aliases.setdefault(ctx.guild.id, self.bot.command_aliases.default_factory()).setdefault('channel_aliases', {}).setdefault(alias_type, {})[alias_name] = int(channel_obj.id)
        cmd = None
        if alias_type == 'ban': cmd = self.create_ban_alias(alias_name)
        elif alias_type == 'flag': cmd = self.create_flag_alias(alias_name)
        elif alias_type == 'unflag': cmd = self.create_unflag_alias(alias_name)
        elif alias_type == 'mute': cmd = self.create_voice_mute_alias(alias_name)
        elif alias_type == 'tmute': cmd = self.create_text_mute_alias(alias_name)
        elif alias_type == 'unban': cmd = self.create_unban_alias(alias_name)
        elif alias_type == 'cow':
            if channel_obj.id != 1222056499959042108: return await self.handler.send_message(ctx, content=f'\U0001F6AB {channel_obj.mention} has no vegan association.')
            cmd = self.create_cow_alias(alias_name)
        elif alias_type == 'uncow':
            if channel_obj.id != 1222056499959042108: return await self.handler.send_message(ctx, content=f'\U0001F6AB {channel_obj.mention} has no vegan association.')
            cmd = self.create_uncow_alias(alias_name)
        elif alias_type == 'unmute': cmd = self.create_unmute_alias(alias_name)
        elif alias_type == 'untmute': cmd = self.create_untextmute_alias(alias_name)
        elif alias_type == 'role': cmd = self.create_role_alias(alias_name)
        elif alias_type == 'unrole': cmd = self.create_unrole_alias(alias_name)
        if cmd: self.bot.add_command(cmd)
        mention = ctx.guild.get_role(int(role_id)).mention if alias_type in ('role','unrole') else channel_obj.mention
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Alias `{alias_name}` ({alias_type}) set to {mention}.', allowed_mentions=discord.AllowedMentions.none())
                
    @app_commands.command(name='backup', description='Creates a backup of the database and uploads it')
    @is_owner_developer_app_predicator()
    async def app_backup(
        self,
        interaction: discord.Interaction
    ):
        try:
            backup_dir = self.setup_backup_directory('/app/backups')
            backup_file = self.perform_backup(
                db_user=os.getenv('POSTGRES_USER'),
                db_name=os.getenv('POSTGRES_DATABASE'),
                db_host=os.getenv('POSTGRES_HOST'),
                db_password=os.getenv('POSTGRES_PASSWORD'),
                backup_dir=backup_dir
            )
            if backup_file:
                await self.send(interaction, file=discord.File(backup_file))
            else:
                await self.send(interaction, content=f'\U0001F6AB Failed to create backup.')
        except Exception as e:
            logger.warning(f'Database error occurred: {e}')
            
    @commands.command(name='backup', help='Creates a backup of the database and uploads it')
    @is_owner_developer_predicator()
    async def backup(
        self,
        ctx: commands.Context
    ):
        try:
            backup_dir = self.setup_backup_directory('/app/backups')
            backup_file = self.perform_backup(
                db_user=os.getenv('POSTGRES_USER'),
                db_name=os.getenv('POSTGRES_DATABASE'),
                db_host=os.getenv('POSTGRES_HOST'),
                db_password=os.getenv('POSTGRES_PASSWORD'),
                backup_dir=backup_dir
            )
            if backup_file:
                await self.handler.send_message(ctx, file=discord.File(backup_file))
            else:
                await selfhandler.send_message(ctx, content=f'\U0001F6AB Failed to create backup.')
        except Exception as e:
            logger.warning(f'Database error occurred: {e}')
    
    @app_commands.command(name='cap', description='Set a duration limit for bans, mutes and text mutes.')
    @is_owner_developer_app_predicator()
    @app_commands.describe(
        channel='Tag a channel or include its snowflake ID',
        moderation_type='One of: `mute`, `ban`, `tmute`',
        duration='(+|-)duration(m|h|d), 0=permanent, default=24h'
    )
    async def cap_app_command(
        self,
        interaction: discord.Interaction,
        channel: str = None,
        moderation_type: str = None,
        duration: str = '24'
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        room_name = ''
        temp_room_obj = None
        is_temp_room = False
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        is_temp_room = True
                        break
                if temp_room_obj:
                    break        
        valid_types = {'mute', 'ban', 'tmute'}
        if moderation_type not in valid_types:
            return await self.send(interaction, content=f'\U0001F6AB Invalid moderation type. Must be one of: {", ".join(valid_types)}')
        expires_at, duration_str = self.parse_duration(duration)
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4',
                interaction.guild.id, channel_obj.id, moderation_type, room_name
            )
            original_duration = row['duration'] if row else None
            await conn.execute('''
                INSERT INTO active_caps (guild_id, channel_id, moderation_type, duration, room_name)
                VALUES ($1, $2, $3, $4, $5)
                ON CONFLICT (guild_id, channel_id, moderation_type, room_name)
                DO UPDATE SET duration = EXCLUDED.duration
            ''', interaction.guild.id, channel_obj.id, moderation_type, duration, room_name)
        if is_temp_room:
            room_display = room_name
        else:
            room_display = channel_obj.mention
        if original_duration:
            msg = f'{self.get_random_emoji()} Cap changed on {room_display} for {moderation_type} from {original_duration} to {duration_str}.'
        else:
            msg = f'{self.get_random_emoji()} Cap set on {room_display} for {moderation_type} for {duration_str}.'
        await self.send(interaction, content=msg, allowed_mentions=discord.AllowedMentions.none())
    
    @commands.command(name='cap', help='Set a duration limit for bans, mutes and text mutes.')
    @is_owner_developer_predicator()
    async def cap_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        moderation_type: Optional[str] = commands.parameter(default=None, description='One of: `mute`, `ban`, `tmute`'),
        *,
        duration: Optional[str] = commands.parameter(default='24', description='Options: (+|-)duration(m|h|d) \n 0 - permanent / 24h - default')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        room_name = ''
        temp_room_obj = None
        is_temp_room = False
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        is_temp_room = True
                        break
                if temp_room_obj:
                    break
        valid_types = {'mute', 'ban', 'tmute'}
        if moderation_type not in valid_types:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid moderation type. Must be one of: {", ".join(valid_types)}')
        expires_at, duration_str = self.parse_duration(duration)
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4',
                ctx.guild.id, channel_obj.id, moderation_type, room_name
            )
            original_duration = row['duration'] if row else None
            await conn.execute('''
                INSERT INTO active_caps (guild_id, channel_id, moderation_type, duration, room_name)
                VALUES ($1, $2, $3, $4, $5)
                ON CONFLICT (guild_id, channel_id, moderation_type, room_name)
                DO UPDATE SET duration = EXCLUDED.duration
            ''', ctx.guild.id, channel_obj.id, moderation_type, duration, room_name)
        if is_temp_room:
            room_display = room_name
        else:
            room_display = channel_obj.mention
        if original_duration:
            msg = f'{self.get_random_emoji()} Cap changed on {room_display} for {moderation_type} from {original_duration} to {duration_str}.'
        else:
            msg = f'{self.get_random_emoji()} Cap set on {room_display} for {moderation_type} for {duration_str}.'
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())
    
    @app_commands.command(name='clear', description='Removes a specific channel ID from all users, including temp-room associations.')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_app_predicator()
    async def clear_channel_access_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $1::bigint),
                    moderator_channel_ids = array_remove(moderator_channel_ids, $1::bigint),
                    updated_at = NOW()
            ''', channel_obj.id)
            await conn.execute('''
                UPDATE users
                SET coordinator_room_names = (
                        SELECT array_agg(name) 
                        FROM unnest(coordinator_room_names) AS name
                        WHERE EXISTS (
                            SELECT 1 FROM command_aliases 
                            WHERE channel_id = $1 AND room_name = name
                        )
                    ),
                    moderator_room_names = (
                        SELECT array_agg(name) 
                        FROM unnest(moderator_room_names) AS name
                        WHERE EXISTS (
                            SELECT 1 FROM command_aliases 
                            WHERE channel_id = $1 AND room_name = name
                        )
                    ),
                    updated_at = NOW()
            ''', channel_obj.id)
        await self.send(interaction, content=f'{self.get_random_emoji()} Removed channel ID `{channel_obj.id}` from all users\' coordinator and moderator access, including temp-room associations.', allowed_mentions=discord.AllowedMentions.none())

   # CHECK
    @commands.command(name='clear', help='Removes a specific channel ID from all users, including temp-room associations and all related records.')
    @is_owner_predicator()
    async def clear_channel_access_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $1::bigint),
                    moderator_channel_ids = array_remove(moderator_channel_ids, $1::bigint),
                    updated_at = NOW()
            ''', channel_obj.id)
            await conn.execute('''
                UPDATE users
                SET coordinator_room_names = ARRAY(
                    SELECT name
                    FROM unnest(coordinator_room_names) AS name
                    WHERE name NOT IN (
                        SELECT room_name FROM temporary_rooms WHERE room_snowflake = $1
                    )
                ),
                moderator_room_names = ARRAY(
                    SELECT name
                    FROM unnest(moderator_room_names) AS name
                    WHERE name NOT IN (
                        SELECT room_name FROM temporary_rooms WHERE room_snowflake = $1
                    )
                ),
                updated_at = NOW()
            ''', channel_obj.id)
            tables_with_channel_id = [
                'command_aliases',
                'active_bans',
                'active_text_mutes',
                'active_voice_mutes',
                'active_stages',
                'stage_coordinators',
                'active_caps'
            ]
            for table in tables_with_channel_id:
                await conn.execute(f'''
                    DELETE FROM {table}
                    WHERE channel_id = $1 OR room_name IN (
                        SELECT room_name FROM temporary_rooms WHERE room_snowflake = $1
                    )
                ''', channel_obj.id)
            await conn.execute('DELETE FROM temporary_rooms WHERE room_snowflake = $1', channel_obj.id)
        guild_aliases = self.bot.command_aliases.setdefault(channel_obj.guild.id, self.bot.command_aliases.default_factory())
        temp_aliases = guild_aliases.get('temp_room_aliases', {})
        for alias_type, aliases in temp_aliases.items():
            for alias_name, data in list(aliases.items()):
                if data.get('channel_id') == channel_obj.id:
                    del aliases[alias_name]
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Removed channel {channel_obj.mention} from all users and deleted all associated records.', allowed_mentions=discord.AllowedMentions.none())
    
    @app_commands.command(name='coord', description='Grants coordinator access for a specific voice channel.')
    @is_owner_developer_app_predicator()
    @app_commands.describe(member='Tag a member or include their snowflake ID', channel='Tag a channel or include its snowflake ID')
    async def create_coordinator_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member_app(interaction, member)
        if not member_obj:
            return await self.send(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        is_owner_or_dev, _ = await check_owner_dev_coord_app(interaction, channel_obj)
        if member_obj.bot:
            return await self.send(interaction, content='\U0001F6AB You cannot make the bot a coordinator.')
        highest_role, success = await check_block_app(interaction, member_obj, channel_obj)
        if not success:
            return await self.send(interaction, content=f'\U0001F6AB You are not allowed to make this `{highest_role}` a coordinator because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if room_name != '':
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, coordinator_room_names)
                    VALUES ($1, ARRAY[$2]::TEXT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET coordinator_room_names = (
                        SELECT ARRAY(SELECT DISTINCT unnest(
                            COALESCE(users.coordinator_room_names, ARRAY[]::TEXT[]) || ARRAY[$2]
                        ))
                    ),
                    updated_at = NOW()
                ''', member_obj.id, room_name)
            else:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, coordinator_channel_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET coordinator_channel_ids = (
                        SELECT ARRAY(SELECT DISTINCT unnest(
                            COALESCE(users.coordinator_channel_ids, ARRAY[]::BIGINT[]) || ARRAY[$2]::BIGINT[]
                        ))
                    ),
                    updated_at = NOW()
                ''', member_obj.id, channel_obj.id)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'create_coordinator', member_obj.id, interaction.user.id, interaction.guild.id, channel_obj.id, 'Created a coordinator')
        return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has been granted coordinator rights in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='coord', help='Grants coordinator access for a specific voice channel.')
    @is_owner_developer_predicator()
    async def create_coordinator_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        room_name = ''
        temp_room_obj = None
        for guild_temp_rooms in self.temp_rooms.values():
            for temp_channel in guild_temp_rooms.values():
                if temp_channel.channel.id == channel_obj.id:
                    temp_room_obj = temp_channel
                    room_name = temp_channel.room_name
                    break
            if temp_room_obj:
                break
        is_owner_or_dev, _ = await check_owner_dev_coord(ctx, channel_obj)
        if member_obj.bot:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot a coordinator.')
        highest_role, success = await check_block(ctx, member_obj, channel_obj)
        if not success:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to make this `{highest_role}` a coordinator because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if room_name != '':
                await conn.execute('''
                INSERT INTO users (discord_snowflake, coordinator_room_names)
                VALUES ($1, ARRAY[$2]::TEXT[])
                ON CONFLICT (discord_snowflake) DO UPDATE SET
                    coordinator_room_names = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(
                                COALESCE(users.coordinator_room_names, ARRAY[]::TEXT[]) || ARRAY[$2]
                            )
                        )
                    ),
                    updated_at = NOW();
                ''', member_obj.id, room_name)
            else:
                await conn.execute('''
                    INSERT INTO users AS u (discord_snowflake, coordinator_channel_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE SET
                        coordinator_channel_ids = (
                            SELECT ARRAY(
                                SELECT DISTINCT unnest(
                                    COALESCE(u.coordinator_channel_ids, ARRAY[]::BIGINT[]) || ARRAY[$2]::BIGINT[]
                                )
                            )
                        ),
                        updated_at = NOW();
                ''', member_obj.id, channel_obj.id)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'create_coordinator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Created a coordinator')
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been granted coordinator rights in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
    
    @app_commands.command(name='del', description='Delete a message by ID (only if you are coordinator/moderator of that temp room).')
    @app_commands.describe(message_id='Message snowflake ID', channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_coordinator_moderator_app_predicator(None)
    async def delete_message_app_command(
        self,
        interaction: discord.Interaction,
        message_id: Optional[str],
        channel: Optional[str]
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.text:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        try:
            msg = await channel_obj.fetch_message(int(message_id))
        except:
            logger.warning('No message with this ID exists.')
        if not msg:
            return await self.send(interaction, content=f'\U0001F6AB No message with ID `{message_id}` found in `{channel_obj.name}`.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        if not (is_owner_or_dev or is_mod_or_coord):
            return await self.send(interaction, content=f'\U0001F6AB You have insufficient privileges in {channel_obj.mention} to delete messages.')
        try:
            await msg.delete()
            return await self.send(interaction, content=f'{self.get_random_emoji()} Message `{message_id}` deleted successfully.')
        except discord.Forbidden:
            logger.warning('Missing permissions to delete the message.')
        return await self.send(interaction, content='\U0001F6AB Failed to delete the message.')

    @commands.command(name='del', help='Delete a message by ID (only if you are coordinator/moderator of that temp room).')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        message_id: Optional[int] = commands.parameter(default=None, description='Message snowflake'),
        *,
        channel: Optional[str] = commands.parameter(default=None, description='Channel or snowflake')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        try:
            msg = await channel_obj.fetch_message(message_id)
            if channel_obj.type != discord.ChannelType.text:
                return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        except:
            logger.warning('No message with this ID exists.')
        if not msg:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB No message with ID `{message_id}` found in `{channel_obj.name}`.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not (is_owner_or_dev or is_mod_or_coord):
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You have insufficient privileges in `{channel_obj.name}` to delete messages.')
        try:
            await msg.delete()
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Message `{message_id}` deleted successfully.')
        except discord.Forbidden:
            logger.warning('Missing permissions to delete the message.')
        return await self.handler.send_message(ctx, content='\U0001F6AB Failed to delete the message.')
        
    @app_commands.command(name='dev', description='Elevates a user\'s permissions to a bot developer.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_app_predicator()
    async def create_developer_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member_app(interaction, member)
        if member_obj.bot:
            return await self.send(interaction, content='\U0001F6AB You cannot make the bot a developer.')
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        highest_role, success = await check_block_app(interaction, member_obj, None)
        if not success:
            return await self.send(interaction, content=f'\U0001F6AB You are not allowed to make this `{highest_role}` a developer because they are a higher/or equivalent role than you in {interaction.guild.name}.')
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
            ''', member_obj.id, interaction.guild.id)
        return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has been granted developer rights in {interaction.guild.name}.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='dev', help='Elevates a user\'s permissions to a bot developer.')
    @is_owner_predicator()
    async def create_developer_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
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
        await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been granted developer rights in {ctx.guild.name}.', allowed_mentions=discord.AllowedMentions.none())

    @app_commands.command(name='mod', description="Elevates a user's permission to VC moderator for a specific channel.")
    @app_commands.describe(member='Tag a member or include their snowflake ID', channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_coordinator_app_predicator(None)
    async def create_moderator_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member_app(interaction, member)
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await self.send(interaction, content='\U0001F6AB You cannot make the bot a moderator.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        room_name = ''
        temp_room_obj = None
        for guild_temp_rooms in self.temp_rooms.values():
            for temp_channel in guild_temp_rooms.values():
                if temp_channel.channel.id == channel_obj.id:
                    temp_room_obj = temp_channel
                    room_name = temp_channel.room_name
                    break
            if temp_room_obj:
                break
        is_owner_or_dev, is_coord = await check_owner_dev_coord_app(interaction, channel_obj)
        if not is_owner_or_dev and not is_coord:
            return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use this command (`mod`) in {channel_obj.mention}')
        highest_role, success = await check_block_app(interaction, member_obj, channel_obj)
        if not success:
            return await self.send(interaction, content=f'\U0001F6AB You are not allowed to make this `{highest_role}` a moderator because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        if not is_owner_or_dev:
            async with self.bot.db_pool.acquire() as conn:
                if room_name != '':
                    row = await conn.fetchrow('''
                        SELECT 1 FROM users
                        WHERE discord_snowflake=$1 AND $2 = ANY(coordinator_room_names)
                    ''', interaction.user.id, room_name)
                else:
                    row = await conn.fetchrow('''
                        SELECT 1 FROM users
                        WHERE discord_snowflake=$1 AND $2 = ANY(coordinator_channel_ids)
                    ''', interaction.user.id, channel_obj.id)
                if not row:
                    return await self.send(interaction, content=f'\U0001F6AB You are not a coordinator in {channel_obj.mention} and cannot assign moderators there.')
        async with self.bot.db_pool.acquire() as conn:
            if room_name != '':
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, moderator_room_names)
                    VALUES ($1, ARRAY[$2]::TEXT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET moderator_room_names = (
                        SELECT ARRAY(SELECT DISTINCT unnest(
                            COALESCE(users.moderator_room_names, ARRAY[]::TEXT[]) || ARRAY[$2]
                        ))
                    ),
                    updated_at = NOW()
                ''', member_obj.id, room_name)
            else:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, moderator_channel_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET moderator_channel_ids = (
                        SELECT ARRAY(SELECT DISTINCT unnest(
                            COALESCE(users.moderator_channel_ids, ARRAY[]::BIGINT[]) || ARRAY[$2]::BIGINT[]
                        ))
                    ),
                    updated_at = NOW()
                ''', member_obj.id, channel_obj.id)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'create_moderator', member_obj.id, interaction.user.id, interaction.guild.id, channel_obj.id, 'Created a moderator')
        return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has been granted moderator rights in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                
    @commands.command(name='mod', help='Elevates a user\'s permission to VC moderator for a specific channel.')
    @is_owner_developer_coordinator_predicator(None)
    async def create_moderator_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot a moderator.')
        channel_obj = await self.resolve_channel(ctx, channel)
        room_name = ''
        temp_room_obj = None
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        for guild_temp_rooms in self.temp_rooms.values():
            for temp_channel in guild_temp_rooms.values():
                if temp_channel.channel.id == channel_obj.id:
                    temp_room_obj = temp_channel
                    room_name = temp_channel.room_name
                    break
            if temp_room_obj:
                break
        is_owner_or_dev, is_coord = await check_owner_dev_coord(ctx, channel_obj)
        if not is_owner_or_dev and not is_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`mod`) in {channel_obj.mention}')
        highest_role, success = await check_block(ctx, member_obj, channel_obj)
        if not success:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to make this `{highest_role}` a moderator because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        if not is_owner_or_dev:
            async with self.bot.db_pool.acquire() as conn:
                if room_name != '':
                    coordinator_row = await conn.fetchrow('''
                        SELECT 1 FROM users
                        WHERE discord_snowflake=$1 AND $2 = ANY(coordinator_room_names)
                    ''', ctx.author.id, room_name)
                else:
                    coordinator_row = await conn.fetchrow('''
                        SELECT 1 FROM users
                        WHERE discord_snowflake=$1 AND $2 = ANY(coordinator_channel_ids)
                    ''', ctx.author.id, channel_obj.id)
                if not coordinator_row:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not a coordinator in {channel_obj.mention} and cannot assign moderators there.')
        async with self.bot.db_pool.acquire() as conn:
            if room_name != '':
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, moderator_room_names)
                    VALUES ($1, ARRAY[$2]::TEXT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET moderator_room_names = (
                        SELECT ARRAY(SELECT DISTINCT unnest(
                            COALESCE(users.moderator_room_names, ARRAY[]::TEXT[]) || ARRAY[$2]
                        ))
                    ),
                    updated_at = NOW()
                ''', member_obj.id, room_name)
            else:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, moderator_channel_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET moderator_channel_ids = (
                        SELECT ARRAY(SELECT DISTINCT unnest(
                            COALESCE(users.moderator_channel_ids, ARRAY[]::BIGINT[]) || ARRAY[$2]::BIGINT[]
                        ))
                    ),
                    updated_at = NOW()
                ''', member_obj.id, channel_obj.id)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'create_moderator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Created a moderator')
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been granted moderator rights in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
        
    @app_commands.command(name='bans', description='Lists ban statistics.')
    @app_commands.describe(target='"all", channel name/ID/mention, or user mention/ID')
    @is_owner_developer_coordinator_moderator_role_app_predicator(None)
    async def list_bans_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in interaction.user.roles)
        member_obj = await self.resolve_member_app(interaction, target)
        if member_obj and member_obj.bot:
            return await self.send(interaction, content='\U0001F6AB You cannot list bans on the bot.')
        channel_obj = await self.resolve_channel_app(interaction, target)
        if member_obj:
            target = None
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall_app(interaction)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use command (`bans`) in {channel_obj.mention}.')
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.send(interaction, content='\U0001F6AB Only owners or developers can list all bans across the server.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(
                    '''SELECT discord_snowflake, room_name, channel_id, expires_at, reason
                       FROM active_bans
                       WHERE guild_id=$1
                       ORDER BY channel_id, room_name, expires_at NULLS LAST''', interaction.guild.id)
            if not rows:
                return await self.send(interaction, content=f'\U0001F6AB No active bans found in {interaction.guild.name}.')
            grouped = defaultdict(list)
            for row in rows: grouped[row['channel_id']].append(row)
            pages = []
            for ch_id, records in grouped.items():
                ch = interaction.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                embed = discord.Embed(
                    title=f'\u26D4 Ban records for {ch_name}',
                    color=discord.Color.red()
                )
                for record in records:
                    user = interaction.guild.get_member(record['discord_snowflake'])
                    reason = record['reason'] or 'No reason provided'
                    if record['expires_at'] is None:
                        duration_str = 'Permanent'
                    else:
                        now = discord.utils.utcnow()
                        delta = record['expires_at'] - now
                        if delta.total_seconds() <= 0: duration_str = 'Expired'
                        else:
                            days, seconds = delta.days, delta.seconds
                            hours = seconds // 3600
                            minutes = (seconds % 3600) // 60
                            duration_str = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                    mention = user.name if user else f'`{record["discord_snowflake"]}`'
                    embed.add_field(name='User', value=f'{mention}\nReason: {reason}\nDuration: {duration_str}', inline=False)
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
        if member_obj:
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch(
                    '''SELECT channel_id, expires_at, reason
                       FROM active_bans
                       WHERE guild_id=$1 AND discord_snowflake=$2 AND room_name = $3''', interaction.guild.id, member_obj.id, room_name)
            bans = [b for b in bans if interaction.guild.get_channel(b['channel_id'])]
            if not bans: return await self.send(interaction, content=f'\U0001F6AB {member_obj.mention} is not banned in any channels.', allowed_mentions=discord.AllowedMentions.none())
            embed = discord.Embed(
                title=f'\u26D4 Ban records for {member_obj.display_name}',
                color=discord.Color.red()
            )
            for record in bans:
                ch_obj = interaction.guild.get_channel(record['channel_id'])
                channel_mention = ch_obj.mention if ch_obj else f'Channel ID `{record["channel_id"]}`'
                reason = record['reason'] or 'No reason provided'
                if record['expires_at'] is None:
                    duration_str = 'Permanent'
                else:
                    now = discord.utils.utcnow()
                    delta = record['expires_at'] - now
                    if delta.total_seconds() <= 0:
                        duration_str = 'Expired'
                    else:
                        days, seconds = delta.days, delta.seconds
                        hours = seconds // 3600
                        minutes = (seconds % 3600) // 60
                        duration_str = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                embed.add_field(name=channel_mention, value=f'Reason: {reason}\nDuration: {duration_str}', inline=False)
            return await self.send(interaction, embed=embed)
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch(
                    '''SELECT discord_snowflake, expires_at, reason
                       FROM active_bans
                       WHERE guild_id=$1 AND channel_id=$2 AND room_name = $3
                       ORDER BY expires_at NULLS LAST''', interaction.guild.id, channel_obj.id, room_name)
            if not bans:
                return await self.send(interaction, content=f'\U0001F6AB No active bans found for {channel_obj.mention}.')
            lines = []
            for record in bans:
                uid = record['discord_snowflake']
                member_obj = interaction.guild.get_member(uid)
                if not member_obj: continue
                name = member_obj.display_name
                if record['expires_at'] is None:
                    time_left = 'Permanent'
                else:
                    now = discord.utils.utcnow()
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
                return await self.send(interaction, content=f'\U0001F6AB No active bans for users currently in {interaction.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'\u26D4 Ban records for {channel_obj.name}',
                    description='\n'.join(chunk),
                    color=discord.Color.red()
                )
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
        return await self.send(interaction, content='\U0001F6AB You must specify a member, a text channel or use "all".')
        
    @commands.command(name='bans', description='Lists ban statistics.')
    @is_owner_developer_coordinator_moderator_role_predicator(None)
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in ctx.author.roles)
        member_obj = await self.resolve_member(ctx, target)
        if member_obj and member_obj.bot:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot list bans on the bot.')
        channel_obj = await self.resolve_channel(ctx, target)
        if member_obj:
            target = None
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall(ctx)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use command (`bans`) for {channel_obj.mention}.')
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners or developers can list all bans across the server.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''SELECT discord_snowflake, channel_id, room_name, expires_at, reason FROM active_bans WHERE guild_id = $1 ORDER BY channel_id, room_name, expires_at NULLS LAST''', ctx.guild.id)
            if not rows:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No active bans found in {ctx.guild.name}.')
            grouped = defaultdict(list)
            for row in rows:
                grouped[row['channel_id']].append(row)
            pages = []
            for ch_id, records in grouped.items():
                ch = ctx.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                embed = discord.Embed(
                    title = f'\u26D4 Ban records for {ch_name}',
                    color = discord.Color.red()
                )
                for record in records:
                    user = ctx.guild.get_member(record['discord_snowflake'])
                    reason = record['reason'] or 'No reason provided'
                    if record['expires_at'] is None:
                        duration_str='Permanent'
                    else:
                        now = discord.utils.utcnow()
                        delta = record['expires_at'] - now
                        if delta.total_seconds() <= 0:
                            duration_str = 'Expired'
                        else:
                            days, seconds = delta.days,delta.seconds
                            hours = seconds // 3600
                            minutes = (seconds % 3600) // 60
                            duration_str = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                    mention = user.name if user else f'`{record["discord_snowflake"]}`'
                    embed.add_field(name='User', value=f'{mention}\nReason: {reason}\nDuration: {duration_str}', inline=False)
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        elif member_obj:
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch('''SELECT channel_id,expires_at,reason FROM active_bans WHERE guild_id=$1 AND discord_snowflake=$2 AND room_name = $3''', ctx.guild.id, member_obj.id, room_name)
            bans = [b for b in bans if ctx.guild.get_channel(b['channel_id'])]
            if not bans:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not banned in any channels.', allowed_mentions=discord.AllowedMentions.none())
            embed = discord.Embed(
                title=f'\u26D4 Ban records for {member_obj.display_name}',
                color=discord.Color.red()
            )
            for record in bans:
                channel_obj = ctx.guild.get_channel(record['channel_id'])
                channel_mention = channel_obj.mention if channel_obj else f'Channel ID `{record["channel_id"]}`'
                reason = record['reason'] or 'No reason provided'
                if record['expires_at'] is None:
                    duration_str = 'Permanent'
                else:
                    now = discord.utils.utcnow()
                    delta = record['expires_at'] - now
                    if delta.total_seconds() <= 0:
                        duration_str = 'Expired'
                    else:
                        days, seconds = delta.days, delta.seconds
                        hours = seconds // 3600
                        minutes = (seconds % 3600) // 60
                        duration_str = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                embed.add_field(name = channel_mention, value = f'Reason: {reason}\nDuration: {duration_str}', inline=False)
            return await self.handler.send_message(ctx, embed=embed)
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
               return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch('''SELECT discord_snowflake,expires_at,reason FROM active_bans WHERE guild_id=$1 AND channel_id=$2 AND room_name = $3 ORDER BY expires_at NULLS LAST''',ctx.guild.id,channel_obj.id, room_name)
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
                    now = discord.utils.utcnow()
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
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No active bans for users currently in {ctx.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i+chunk_size]
                embed = discord.Embed(
                    title = f'\u26D4 Ban records for {channel_obj.name}',
                    description = '\n'.join(chunk),
                    color = discord.Color.red()
                )
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F6AB You must specify a member, a text channel or use "all".')

    @app_commands.command(name='caps', description='List active caps for a channel or all channels if "all" is provided.')
    @is_owner_developer_coordinator_moderator_app_predicator(None)
    @app_commands.describe(target='"all", channel name/ID/mention')
    async def list_caps_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in interaction.user.roles)
        is_owner_or_dev, _ = await check_owner_dev_coord_app(interaction, None)
        if target and target.lower() == "all":
            if not is_owner_or_dev:
                return await self.send(interaction, content='\U0001F6AB Only owners or developers can list all caps.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(
                    "SELECT channel_id, moderation_type, duration FROM active_caps WHERE guild_id=$1",
                    interaction.guild.id
                )
            if not rows:
                return await self.send(interaction, content='\U0001F6AB No caps found server-wide.')
            lines = []
            for row in rows:
                ch = interaction.guild.get_channel(row["channel_id"])
                ch_name = ch.mention if ch else f'Channel ID `{row["channel_id"]}`'
                lines.append(f'**{row["moderation_type"]} in {ch_name}** → `{row["duration"]}`')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                embed = discord.Embed(
                    title="All Active Caps in Server",
                    description="\n".join(lines[i:i+chunk_size]),
                    color=discord.Color.red()
                )
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
        channel_obj = await self.resolve_channel_app(interaction, target)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.send(interaction, content=f'\U0001F6AB You do not have permission to view caps for {channel_obj.mention}.')
        caps = await self.get_caps_for_channel(interaction.guild.id, channel_obj.id)
        if not caps:
            return await self.send(interaction, content='\U0001F6AB No caps found for this channel.')
        lines = [
            f'**{moderation_type} in {channel_obj.mention}** → `{duration}`'
            for moderation_type, duration in caps
        ]
        embed = discord.Embed(
            title=f"Active Caps for {channel_obj.mention}",
            description="\n".join(lines),
            color=discord.Color.red()
        )
        return await self.send(interaction, embed=embed)
        
    @commands.command(name='caps', help='List active caps for a channel or all channels if "all" is provided.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in ctx.author.roles)
        if target and target.lower() == "all":
            is_owner_or_dev, _ = await check_owner_dev_coord_mod(ctx, None)
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners or developers can list all caps.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(
                    "SELECT channel_id, moderation_type, duration FROM active_caps WHERE guild_id=$1",
                    ctx.guild.id
                )
            if not rows:
                return await self.handler.send_message(ctx, content='\U0001F6AB No caps found server-wide.')
            lines = []
            for row in rows:
                ch = ctx.guild.get_channel(row["channel_id"])
                ch_name = ch.mention if ch else f'Channel ID `{row["channel_id"]}`'
                lines.append(f'**{row["moderation_type"]} in {ch_name}** → `{row["duration"]}`')
            embed = discord.Embed(
                title='All Active Caps in Server',
                description='\n'.join(lines),
                color=discord.Color.red()
            )
            return await self.handler.send_message(ctx, embed=embed)
        channel_obj = await self.resolve_channel(ctx, target)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to list caps in {channel_obj.mention}.')
        caps = await self.get_caps_for_channel(ctx.guild.id, channel_obj.id)
        if not caps:
            return await self.handler.send_message(ctx, content='\U0001F6AB No caps found for this channel.')
        lines = [f'**{mtype} in {channel_obj.mention}** → `{duration}`' for mtype, duration in caps]
        embed = discord.Embed(
            title=f'Active Caps for {channel_obj.mention}',
            description='\n'.join(lines),
            color=discord.Color.red()
        )
        await self.handler.send_message(ctx, embed=embed)
    
    @app_commands.command(name='cmds', description='List command aliases routed to a specific channel, temp room, or all channels if "all" is provided.')
    @app_commands.describe(target='"all", channel name/ID/mention')
    @is_owner_developer_coordinator_moderator_app_predicator(None)
    async def list_commands_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if interaction.guild is None:
            return await self.send(interaction, content="This command must be used in a server.")
        aliases = self.bot.command_aliases.get(interaction.guild.id, self.bot.command_aliases.default_factory())
        lines = []
        channel_obj = await self.resolve_channel_app(interaction, target)
        temp_room_obj = None
        found_aliases = False
        for guild_temp_rooms in self.temp_rooms.values():
            for temp_channel in guild_temp_rooms.values():
                if temp_channel.room_name.lower() == channel_obj.name.lower():
                    temp_room_obj = temp_channel
                    break
            if temp_room_obj:
                break
        found_aliases = False
        if target and target.lower() == 'all':
            is_owner_or_dev, _ = await check_owner_dev_coord_mod_app(interaction, None)
            if not is_owner_or_dev:
                return await self.send(interaction, content='\U0001F6AB Only owners or developers can list all aliases across the server.')
            for kind, type_map in aliases.get('channel_aliases', {}).items():
                grouped_by_channel = defaultdict(list)
                for name, cid in type_map.items():
                    grouped_by_channel[cid].append(name)
                for ch_id, names in grouped_by_channel.items():
                    ch = interaction.guild.get_channel(ch_id)
                    ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                    lines.append(f'**{kind.capitalize()} in {ch_name}**')
                    lines.extend(f'`{name}`' for name in names)
                    found_aliases = True
            for kind, type_map in aliases.get('role_aliases', {}).items():
                grouped_by_channel = defaultdict(list)
                for name, data in type_map.items():
                    if isinstance(data, dict) and 'channel_id' in data:
                        grouped_by_channel[data['channel_id']].append((name, data.get('role_id')))
                for ch_id, entries in grouped_by_channel.items():
                    ch = ctx.guild.get_channel(ch_id)
                    ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                    lines.append(f'**{kind.capitalize()} Role Aliases**')
                    for name, rid in entries:
                        role = ctx.guild.get_role(rid)
                        mention = role.mention if role else f'<@&{rid}>'
                        lines.append(f'`{name}` → {mention}')
                    found_aliases = True
            for alias_type, room_map in aliases.get('temp_room_aliases', {}).items():
                for alias_name, data in room_map.items():
                    room_name = data.get('room_name') or '(unnamed room)'
                    lines.append(f'**{alias_type.capitalize()}**')
                    lines.append(f'`{alias_name}`')
                    found_aliases = True
        else:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
            for kind, type_map in aliases.get('channel_aliases', {}).items():
                channel_entries = {name: cid for name, cid in type_map.items() if cid == channel_obj.id}
                if channel_entries:
                    found_aliases = True
                    lines.append(f'**{kind.capitalize()}**')
                    lines.extend(f'`{name}`' for name, cid in channel_entries.items())
            for kind, type_map in aliases.get('role_aliases', {}).items():
                role_entries = {name: data for name, data in type_map.items() if isinstance(data, dict) and data.get('channel_id') == channel_obj.id}
                if role_entries:
                    found_aliases = True
                    lines.append(f'**{kind.capitalize()}**')
                    for name, data in role_entries.items():
                        rid = data.get('role_id')
                        role = ctx.guild.get_role(rid)
                        mention = role.mention if role else f'<@&{rid}>'
                        lines.append(f'`{name}` → {mention}')
            if temp_room_obj:
                for alias_type, room_map in aliases.get('temp_room_aliases', {}).items():
                    for alias_name, data in room_map.items():
                        if temp_room_obj and data.get('room_name', '') == temp_room_obj.room_name:
                            lines.append(f'**{alias_type.capitalize()}`**')
                            lines.append(f'`{alias_name}`')
                            found_aliases = True
        if not found_aliases:
            return await self.send(interaction, content=f'\U0001F6AB No aliases found for the requested target: {target if target != None else channel_obj.mention}.')
        if target and target.lower() == 'all':
            pages = []
            chunk_size = 18
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title='All Aliases in Server',
                    description='\n'.join(chunk),
                    color=discord.Color.blue()
                )
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            await paginator.start()
        else:
            embed_title = f'Aliases for {channel_obj.mention}'
            embed = discord.Embed(
                title=embed_title,
                description='\n'.join(lines),
                color=discord.Color.blue()
            )
            await self.send(interaction, embed=embed)

    @commands.command(name='cmds', help='List command aliases routed to a specific channel, temp room, or all channels if "all" is provided.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_commands_text_command(
            self,
            ctx: commands.Context,
            target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or temp room name')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        aliases = self.bot.command_aliases.get(ctx.guild.id, self.bot.command_aliases.default_factory())
        lines = []
        channel_obj = await self.resolve_channel(ctx, target)
        temp_room_obj = None
        found_aliases = False
        for guild_temp_rooms in self.temp_rooms.values():
            for temp_channel in guild_temp_rooms.values():
                if temp_channel.room_name.lower() == channel_obj.name.lower():
                    temp_room_obj = temp_channel
                    break
            if temp_room_obj:
                break
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`cmds`) here.')
        if target and target.lower() == 'all':
            is_owner_or_dev, _ = await check_owner_dev_coord_mod(ctx, None)
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
                    lines.extend(f'`{name}`' for name in names)
                    found_aliases = True
            for kind, type_map in aliases.get('role_aliases', {}).items():
                grouped_by_channel = defaultdict(list)
                for name, data in type_map.items():
                    if isinstance(data, dict) and 'channel_id' in data:
                        grouped_by_channel[data['channel_id']].append((name, data.get('role_id')))
                for ch_id, entries in grouped_by_channel.items():
                    ch = ctx.guild.get_channel(ch_id)
                    ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                    lines.append(f'**{kind.capitalize()} Role Aliases**')
                    for name, rid in entries:
                        role = ctx.guild.get_role(rid)
                        mention = role.mention if role else f'<@&{rid}>'
                        lines.append(f'`{name}` → {mention}')
                    found_aliases = True
            for alias_type, room_map in aliases.get('temp_room_aliases', {}).items():
                for alias_name, data in room_map.items():
                    room_name = data.get('room_name') or '(unnamed room)'
                    lines.append(f'**{alias_type.capitalize()}**')
                    lines.append(f'`{alias_name}`')
                    found_aliases = True
        else:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            for kind, type_map in aliases.get('channel_aliases', {}).items():
                channel_entries = {name: cid for name, cid in type_map.items() if cid == channel_obj.id}
                if channel_entries:
                    found_aliases = True
                    lines.append(f'**{kind.capitalize()}**')
                    lines.extend(f'`{name}`' for name, cid in channel_entries.items())
            for kind, type_map in aliases.get('role_aliases', {}).items():
                role_entries = {name: data for name, data in type_map.items() if isinstance(data, dict) and data.get('channel_id') == channel_obj.id}
                if role_entries:
                    found_aliases = True
                    lines.append(f'**{kind.capitalize()}**')
                    for name, data in role_entries.items():
                        rid = data.get('role_id')
                        role = ctx.guild.get_role(rid)
                        mention = role.mention if role else f'<@&{rid}>'
                        lines.append(f'`{name}` → {mention}')
            if temp_room_obj:
                for alias_type, room_map in aliases.get('temp_room_aliases', {}).items():
                    for alias_name, data in room_map.items():
                        if temp_room_obj and data.get('room_name', '') == temp_room_obj.room_name:
                            lines.append(f'**{alias_type.capitalize()}**')
                            lines.append(f'`{alias_name}`')
                            found_aliases = True
        if not found_aliases:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB No aliases found for the requested target: {target if target != None else channel_obj.mention}.')
        if target and target.lower() == 'all':
            pages = []
            chunk_size = 18
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title='All Aliases in Server',
                    description='\n'.join(chunk),
                    color=discord.Color.blue()
                )
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            await paginator.start()
        else:
            embed_title = f'Aliases for {channel_obj.mention}'
            embed = discord.Embed(
                title=embed_title,
                description='\n'.join(lines),
                color=discord.Color.blue()
            )
            await self.handler.send_message(ctx, embed=embed)
 
    @app_commands.command(name='coords', description='Lists coordinators for a specific voice channel, all, or a member.')
    @app_commands.describe(target='"all", member or channel name/ID/mention')
    @is_owner_developer_coordinator_moderator_app_predicator(None)
    async def list_coordinators_app_command(
        self,
        interaction : discord.Interaction,
        target : Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content="This command must be used in a server.")
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in interaction.user.roles)
        member_obj = await self.resolve_member_app(interaction, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel_app(interaction, target)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        if member_obj:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall_app(interaction)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use (`coords`) in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if target and target.lower() == 'all':
                if not is_owner_or_dev:
                    return await self.send(interaction, content='\U0001F6AB You are not authorized to list all coordinators.')
                query = 'SELECT unnest(coordinator_channel_ids) AS channel_id, discord_snowflake FROM users WHERE coordinator_channel_ids IS NOT NULL'
                rows = await conn.fetch(query)
                if not rows:
                    return await self.send(interaction, content='\U0001F6AB No coordinators found in any voice channels.')
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
                if len(pages) == 1:
                    return await self.send(interaction, embed=pages[0], allowed_mentions=discord.AllowedMentions.none())
                paginator = UserPaginator(self.bot, interaction, pages)
                return await paginator.start()
            elif member_obj:
                query = 'SELECT coordinator_channel_ids, coordinator_room_names FROM users WHERE discord_snowflake=$1'
                row = await conn.fetchrow(query, member_obj.id)
                channels = []
                if row:
                    channels.extend(row.get('coordinator_channel_ids') or [])
                    if room_name and row.get('coordinator_room_names') and room_name in row['coordinator_room_names']:
                        channels.append(channel_obj.id if channel_obj else None)
                if not channels:
                    return await self.send(interaction, content=f'\U0001F6AB {member_obj.display_name} is not a coordinator in any channels.')
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
                if len(embeds) == 1:
                    return await self.send(interaction, embed=embeds[0], allowed_mentions=discord.AllowedMentions.none())
                paginator = UserPaginator(self.bot, interaction, embeds)
                return await paginator.start()
            elif channel_obj:
                if channel_obj.type != discord.ChannelType.voice:
                    return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
                query = 'SELECT discord_snowflake FROM users WHERE $1=ANY(coordinator_channel_ids)'
                rows = await conn.fetch(query, channel_obj.id)
                if room_name:
                    rows2 = await conn.fetch('SELECT discord_snowflake FROM users WHERE $1=ANY(coordinator_room_names)', room_name)
                    rows = list({r['discord_snowflake']: r for r in rows + rows2}.values())
                if not rows:
                    return await self.send(interaction, content=f'\U0001F6AB No coordinators found for {channel_obj.mention}.')
                lines = []
                for r in rows:
                    uid = r['discord_snowflake']
                    m = interaction.guild.get_member(uid)
                    if m:
                        lines.append(f'• {m.display_name} — <@{uid}>')
                if not lines:
                    return await self.send(interaction, content=f'\U0001F6AB No coordinators currently in {guild.name}.')
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
            
    @commands.command(name='coords', help='Lists coordinators for a specific voice channel, all, or a member.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Voice channel name, mention, ID, "all", or member ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in ctx.author.roles)
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target=None
        channel_obj = await self.resolve_channel(ctx, target)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        if member_obj:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall(ctx)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`coords`) in {channel_obj.mention}.')
        async with ctx.bot.db_pool.acquire() as conn:
            if target and target.lower()=='all':
                if not is_owner_or_dev:
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
                if len(pages) == 1:
                    return await self.handler.send_message(ctx, embed=pages[0], allowed_mentions=discord.AllowedMentions.none())
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            elif member_obj:
                query = 'SELECT coordinator_channel_ids, coordinator_room_names FROM users WHERE discord_snowflake=$1'
                row = await conn.fetchrow(query, member_obj.id)
                channels = []
                if row:
                    channels.extend(row.get('coordinator_channel_ids') or [])
                    if room_name and row.get('coordinator_room_names') and room_name in row['coordinator_room_names']:
                        channels.append(channel_obj.id if channel_obj else None)
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
                if len(embeds) == 1:
                    return await self.handler.send_message(ctx, embed=embeds[0], allowed_mentions=discord.AllowedMentions.none())
                paginator = Paginator(self.bot, ctx, embeds)
                return await paginator.start()
            elif channel_obj:
                if channel_obj.type != discord.ChannelType.voice:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
                query = 'SELECT discord_snowflake FROM users WHERE $1=ANY(coordinator_channel_ids)'
                rows = await conn.fetch(query, channel_obj.id)
                if room_name:
                    rows2 = await conn.fetch('SELECT discord_snowflake FROM users WHERE $1=ANY(coordinator_room_names)',room_name)
                    rows = list({r['discord_snowflake']: r for r in rows+rows2}.values())
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
    
    @app_commands.command(name='cstage', description='Create a stage in the current or specified channel.')
    @app_commands.describe(channel='Tag a voice/stage channel', duration='Duration of the stage (e.g., 1h, 30m)')
    @is_owner_developer_app_predicator()
    async def stage_create_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None,
        duration: str = '1'
    ):
        if not interaction.guild:
            return await self.send(interaction, content="This command must be used in a server.")
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        expires_at, duration_display = self.parse_duration(duration)
        skipped, muted, failed = [], [], []
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (guild_id, channel_id, initiator_id, expires_at, room_name)
                VALUES ($1, $2, $3, $4, $5)
                ON CONFLICT (guild_id, channel_id, room_name)
                DO UPDATE SET expires_at=EXCLUDED.expires_at, initiator_id=EXCLUDED.initiator_id
            ''', interaction.guild.id, channel_obj.id, interaction.user.id, expires_at, room_name)
            await conn.execute('''
                INSERT INTO stage_coordinators (guild_id, channel_id, discord_snowflake, room_name)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT DO NOTHING
            ''', interaction.guild.id, channel_obj.id, interaction.user.id, room_name)
            for user in channel_obj.members:
                if await is_owner_member(user, self.bot) or await is_developer_member(user, self.bot) \
                   or await is_coordinator_via_objects(user, channel_obj) or await is_moderator_via_objects(user, channel_obj) \
                   or user.id == interaction.user.id:
                    skipped.append(user)
                    continue
                try:
                    await conn.execute('''
                        INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, reason, target, room_name)
                        VALUES ($1,$2,$3,$4,$5,'room',$6)
                        ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
                        DO UPDATE SET expires_at=EXCLUDED.expires_at
                    ''', interaction.guild.id, user.id, channel_obj.id, expires_at, 'Stage active', room_name)
                    if user.voice and user.voice.channel.id == channel_obj.id:
                        await user.edit(mute=True)
                    muted.append(user)
                except Exception as e:
                    logger.warning(f'Failed to mute {user.name}: {e}')
                    failed.append(user)
        msg = f'{self.get_random_emoji()} Stage created for {duration_display} in {channel_obj.mention}.\nMuted {len(muted)} user(s). Skipped {len(skipped)}.'
        if failed:
            msg += f'\U0001F6AB Failed: {len(failed)}.'
        await self.send(interaction, content=msg, allowed_mentions=discord.AllowedMentions.none())
    
    @commands.command(name='cstage', help='Create a stage in the current or specified channel.')
    @is_owner_developer_predicator()
    async def stage_create_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        *,
        duration: str = '1'
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        expires_at, duration_display = self.parse_duration(duration)
        skipped, muted, failed = [], [], []
        room_name = ''
        temp_room_obj = None
        for guild_temp_rooms in self.temp_rooms.values():
            for temp_channel in guild_temp_rooms.values():
                if temp_channel.channel.id == channel_obj.id:
                    temp_room_obj = temp_channel
                    room_name = temp_channel.room_name
                    break
            if temp_room_obj:
                break
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (guild_id, channel_id, initiator_id, expires_at, room_name)
                VALUES ($1, $2, $3, $4, $5)
                ON CONFLICT (guild_id, channel_id, room_name)
                DO UPDATE SET expires_at=EXCLUDED.expires_at, initiator_id=EXCLUDED.initiator_id
            ''', ctx.guild.id, channel_obj.id, ctx.author.id, expires_at, room_name)
            await conn.execute('''
                INSERT INTO stage_coordinators (guild_id, channel_id, discord_snowflake, room_name)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT DO NOTHING
            ''', ctx.guild.id, channel_obj.id, ctx.author.id, room_name)
            for user in channel_obj.members:
                if await is_owner_member(user, ctx.bot) or await is_developer_member(user, ctx.bot) \
                   or await is_coordinator_via_objects(user, channel_obj) or await is_moderator_via_objects(user, channel_obj) \
                   or user.id == ctx.author.id:
                    skipped.append(user)
                    continue
                try:
                    await conn.execute('''
                        INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, reason, target, room_name)
                        VALUES ($1,$2,$3,$4,$5,'room',$6)
                        ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
                        DO UPDATE SET expires_at=EXCLUDED.expires_at
                    ''', ctx.guild.id, user.id, channel_obj.id, expires_at, 'Stage active', room_name)
                    if user.voice and user.voice.channel.id == channel_obj.id:
                        await user.edit(mute=True)
                    muted.append(user)
                except Exception as e:
                    logger.warning(f'Failed to mute {user.name}: {e}')
                    failed.append(user)
        msg = f'{self.get_random_emoji()} Stage created for {duration_display} in {channel_obj.mention}.\nMuted {len(muted)} user(s). Skipped {len(skipped)}.'
        if failed:
            msg += f'\n\U000026A0\U0000FE0F Failed: {len(failed)}.'
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())
        
    @app_commands.command(name='devs', description='Lists developers.')
    @app_commands.describe(target='"all" or a member snowflake ID/mention')
    @is_owner_developer_app_predicator()
    async def list_developers_app_command(
        self,
        interaction : discord.Interaction,
        target : Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content="This command must be used in a server.")
        async with self.bot.db_pool.acquire() as conn:
            pages = []
            member_obj = await self.resolve_member_app(interaction, target)
            if target is None:
                rows = await conn.fetch(
                    'SELECT discord_snowflake FROM users WHERE $1 = ANY(developer_guild_ids)',
                    interaction.guild.id
                )
                if not rows:
                    return await self.send(interaction, content=f'\U0001F6AB No developers are configured in {interaction.guild.name}.')
                for row in rows:
                    user = interaction.guild.get_member(row['discord_snowflake'])
                    name = user.display_name if user else f'User ID {row["discord_snowflake"]}'
                    embed = discord.Embed(
                        title = f'Developer: {name}',
                        color = discord.Color.blurple()
                    )
                    pages.append(embed)
            elif target.lower() == 'all':
                rows = await conn.fetch(
                    'SELECT discord_snowflake, developer_guild_ids FROM users WHERE array_length(developer_guild_ids, 1) > 0'
                )
                if not rows:
                    return await self.send(interaction, content='\U0001F6AB No developers are configured.')
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
                    return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from input: {target}.')
                row = await conn.fetchrow(
                    'SELECT developer_guild_ids FROM users WHERE discord_snowflake=$1',
                    member_obj.id
                )
                if not row or not row['developer_guild_ids']:
                    return await self.send(interaction, content=f'\U0001F6AB {member_obj.mention} is not a developer in any guilds.', allowed_mentions=discord.AllowedMentions.none())
                guilds = [self.bot.get_guild(gid).name for gid in row['developer_guild_ids'] if self.bot.get_guild(gid)]
                embed = discord.Embed(
                    title = f'Developer guilds for {member_obj.display_name}',
                    description = ', '.join(guilds) if guilds else 'No known guilds',
                    color = discord.Color.blurple()
                )
                pages.append(embed)
        paginator = UserPaginator(self.bot, interaction, pages)
        return await paginator.start()
        
    @commands.command(name='devs', help='Lists developers.')
    @is_owner_developer_predicator()
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", or user mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        async with self.bot.db_pool.acquire() as conn:
            pages = []
            member_obj = await self.resolve_member(ctx, target)
            if member_obj and member_obj.bot:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot list in the guilds the bot is a developer.')
            elif target and target.lower() == 'all':
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

    @app_commands.command(name='flags', description='List flag statistics.')
    @is_owner_developer_coordinator_moderator_app_predicator(None)
    async def list_flags_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in interaction.user.roles)
        member_obj = await self.resolve_member_app(interaction, target)
        if member_obj and member_obj.bot:
            return await self.send(interaction, content='\U0001F6AB You cannot list flags on the bot.')
        channel_obj = await self.resolve_channel_app(interaction, target)
        if member_obj:
            target=None
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall_app(interaction)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use this command (`flags`) for {channel_obj.mention}.')
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.send(interaction, content='\U0001F6AB Only owners or developers can list flags across all channels.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('SELECT channel_id, discord_snowflake FROM active_flags WHERE guild_id=$1', interaction.guild.id)
            if not rows: return await self.send(interaction, content='\U0001F6AB No users are flagged in any voice channels.')
            channel_map = defaultdict(list)
            for row in rows: channel_map[row['channel_id']].append(row['discord_snowflake'])
            pages = []
            for ch_id, user_ids in sorted(channel_map.items()):
                ch = interaction.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(
                    title=f'\U0001F6A9 Flag records for {ch_name}',
                    color=discord.Color.yellow()
                )
                for uid in user_ids:
                    m = interaction.guild.get_member(uid)
                    mention = m.mention if m else f'<@{uid}>'
                    embed.add_field(name=f'{interaction.guild.name}', value=f'• {mention}', inline=False)
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
        if member_obj:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('SELECT channel_id, reason FROM active_flags WHERE guild_id=$1 AND discord_snowflake=$2', interaction.guild.id, member_obj.id)
            rows = [r for r in rows if interaction.guild.get_channel(r['channel_id'])]
            if not rows:
                return await self.send(interaction, content=f'\U0001F6AB {member_obj.mention} is not flagged in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
            lines = []
            for r in rows:
                ch = interaction.guild.get_channel(r['channel_id'])
                ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                reason = r['reason'] or "No reason given"
                lines.append(f'• {ch_name} — {reason}')
            embed=discord.Embed(
                title=f'\U0001F6A9 Flag records for {member_obj.display_name}',
                description='\n'.join(lines),
                color=discord.Color.orange()
            )
            return await self.send(interaction, embed=embed,allowed_mentions=discord.AllowedMentions.none())
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('SELECT discord_snowflake FROM active_flags WHERE guild_id=$1 AND channel_id=$2', interaction.guild.id, channel_obj.id)
            if not rows:
                return await self.send(interaction, content=f'\U0001F6AB No users are flagged for {channel_obj.mention}.')
            pages = []
            chunk_size = 18
            for i in range(0, len(rows), chunk_size):
                chunk = rows[i:i+chunk_size]
                formatted_lines = []
                for record in chunk:
                    uid = record['discord_snowflake']
                    m = interaction.guild.get_member(uid)
                    formatted_lines.append(f'• {m.display_name} — <@{uid}>')
                if formatted_lines:
                    embed=discord.Embed(
                        title=f'\U0001F6A9 Flag records for {channel_obj.name}',
                        color=discord.Color.red()
                    )
                    embed.add_field(name=f'{guild.name}', value='\n'.join(formatted_lines), inline=False)
                    pages.append(embed)
            if not pages:
                return await self.send(interaction, content=f'\U0001F6AB No flagged users currently in {interaction.guild.name}.')
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
        
    @commands.command(name='flags', help='List flag statistics.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in ctx.author.roles)
        member_obj = await self.resolve_member(ctx, target)
        if member_obj and member_obj.bot:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot list flags on the bot.')
        channel_obj = await self.resolve_channel(ctx, target)
        if member_obj:
            target = None
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall(ctx)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`flags`) for {channel_obj.mention}.')
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
                return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} No users are flagged in any voice channels.')
            channel_map = defaultdict(list)
            for row in rows:
                channel_map[row['channel_id']].append(row['discord_snowflake'])
            pages = []
            for ch_id, user_ids in sorted(channel_map.items()):
                ch = ctx.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(
                    title=f'\U0001F6A9 Flag records for {ch_name}',
                    color=discord.Color.yellow()
                )
                for uid in user_ids:
                    m = ctx.guild.get_member(uid)
                    mention = m.mention if m else f'<@{uid}>'
                    embed.add_field(name=f'{ctx.guild.name}', value=f'• {mention}', inline=False)
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        if member_obj:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT channel_id, reason
                    FROM active_flags
                    WHERE guild_id = $1 AND discord_snowflake = $2
                ''', ctx.guild.id, member_obj.id)
                rows = [r for r in rows if ctx.guild.get_channel(r['channel_id'])]
                if not rows:
                    return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} is not flagged in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in rows:
                    ch = ctx.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                    reason = r['reason'] or "No reason given"
                    lines.append(f'• {ch_name} — {reason}')
                embed = discord.Embed(
                    title=f'\U0001F6A9 Flag records for {member_obj.display_name}',
                    description='\n'.join(lines),
                    color=discord.Color.orange()
                )
                return await self.handler.send_message(ctx, embed=embed, allowed_mentions=discord.AllowedMentions.all())
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT discord_snowflake
                    FROM active_flags
                    WHERE guild_id = $1 AND channel_id = $2
                ''', ctx.guild.id, channel_obj.id)
            if not rows:
                return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} No users are flagged for {channel_obj.mention}.')
            pages = []
            chunk_size = 18
            for i in range(0, len(rows), chunk_size):
                chunk = rows[i:i+chunk_size]
                formatted_lines = []
                for record in chunk:
                    uid = record['discord_snowflake']
                    member_obj = ctx.guild.get_member(uid)
                    formatted_lines.append(f'• {member_obj.display_name} — <@{uid}>')
                if formatted_lines:
                    embed = discord.Embed(
                        title=f'\U0001F6A9 Flag records for {channel_obj.name}',
                        color=discord.Color.red()
                    )
                    embed.add_field(name=f'{ctx.guild.name}', value='\n'.join(formatted_lines), inline=False)
                    pages.append(embed)
            if not pages:
                return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} No flagged users currently in {ctx.guild.name}.')
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
     
    @app_commands.command(name='hero', description='Toggle the special feature for a member.')
    @is_owner_predicator()
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    async def vegan_hero_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        member_obj = await self.resolve_member(interaction, member)
        if not member_obj:
            return await self.send(interaction, content='\U0001F6AB Could not resolve the member.')
        if member_obj.bot:
            return await self.send(interaction, content='\U0001F6AB You cannot make the bot a superhero.')
        if not hasattr(self, "super"):
            self.super = {"state": False, "members": set()}
        self.super["state"] = not self.super["state"]
        if self.super["state"]:
            self.super["members"].add(member_obj.id)
        else:
            self.super["members"].discard(member_obj.id)
            for channel in interaction.guild.channels:
                await self.unrestrict(interaction.guild, member_obj)
        state = 'ON' if self.super["state"] else 'OFF'
        await self.send(interaction, content=f'{self.get_random_emoji()} Superhero mode {state}. All mutes, bans, flags, and text-mutes against you have been forgiven.')
           
    @commands.command(name='hero', hidden=True)
    @is_owner_predicator()
    async def vegan_hero_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(description='Tag a member or include their snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve the member.')
        if member_obj.bot:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot a superhero.')
        if not hasattr(self, "super"):
            self.super = {"state": False, "members": set()}
        self.super["state"] = not self.super["state"]
        if self.super["state"]:
            self.super["members"].add(member_obj.id)
        else:
            self.super["members"].discard(member_obj.id)
            for channel in ctx.guild.channels:
                await self.unrestrict(ctx.guild, member_obj)
        state = f'ON' if self.super["state"] else f'OFF'
        await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Superhero mode {state}. All mutes, bans, flags, and text-mutes against you have been forgiven.')
                
    @app_commands.command(name='logs', description='Lists log channels for this guild or a specified guild ID.')
    @is_owner_developer_app_predicator()
    async def list_logs_app_command(
        self,
        interaction: discord.Interaction
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        entries = self.log_channels.get(interaction.guild.id, [])
        if not entries:
            return await self.send(interaction, content=f'\U0001F6AB No log channels configured in {interaction.guild.name}.')
        embed = discord.Embed(
            title=f'{self.get_random_emoji()} Log Channels',
            color=discord.Color.blue()
        )
        embed.set_footer(text=f'Guild ID: {interaction.guild.id}')
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT * FROM log_channels WHERE guild_id=$1;', interaction.guild.id)
            for row in rows:
                channel_obj = self.bot.get_channel(row['channel_id'])
                mention = channel_obj.mention if channel_obj else f'`{row["channel_id"]}`'
                enabled = row.get('enabled',False)
                log_type = row.get('type') or 'general'
                snowflakes = row.get('snowflakes') or []
                if log_type == 'general':
                    detail = f"Logs all events in {interaction.guild.name}"
                elif log_type == 'channel':
                    detail=f"Logs only events from channels: {', '.join(f'<#{s}>' for s in snowflakes)}"
                elif log_type == 'member':
                    detail=f"Logs only events for members: {', '.join(f'<@{s}>' for s in snowflakes)}"
                else:
                    detail="Unknown filter"
                embed.add_field(name=f"{mention} {'✅' if enabled else '\U0001F6AB'}", value=f"Type: **{log_type}**\n{detail}", inline=False)
        await self.send(interaction, embed=embed)
       
    @commands.command(name='logs', help='Lists log channels for this guild or a specified guild ID.')
    @is_owner_developer_predicator()
    async def list_logs_text_command(
        self,
        ctx: commands.Context
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        entries = self.log_channels.get(ctx.guild.id, [])
        if not entries:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB No log channels configured in {ctx.guild.name}.')
        embed = discord.Embed(
            title=f'{self.get_random_emoji()} Log Channels',
            color=discord.Color.blue()
        )
        embed.set_footer(text=f'Guild ID: {ctx.guild.id}')
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT * FROM log_channels WHERE guild_id=$1;', ctx.guild.id)
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
                    name=f"{mention} {'✅' if enabled else '\U0001F6AB'}",
                    value=f"Type: **{log_type}**\n{detail}",
                    inline=False
                )
        await self.handler.send_message(ctx, embed=embed)

    @app_commands.command(name='ls', description='List users cowed as going vegan in this guild.')
    @app_commands.describe(target='A member or channel snowflake ID/mention')
    @is_owner_developer_coordinator_moderator_app_predicator(None)
    async def list_members_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in interaction.user.roles)
        member_obj = await self.resolve_member_app(interaction, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel_app(interaction, target)
        if member_obj:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall(ctx)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use this command (`ls`) in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if member_obj:
                rows = await conn.fetch('''SELECT channel_id, created_at FROM active_cows WHERE guild_id=$1 AND discord_snowflake=$2''', interaction.guild.id, member_obj.id)
                if not rows:
                    return await self.send(interaction, content=f'\U0001F6AB {member_obj.mention} is not cowed in any channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in rows:
                    ch = interaction.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                    created_at = discord.utils.format_dt(r['created_at'],style='R') if r['created_at'] else ''
                    lines.append(f'• {ch_name} — {created_at}')
                embed=discord.Embed(
                    title=f'🐮 {member_obj.display_name}',
                    description='\n'.join(lines),
                    color=discord.Color.green()
                )
                return await self.send(interaction, embed=embed, allowed_mentions=discord.AllowedMentions.all())
            elif channel_obj:
                if channel_obj.type != discord.ChannelType.voice:
                    return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
                rows = await conn.fetch('''SELECT discord_snowflake, created_at FROM active_cows WHERE guild_id=$1 AND channel_id=$2''', interaction.guild.id, channel_obj.id)
                if not rows:
                    return await self.send(interaction, content=f'\U0001F6AB No users are cowed in {channel_obj.mention}.')
                lines = []
                for row in rows:
                    uid = row['discord_snowflake']
                    m = interaction.guild.get_member(uid)
                    if not m:
                        continue
                    created_at = discord.utils.format_dt(row['created_at'],style='R') if row['created_at'] else ''
                    lines.append(f'• {m.display_name} — <@{uid}> — {created_at}')
                if not lines:
                    return await self.send(interaction, content=f'\U0001F6AB No new vegans currently in {guild.name}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(lines), chunk_size):
                    chunk=lines[i:i+chunk_size]
                    embed=discord.Embed(
                        title=f'🐮 Vegan records for {channel_obj.name}',
                        description='\n'.join(chunk),
                        color=discord.Color.green()
                    )
                    pages.append(embed)
                paginator = UserPaginator(self.bot, interaction, pages)
                return await paginator.start()
                
    @commands.command(name='ls', help='List users cowed as going vegan in this guild.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_members_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in ctx.author.roles)
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        if member_obj and member_obj.bot:
            return await self.send(interaction, content='\U0001F6AB You cannot get a vegan status for the bot.')
        channel_obj = await self.resolve_channel(ctx, target)
        if member_obj:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall(ctx)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`ls`) in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if member_obj:
                rows = await conn.fetch('''
                    SELECT channel_id, created_at
                    FROM active_cows
                    WHERE guild_id = $1 AND discord_snowflake = $2
                ''', ctx.guild.id, member_obj.id)
                if not rows:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not cowed in any channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in rows:
                    ch = ctx.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                    created_at = discord.utils.format_dt(r['created_at'], style='R') if r['created_at'] else ''
                    lines.append(f'• {ch_name} — {created_at}')
                embed = discord.Embed(
                    title=f'🐮 {member_obj.display_name}',
                    description='\n'.join(lines),
                    color=discord.Color.green()
                )
                return await self.handler.send_message(ctx, embed=embed, allowed_mentions=discord.AllowedMentions.all())
            elif channel_obj:
                if channel_obj.type != discord.ChannelType.voice:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
                rows = await conn.fetch('''
                    SELECT discord_snowflake, created_at
                    FROM active_cows
                    WHERE guild_id = $1 AND channel_id = $2
                ''', ctx.guild.id, channel_obj.id)
                if not rows:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are cowed in {channel_obj.mention}.')
                lines = []
                for row in rows:
                    uid = row['discord_snowflake']
                    m = ctx.guild.get_member(uid)
                    if not m:
                        continue
                    created_at = discord.utils.format_dt(row['created_at'], style='R') if row['created_at'] else ''
                    lines.append(f'• {m.display_name} — <@{uid}> — {created_at}')
                if not lines:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No new vegans currently in {guild.name}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i + chunk_size]
                    embed = discord.Embed(
                        title=f'🐮 Vegan records for {channel_obj.name}',
                        description='\n'.join(chunk),
                        color=discord.Color.green()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
    
    @app_commands.command(name='mods', description='Lists moderator statistics.')
    @app_commands.describe(target='A member or channel snowflake ID/mention')
    @is_owner_developer_coordinator_moderator_app_predicator(None)
    async def list_moderators_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in interaction.user.roles)
        member_obj = await self.resolve_member_app(interaction, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel_app(interaction, target)
        if member_obj:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall_app(interaction)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use this command (`mods`) in {channel_obj.mention}.')
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.send(interaction, content='\U0001F6AB You are not authorized to list all moderators.')
            query = '''SELECT unnest(moderator_channel_ids) AS channel_id, discord_snowflake FROM users WHERE moderator_channel_ids IS NOT NULL'''
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(query)
            if not rows:
                return await self.send(interaction, content='\U0001F6AB No moderators found in any voice channels.')
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
            if member_obj.id != interaction.user.id and not is_mod_or_coord and not is_owner_or_dev:
                return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use this command (`mods`) for {member_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            query = '''SELECT moderator_channel_ids FROM users WHERE discord_snowflake=$1'''
            async with self.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(query, member_obj.id)
            if not row or not row['moderator_channel_ids']:
                return await self.send(interaction, content=f'\U0001F6AB {member_obj.display_name} is not a moderator in any channels.')
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
                return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
            query = '''SELECT discord_snowflake FROM users WHERE $1=ANY(moderator_channel_ids)'''
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(query, channel_obj.id)
            if not rows:
                return await self.send(interaction, content=f'\U0001F6AB No moderators found for {channel_obj.mention}.')
            lines = []
            for row in rows:
                uid = row['discord_snowflake']
                m = interaction.guild.get_member(uid)
                if not m:
                    continue
                lines.append(f'• {m.display_name} — <@{uid}>')
            if not lines:
                return await self.send(interaction, content=f'\U0001F6AB No moderators currently in {interaction.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(title=f'\U0001F6E1 Moderators for {channel_obj.name}', description='\n'.join(chunk), color=discord.Color.magenta())
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
        
    @commands.command(name='mods',help='Lists moderator statistics.')
    @is_owner_developer_coordinator_moderator_predicator(None)
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Voice channel name/mention/ID, "all", or member mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in ctx.author.roles)
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        if member_obj:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall(ctx)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`mods`) in {channel_obj.mention}.')
        async with ctx.bot.db_pool.acquire() as conn:
            if target and target.lower()=='all':
                if not is_owner_or_dev:
                    return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all moderators.')
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
                if len(pages)==1: return await self.handler.send_message(ctx, embed=pages[0],allowed_mentions=discord.AllowedMentions.none())
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            elif member_obj:
                query = 'SELECT moderator_channel_ids, moderator_room_names FROM users WHERE discord_snowflake=$1'
                row = await conn.fetchrow(query,member_obj.id)
                channels = []
                if row:
                    channels.extend(row.get('moderator_channel_ids') or [])
                    if row.get('moderator_room_names') and room_name in row['moderator_room_names']:
                        channels.append(channel_obj.id if channel_obj else None)
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
                if len(embeds) == 1:
                    return await self.handler.send_message(ctx, embed=embeds[0], allowed_mentions=discord.AllowedMentions.none())
                paginator = Paginator(self.bot, ctx, embeds)
                return await paginator.start()
            elif channel_obj:
                if channel_obj.type != discord.ChannelType.voice:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
                query = 'SELECT discord_snowflake FROM users WHERE $1=ANY(moderator_channel_ids)'
                rows = await conn.fetch(query,channel_obj.id)
                if room_name:
                    rows2 = await conn.fetch('SELECT discord_snowflake FROM users WHERE $1=ANY(moderator_room_names)',room_name)
                    rows = list({r['discord_snowflake']:r for r in rows + rows2}.values())
                if not rows: return await self.handler.send_message(ctx, content=f'\U0001F6AB No moderators found for {channel_obj.mention}.')
                lines = []
                for r in rows:
                    uid = r['discord_snowflake']
                    m = ctx.guild.get_member(uid)
                    if not m: continue
                    lines.append(f'• {m.display_name} — <@{uid}>')
                if not lines: return await self.handler.send_message(ctx, content=f'\U0001F6AB No moderators currently in {ctx.guild.name}.')
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
                paginator = Paginator(self.bot,ctx,pages)
                return await paginator.start()

    @app_commands.command(name='mutes', description='Lists mute statistics.')
    @is_owner_developer_coordinator_moderator_role_app_predicator(None)
    async def list_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in interaction.user.roles)
        member_obj = await self.resolve_member_app(interaction, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel_app(interaction, target)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        if member_obj:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall_app(interaction)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use this command (`mutes`) in {channel_obj.mention}.')
        if target and target.lower() == 'all' and is_owner_or_dev:
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''SELECT discord_snowflake, channel_id, room_name, expires_at, COALESCE(reason,'No reason provided') AS reason FROM active_voice_mutes WHERE guild_id=$1 AND target='user' ORDER BY channel_id, room_name, discord_snowflake''', interaction.guild.id)
            if not records:
                return await self.send(interaction, content=f'\U0001F6AB No muted users currently in {interaction.guild.name}.')
            grouped = defaultdict(list)
            for r in records:
                grouped[r['channel_id']].append(r)
            pages = []
            for channel_id, user_entries in sorted(grouped.items()):
                channel = interaction.guild.get_channel(channel_id)
                channel_name = channel.mention if channel else f'Unknown Channel ({channel_id})'
                chunk_size = 18
                for i in range(0, len(user_entries), chunk_size):
                    embed = discord.Embed(title=f'\U0001F507 Mutes records for {channel_name}', color=discord.Color.orange())
                    for r in user_entries[i:i + chunk_size]:
                        user_id = r['discord_snowflake']
                        member = interaction.guild.get_member(user_id)
                        name = member.display_name if member else f'User ID {user_id}'
                        mention = member.mention if member else f'`{user_id}`'
                        duration_str = self.fmt_duration(r['expires_at'])
                        embed.add_field(name=name, value=f'{mention}\nReason: {r["reason"]}\nDuration: {duration_str}', inline=False)
                    pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
        if member_obj:
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''SELECT guild_id, channel_id, expires_at, reason FROM active_voice_mutes WHERE discord_snowflake=$1 AND guild_id=$2 AND target='user' AND room_name=$3''', member_obj.id, interaction.guild.id, room_name)
            records = [r for r in records if interaction.guild.get_channel(r['channel_id'])]
            if not records:
                return await self.send(interaction, content=f'\U0001F6AB {member_obj.mention} is not muted in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
            lines = []
            for r in records:
                ch = interaction.guild.get_channel(r['channel_id'])
                ch_mention = ch.mention if ch else f'`{r["channel_id"]}`'
                lines.append(f'• {ch_mention} — {r["reason"]} — {self.fmt_duration(r["expires_at"])}')
            embed = discord.Embed(
                title=f'\U0001F507 Mute records for {member_obj.display_name}',
                description='\n'.join(lines),
                color=discord.Color.orange()
            )
            return await self.send(interaction, embed=embed, allowed_mentions=discord.AllowedMentions.none())
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''SELECT guild_id, channel_id, expires_at, reason, discord_snowflake FROM active_voice_mutes WHERE channel_id=$1 AND guild_id=$2 AND target='user' AND room_name=$3''', channel_obj.id, interaction.guild.id, room_name)
            if not records:
                return await self.send(interaction, content=f'\U0001F6AB No users are currently muted in {channel_obj.mention}.')
            lines = []
            for r in records:
                uid = r['discord_snowflake']
                member_obj = interaction.guild.get_member(uid)
                if not member_obj:
                    continue
                lines.append(f'• {member_obj.display_name} — <@{uid}> — {self.fmt_duration(r["expires_at"])}')
            if not lines:
                return await self.send(interaction, content=f'\U0001F6AB No muted users currently in {guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'\U0001F507 Mute records for {channel_obj.name}',
                    color=discord.Color.orange()
                )
                embed.add_field(name=f'{interaction.guild.name}', value='\n'.join(chunk), inline=False)
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()

    @commands.command(name='mutes', help='Lists mute statistics.')
    @is_owner_developer_coordinator_moderator_role_predicator(None)
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in ctx.author.roles)
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        if member_obj:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall(ctx)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`mutes`) in {channel_obj.mention}.')
        if target and target.lower() == 'all':
            if is_owner_or_dev:
                async with self.bot.db_pool.acquire() as conn:
                    records = await conn.fetch('''
                        SELECT discord_snowflake, channel_id, room_name, expires_at, COALESCE(reason, 'No reason provided') AS reason
                        FROM active_voice_mutes
                        WHERE guild_id = $1
                          AND target = 'user'
                        ORDER BY channel_id, room_name, discord_snowflake
                    ''', ctx.guild.id)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No muted users currently in {ctx.guild.name}.')
                grouped = defaultdict(list)
                for record in records:
                    grouped[record['channel_id']].append(record)
                pages = []
                for channel_id, user_entries in sorted(grouped.items()):
                    channel = ctx.guild.get_channel(channel_id)
                    channel_name = channel.mention if channel else f'Unknown Channel ({channel_id})'
                    chunk_size = 18
                    for i in range(0, len(user_entries), chunk_size):
                        embed = discord.Embed(
                            title=f'\U0001F507 Mutes records for {channel_name}',
                            color=discord.Color.orange()
                        )
                        for record in user_entries[i:i + chunk_size]:
                            user_id = record['discord_snowflake']
                            member = ctx.guild.get_member(user_id)
                            name = member.display_name if member else f'User ID {user_id}'
                            mention = member.mention if member else f'`{user_id}`'
                            reason = record['reason'] or 'No reason provided'
                            duration_str = self.fmt_duration(record['expires_at'])
                            embed.add_field(name=name, value=f'{mention}\nReason: {reason}\nDuration: {duration_str}', inline=False)
                        pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
        if member_obj:
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''
                    SELECT guild_id, channel_id, expires_at, reason
                    FROM active_voice_mutes
                    WHERE discord_snowflake = $1
                      AND guild_id = $2
                      AND target = 'user'
                      AND room_name = $3
                ''', member_obj.id, ctx.guild.id, room_name)
            records = [r for r in records if ctx.guild.get_channel(r['channel_id'])]
            if not records:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not muted in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
            description_lines = []
            for record in records:
                channel_obj = ctx.guild.get_channel(record['channel_id'])
                channel_mention = channel_obj.mention if channel_obj else f'`{record['channel_id']}`'
                reason = record['reason']
                duration_str = self.fmt_duration(record['expires_at'])
                description_lines.append(f'• {channel_mention} — {reason} — {duration_str}')
            embed = discord.Embed(title=f'\U0001F507 Mute records for {member_obj.display_name}', description='\n'.join(description_lines), color=discord.Color.orange())
            return await self.handler.send_message(ctx, embed=embed, allowed_mentions=discord.AllowedMentions.none())
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''
                    SELECT guild_id, channel_id, expires_at, reason, discord_snowflake
                    FROM active_voice_mutes
                    WHERE channel_id = $1
                      AND guild_id = $2
                      AND target = 'user'
                      AND room_name = $3
                ''', channel_obj.id, ctx.guild.id, room_name)
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
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No muted users currently in {channel_obj.mention}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(description_lines), chunk_size):
                    chunk = description_lines[i:i + chunk_size]
                    embed = discord.Embed(title=f'\U0001F507 Mute records for {channel_obj.name}', color=discord.Color.orange())
                    embed.add_field(name=f'{ctx.guild.name}', value='\n'.join(chunk), inline=False)
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()

    @app_commands.command(name='mlog', description='Create, modify, or delete a log channel.')
    @app_commands.describe(
        channel='Tag a channel or include its snowflake ID',
        action='create | modify | delete',
        log_type='Type of logs: member, channel, general',
        snowflakes='Optional list of member IDs to include in logs'
    )
    @is_owner_developer_app_predicator()
    async def modify_log_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None,
        action: Optional[str] = None,
        log_type: Optional[str] = None,
        snowflakes: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        if action not in ['create', 'modify', 'delete']:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid action `create`, `modify` or `delete`.')
        sf = [int(s) for s in snowflakes.split()] if snowflakes else []
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.text:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid text channel.')
        current_entries = self.log_channels.setdefault(interaction.guild.id, [])
        async with self.bot.db_pool.acquire() as conn:
            if action.lower() == 'delete':
                await conn.execute('DELETE FROM log_channels WHERE guild_id=$1 AND channel_id=$2;', interaction.guild.id, channel_obj.id)
                current_entries[:] = [e for e in current_entries if e['channel_id'] != channel_obj.id]
                return await self.send(interaction, content=f'{self.get_random_emoji()} Log channel {channel_obj.mention} deleted.')
            existing = await conn.fetchrow('SELECT * FROM log_channels WHERE guild_id=$1 AND channel_id=$2;', interaction.guild.id, channel_obj.id)
            if existing:
                await conn.execute('UPDATE log_channels SET type=$1, snowflakes=$2, enabled=TRUE WHERE guild_id=$3 AND channel_id=$4;', log_type, sf if sf else None, interaction.guild.id, channel_obj.id)
                for e in current_entries:
                    if e['channel_id'] == channel_obj.id:
                        e.update({'type': log_type, 'snowflakes': sf if sf else None, 'enabled': True})
                msg = f'{self.get_random_emoji()} Log channel {channel_obj.mention} updated with type `{log_type or "general"}`.'
            else:
                await conn.execute('INSERT INTO log_channels (guild_id, channel_id, type, snowflakes, enabled) VALUES ($1, $2, $3, $4, TRUE);', interaction.guild.id, channel_obj.id, log_type, sf if sf else None)
                current_entries.append({'guild_id': interaction.guild.id, 'channel_id': channel_obj.id, 'type': log_type, 'snowflakes': sf if sf else None, 'enabled': True})
                msg = f'{self.get_random_emoji()} Log channel {channel_obj.mention} created with type `{log_type or "general"}`.'
        self.log_channels[interaction.guild.id] = current_entries
        await self.send(interaction, content=msg, allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='mlog', help='Create, modify, or delete a log channel.')
    @is_owner_developer_predicator()
    async def modify_log_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        action: Optional[str] = commands.parameter(default=None, description='create | modify | delete'),
        log_type: Optional[str] = commands.parameter(default=None, description='Type of logs: member, channel, general'),
        *snowflakes: Optional[int]
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        if action not in ['create', 'modify', 'delete']:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid action `create`, `modify` or `delete`.')
        sf = [int(s) for s in snowflakes] if snowflakes else []
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.text:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid text channel.')
        current_entries = self.log_channels.setdefault(ctx.guild.id, [])
        async with self.bot.db_pool.acquire() as conn:
            if action.lower() == 'delete':
                await conn.execute('DELETE FROM log_channels WHERE guild_id=$1 AND channel_id=$2;', ctx.guild.id, channel_obj.id)
                current_entries[:] = [e for e in current_entries if e['channel_id'] != channel_obj.id]
                await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Log channel {channel_obj.mention} deleted.')
                return
            existing = await conn.fetchrow('SELECT * FROM log_channels WHERE guild_id=$1 AND channel_id=$2;', ctx.guild.id, channel_obj.id)
            if existing:
                await conn.execute('UPDATE log_channels SET type=$1, snowflakes=$2, enabled=TRUE WHERE guild_id=$3 AND channel_id=$4;', log_type, sf if sf else None, ctx.guild.id, channel_obj.id)
                for e in current_entries:
                    if e['channel_id'] == channel_obj.id:
                        e.update({'type': log_type, 'snowflakes': sf if sf else None, 'enabled': True})
                msg = f'{self.get_random_emoji()} Log channel {channel_obj.mention} updated with type `{log_type or "general"}`.'
            else:
                await conn.execute('INSERT INTO log_channels (guild_id, channel_id, type, snowflakes, enabled) VALUES ($1, $2, $3, $4, TRUE);', ctx.guild.id, channel_obj.id, log_type, sf if sf else None)
                current_entries.append({
                    'guild_id': ctx.guild.id,
                    'channel_id': channel_obj.id,
                    'type': log_type,
                    'snowflakes': sf if sf else None,
                    'enabled': True
                })
                msg = f'{self.get_random_emoji()} Log channel {channel_obj.mention} created with type `{log_type or "general"}`.'
        self.log_channels[ctx.guild.id] = current_entries
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())

    @app_commands.command(name='log', description='Toggle logging for a channel on or off.')
    @is_owner_developer_coordinator_app_predicator(None)
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    async def toggle_log_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.text:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid text channel.')
        target_role, allowed = await check_block_app(interaction, interaction.user, channel_obj)
        if not allowed:
            return await self.send(interaction, content=f'\U0001F6AB You must be {target_role} or higher to toggle logging for {channel_obj.mention}.')
        current_channels = self.log_channels.setdefault(interaction.guild.id, [])
        async with self.bot.db_pool.acquire() as conn:
            existing = await conn.fetchrow(
                'SELECT * FROM log_channels WHERE guild_id=$1 AND channel_id=$2;',
                interaction.guild.id, channel_obj.id
            )
            if existing:
                new_status = not existing['enabled']
                await conn.execute(
                    'UPDATE log_channels SET enabled=$1 WHERE guild_id=$2 AND channel_id=$3;',
                    new_status, interaction.guild.id, channel_obj.id
                )
                msg = f'{self.get_random_emoji()} Logging for {channel_obj.mention} toggled {"on" if new_status else "off"}.'
            else:
                msg = f'\U0001F6AB No logging for {channel_obj.mention} exists'
        self.log_channels[interaction.guild.id] = current_channels
        await self.send(interaction, content=msg, allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='log', help='Toggle logging for a channel on or off.')
    @is_owner_developer_coordinator_predicator(None)
    async def toggle_log_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.text:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid text channel.')
        target_role, allowed = await check_block(ctx, ctx.author, channel_obj)
        if not allowed:
            await self.handler.send_message(ctx, content=f'\U0001F6AB You must be {target_role} or higher to toggle logging for {channel_obj.mention}.')
            return
        current_channels = self.log_channels.setdefault(ctx.guild.id, [])
        async with self.bot.db_pool.acquire() as conn:
            existing = await conn.fetchrow(
                'SELECT * FROM log_channels WHERE guild_id=$1 AND channel_id=$2;',
                ctx.guild.id, channel_obj.id
            )
            if existing:
                new_status = not existing['enabled']
                await conn.execute(
                    'UPDATE log_channels SET enabled=$1 WHERE guild_id=$2 AND channel_id=$3;',
                    new_status, ctx.guild.id, channel_obj.id
                )
                msg = f'{self.get_random_emoji()} Logging for {channel_obj.mention} toggled {"on" if new_status else "off"}.'
            else:
                msg = f'\U0001F6AB No logging for {channel_obj.mention} exists'
        self.log_channels[ctx.guild.id] = current_channels
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())
        
    @app_commands.command(name='mstage', description='Mute/unmute a member in the active stage.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_developer_coordinator_moderator_app_predicator()
    async def stage_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        member_obj = await self.resolve_member_app(interaction, member)
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Invalid member for {member}.')
        room_name = ''
        is_temp_room = False
        channel_obj = await self.resolve_channel_app(interaction, None)
        async with self.bot.db_pool.acquire() as conn:
            records = await conn.fetch('''
                SELECT channel_id, room_name
                FROM active_stages
                WHERE guild_id=$1
            ''', interaction.guild.id)
            if not records:
                return await self.send(interaction, content='\U0001F6AB No active stage found.')
            temp_rooms = self.temp_rooms.get(interaction.guild.id, {})
            record = None
            for r in records:
                for t in temp_rooms.values():
                    if r['room_name'] == t.room_name:
                        record = r
                        room_name = t.room_name
                        is_temp_room = True
                        break
                if record:
                    break
            stage = await conn.fetchrow('''
                SELECT initiator_id
                FROM active_stages
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', interaction.guild.id, channel_obj.id, room_name)
            if not stage:
                return await self.send(interaction, content='\U0001F6AB No active stage found.')
            if member_obj.id == stage['initiator_id']:
                return await self.send(interaction, content='\U0001F6AB You cannot mute the stage initiator.')
            is_owner_or_dev, _ = await check_owner_dev_coord_app(interaction, channel_obj)
            is_coordinator = await conn.fetchval('''
                SELECT 1
                FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
            ''', interaction.guild.id, channel_obj.id, room_name, interaction.user.id)
            if not is_coordinator and not is_owner_or_dev:
                return await self.send(interaction, content='\U0001F6AB Only stage coordinators or above can use this command.')
            try:
                await member_obj.edit(mute=not member_obj.voice.mute)
                return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has been {"muted" if member_obj.voice.mute else "unmuted"}.')
            except Exception as e:
                logger.warning(f'Failed to toggle mute: {e}')
            return await self.send(interaction, content=f'\U0001F6AB Failed to toggle mute for {member_obj.mention}.')
                
    @commands.command(name='mstage', help='Mute/unmute a member in the active stage.')
    @is_owner_developer_coordinator_moderator_predicator()
    async def stage_mute_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from target: `{member}`.')
        channel_obj = await self.resolve_channel(ctx, None)
        room_name = ''
        is_temp_room = False
        async with self.bot.db_pool.acquire() as conn:
            records = await conn.fetch('''
                SELECT channel_id, room_name
                FROM active_stages
                WHERE guild_id=$1
            ''', ctx.guild.id)
            if not records:
                return await self.handler.send_message(ctx, content='\U0001F6AB No active stage found.')
            temp_rooms = self.temp_rooms.get(ctx.guild.id, {})
            record = None
            for r in records:
                for t in temp_rooms.values():
                    if r['room_name'] == t.room_name:
                        record = r
                        is_temp_room = True
                        room_name = t.room_name
                        break
                if record:
                    break
            stage = await conn.fetchrow('''
                SELECT initiator_id
                FROM active_stages
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', ctx.guild.id, channel_obj.id, room_name)
            if not stage:
                return await self.handler.send_message(ctx, content='\U0001F6AB No active stage found.')
            if member_obj.id == stage['initiator_id']:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot mute the stage initiator.')
            is_owner_or_dev, _ = await check_owner_dev_coord(ctx, channel_obj)
            is_coordinator = await conn.fetchval('''
                SELECT 1
                FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
            ''', ctx.guild.id, channel_obj.id, room_name, ctx.author.id)
            if not is_coordinator and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB Only stage coordinators can use this command.')
            try:
                await member_obj.edit(mute=not member_obj.voice.mute)
                return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been {"muted" if member_obj.voice.mute else "unmuted"}.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                logger.warning(f'Failed to toggle mute: {e}')
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Failed to toggle mute for {member_obj.display_name}.')
    
    @app_commands.command(name='pstage', description='Promote/demote a member as stage coordinator.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_developer_coordinator_moderator_app_predicator()
    async def stage_promote_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        member_obj = await self.resolve_member_app(interaction, member)
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from target: `{member}`.')
        channel_obj = await self.resolve_channel_app(interaction, None)
        temp_rooms = self.temp_rooms.get(interaction.guild.id, {})
        room_name = ''
        is_temp_room = False
        for temp_channel in temp_rooms.values():
            if temp_channel.room_name == channel_obj.name:
                channel_obj = temp_channel
                is_temp_room = True
                room_name = temp_channel.room_name
                break
        async with self.bot.db_pool.acquire() as conn:
            stage = await conn.fetchrow('SELECT initiator_id, channel_id, room_name FROM active_stages WHERE guild_id=$1', interaction.guild.id)
            if not stage:
                return await self.send(interaction, content='\U000026A0\U0000FE0F No active stage found.')
            is_owner_or_dev, is_coord_or_mod = await check_owner_dev_coord_mod_app(interaction, channel_obj)
            if not is_owner_or_dev and not is_coord_or_mod:
                return await self.send(interaction, content=f'\U0001F6AB You do not have permission to moderate this stage.')
            if member_obj.id == stage['initiator_id']:
                return await self.send(interaction, content='\U0001F6AB Cannot change the initiator role.')
            is_coordinator = await conn.fetchval('''
                SELECT 1 FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
            ''', interaction.guild.id, channel_obj.id, room_name, member_obj.id)
            try:
                if is_coordinator:
                    await conn.execute('''
                        DELETE FROM stage_coordinators
                        WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
                    ''', interaction.guild.id, channel_obj.id, room_name, member_obj.id)
                    if member_obj.voice and not member_obj.voice.mute:
                        await member_obj.edit(mute=True)
                    return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has been demoted from stage coordinator.', allowed_mentions=discord.AllowedMentions.none())
                else:
                    await conn.execute('''
                        INSERT INTO stage_coordinators (guild_id, channel_id, room_name, discord_snowflake)
                        VALUES ($1,$2,$3,$4) ON CONFLICT DO NOTHING
                    ''', interaction.guild.id, channel_obj.id, room_name, member_obj.id)
                    await conn.execute('''
                        DELETE FROM active_voice_mutes
                        WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target='room'
                    ''', interaction.guild.id, member_obj.id, channel_obj.id)
                    if member_obj.voice and member_obj.voice.mute:
                        await member_obj.edit(mute=False)
                    return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has been promoted to stage coordinator.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                logger.warning(f'Failed to toggle promotion: {e}')
            return await self.send(interaction, content=f'\U0001F6AB Failed to toggle promotion for {member_obj.mention}.')
    
    @commands.command(name='pstage', help='Promote/demote a member as stage coordinator.')
    @is_owner_developer_coordinator_moderator_predicator()
    async def stage_promote_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from target: `{member}`.')
        channel_obj = await self.resolve_channel(ctx, None)
        temp_rooms = self.temp_rooms.get(ctx.guild.id, {})
        room_name = ''
        is_temp_room = False
        for temp_channel in temp_rooms.values():
            if temp_channel.room_name == channel_obj.name:
                channel_obj = temp_channel
                is_temp_room = True
                room_name = temp_channel.room_name
                break
        async with self.bot.db_pool.acquire() as conn:
            stage = await conn.fetchrow('SELECT initiator_id, channel_id, room_name FROM active_stages WHERE guild_id=$1', ctx.guild.id)
            if not stage:
                return await self.handler.send_message(ctx, content='\U0001F6AB No active stage found.')
            is_owner_or_dev, is_coord_or_mod = await check_owner_dev_coord_mod(ctx, channel_obj)
            if not is_owner_or_dev and not is_coord_or_mod:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to moderate this stage.')
            if member_obj.id == stage['initiator_id']:
                return await self.handler.send_message(ctx, content='\U0001F6AB Cannot change the initiator role.')
            is_coordinator = await conn.fetchval('''
                SELECT 1 FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
            ''', ctx.guild.id, channel_obj.id, room_name, member_obj.id)
            try:
                if is_coordinator:
                    await conn.execute('''
                        DELETE FROM stage_coordinators
                        WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
                    ''', ctx.guild.id, channel_obj.id, room_name, member_obj.id)
                    if member_obj.voice and not member_obj.voice.mute:
                        await member_obj.edit(mute=True)
                    return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been demoted from stage coordinator.', allowed_mentions=discord.AllowedMentions.none())
                    
                else:
                    await conn.execute('''
                        INSERT INTO stage_coordinators (guild_id, channel_id, room_name, discord_snowflake)
                        VALUES ($1,$2,$3,$4) ON CONFLICT DO NOTHING
                    ''', ctx.guild.id, channel_obj.id, room_name, member_obj.id)
                    await conn.execute('''
                        DELETE FROM active_voice_mutes
                        WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target='room'
                    ''', ctx.guild.id, member_obj.id, channel_obj.id)
                    if member_obj.voice and member_obj.voice.mute:
                        await member_obj.edit(mute=False)
                    return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been promoted to stage coordinator.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                logger.warning(f'Failed to toggle promotion: {e}')
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Failed to toggle promotion for {member_obj.mention}.')

    @app_commands.command(name='migrate', description='Migrate a temporary room to a new channel.')
    @app_commands.describe(old_name='Old temporary room name', new_channel='New channel to migrate to')
    @is_owner_developer_coordinator_app_predicator()
    async def rename_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        old_name: str,
        new_channel: discord.abc.GuildChannel
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        async with self.bot.db_pool.acquire() as conn:
            rooms = await conn.fetch(
                'SELECT room_name, owner_snowflake, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_name=$2',
                interaction.guild.id, old_name
            )
            if not rooms:
                return await self.send(interaction, content=f'\U0001F6AB No temporary room named `{old_name}` found.')
            if len(rooms) > 1:
                return await self.send(interaction, content=f'\U0001F6AB Multiple temporary rooms named `{old_name}` exist. Migration failed.')
            temp = rooms[0]
            channel_obj = await self.resolve_channel_app(interaction, new_channel.id)
            if channel_obj.type != discord.ChannelType.voice:
                return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
            is_owner = temp['owner_snowflake'] == interaction.user.id
            is_owner_or_dev, _ = await check_owner_dev_coord_app(interaction, channel_obj)
            if not (is_owner_or_dev or is_owner):
                return await self.send(interaction, content='\U0001F6AB Only the owner or developers can migrate this room.')
            await conn.execute(
                'UPDATE temporary_rooms SET room_name=$3, room_snowflake=$4 WHERE guild_snowflake=$1 AND room_name=$2',
                interaction.guild.id, old_name, channel_obj.name, channel_obj.id
            )
            tables = ['active_bans','active_text_mutes','active_voice_mutes','active_stages','stage_coordinators','active_caps','command_aliases']
            for table in tables:
                await conn.execute(
                    f'UPDATE {table} SET room_name=$3, channel_id=$4 WHERE guild_id=$1 AND room_name=$2',
                    interaction.guild.id, old_name, channel_obj.name, channel_obj.id
                )
            await conn.execute(
                'UPDATE users SET coordinator_room_names=array_replace(coordinator_room_names, $1, $2) WHERE $1=ANY(coordinator_room_names)',
                old_name, channel_obj.name
            )
            await conn.execute(
                'UPDATE users SET coordinator_channel_ids=array_replace(coordinator_channel_ids, $1::bigint, $2::bigint) WHERE $1=ANY(coordinator_channel_ids)',
                temp['room_snowflake'], channel_obj.id
            )
            await conn.execute(
                'UPDATE users SET moderator_room_names=array_replace(moderator_room_names, $1, $2) WHERE $1=ANY(moderator_room_names)',
                old_name, channel_obj.name
            )
            await conn.execute(
                'UPDATE users SET moderator_channel_ids=array_replace(moderator_channel_ids, $1::bigint, $2::bigint) WHERE $1=ANY(moderator_channel_ids)',
                temp['room_snowflake'], channel_obj.id
            )
            if interaction.guild.id in self.bot.command_aliases:
                guild_alias_root = self.bot.command_aliases[interaction.guild.id]
                if 'temp_room_aliases' in guild_alias_root:
                    temp_room_aliases = guild_alias_root['temp_room_aliases']
                    for alias_type, aliases in temp_room_aliases.items():
                        for alias_name, alias_data in aliases.items():
                            if alias_data.get('room_name') == old_name:
                                alias_data['channel_id'] = channel_obj.id
                                alias_data['room_name'] = channel_obj.name
            if interaction.guild.id in self.temp_rooms:
                if old_name in self.temp_rooms[interaction.guild.id]:
                    temp_channel_obj = self.temp_rooms[interaction.guild.id].pop(old_name)
                    temp_channel_obj.room_name = channel_obj.name
                    temp_channel_obj.channel = channel_obj
                    self.temp_rooms[interaction.guild.id][channel_obj.name] = temp_channel_obj
            return await self.send(interaction, content=f'{self.get_random_emoji()} Temporary room `{old_name}` migrated to {channel_obj.mention} and renamed to `{channel_obj.name}`.')

    @commands.command(name='migrate', help='Migrate a temporary room to a new channel by snowflake.')
    @is_owner_developer_coordinator_predicator()
    async def migrate_temp_room_command(
        self,
        ctx,
        old_name: Optional[str] = commands.parameter(default=None, description='Provide a channel name'),
        new_target: Optional[int] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        async with self.bot.db_pool.acquire() as conn:
            rooms = await conn.fetch(
                'SELECT room_name, owner_snowflake, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_name=$2',
                ctx.guild.id, old_name
            )
            if not rooms:
                return await self.handler.send_message(ctx, content=f'No temporary room named `{old_name}` found.')
            if len(rooms) > 1:
                return await self.handler.send_message(ctx, content=f'Multiple temporary rooms named `{old_name}` exist. Migration failed.')
            temp = rooms[0]
            channel_obj = await self.resolve_channel(ctx, new_target)
            if channel_obj.type != discord.ChannelType.voice:
                return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            is_owner = temp['owner_snowflake'] == ctx.author.id
            is_owner_or_dev, _ = await check_owner_dev_coord(ctx, channel_obj)
            if not (is_owner_or_dev or is_owner):
                return await self.handler.send_message(ctx, content='Only the owner or developers can migrate this room.')
            await conn.execute(
                'UPDATE temporary_rooms SET room_name=$3, room_snowflake=$4 WHERE guild_snowflake=$1 AND room_name=$2',
                ctx.guild.id, old_name, channel_obj.name, channel_obj.id
            )
            tables = ['active_bans','active_text_mutes','active_voice_mutes','active_stages','stage_coordinators','active_caps','command_aliases']
            for table in tables:
                await conn.execute(
                    f'UPDATE {table} SET room_name=$3, channel_id=$4 WHERE guild_id=$1 AND room_name=$2',
                    ctx.guild.id, old_name, channel_obj.name, channel_obj.id
                )
            await conn.execute(
                'UPDATE users SET coordinator_room_names=array_replace(coordinator_room_names, $1, $2) WHERE $1=ANY(coordinator_room_names)',
                old_name, channel_obj.name
            )
            await conn.execute(
                'UPDATE users SET coordinator_channel_ids=array_replace(coordinator_channel_ids, $1::bigint, $2::bigint) WHERE $1=ANY(coordinator_channel_ids)',
                temp['room_snowflake'], channel_obj.id
            )
            await conn.execute(
                'UPDATE users SET moderator_room_names=array_replace(moderator_room_names, $1, $2) WHERE $1=ANY(moderator_room_names)',
                old_name, channel_obj.name
            )
            await conn.execute(
                'UPDATE users SET moderator_channel_ids=array_replace(moderator_channel_ids, $1::bigint, $2::bigint) WHERE $1=ANY(moderator_channel_ids)',
                temp['room_snowflake'], channel_obj.id
            )
            if ctx.guild.id in self.bot.command_aliases:
                guild_alias_root = self.bot.command_aliases[ctx.guild.id]
                if 'temp_room_aliases' in guild_alias_root:
                    temp_room_aliases = guild_alias_root['temp_room_aliases']
                    for alias_type, aliases in temp_room_aliases.items():
                        for alias_name, alias_data in aliases.items():
                            if alias_data.get('room_name') == old_name:
                                alias_data['channel_id'] = channel_obj.id
                                alias_data['room_name'] = channel_obj.name
            if ctx.guild.id in self.temp_rooms:
                if old_name in self.temp_rooms[ctx.guild.id]:
                    temp_channel_obj = self.temp_rooms[ctx.guild.id].pop(old_name)
                    temp_channel_obj.room_name = channel_obj.name
                    temp_channel_obj.channel = channel_obj
                    self.temp_rooms[ctx.guild.id][channel_obj.name] = temp_channel_obj
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Temporary room `{old_name}` migrated to {channel_obj.mention} and renamed to `{channel_obj.name}`.')

    @app_commands.command(name='rmv', description='Move all the members in one room to another.')
    @app_commands.describe(source_id='Tag the source channel or include its snowflake ID',target_id='Tag the target channel or include its snowflake ID')
    @is_owner_app_predicator()
    async def room_move_all_app_command(
        self,
        interaction: discord.Interaction,
        source_id: Optional[str] = None,
        target_id: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        source_channel = interaction.guild.get_channel(int(source_id))
        target_channel = interaction.guild.get_channel(int(target_id))
        if not source_channel or not target_channel:
            return await self.send(interaction, content='\U0001F6AB One or both channels are invalid.')
        if not isinstance(source_channel, discord.VoiceChannel) or not isinstance(target_channel, discord.VoiceChannel):
            raise ValueError('\U0001F6AB Both source and target must be voice channels.')
        for member in source_channel.members:
            try:
                await member.move_to(target_channel)
            except discord.Forbidden:
                logger.warning(f'\U0001F6AB Missing permissions to move {member}.')
            except discord.HTTPException:
                logger.warning(f'\U0001F6AB Failed to move {member} due to a network error.')
        await self.send(interaction, content=f'{self.get_random_emoji()} Moved all members from {source_channel.mention} to {target_channel.mention}.')
        
    @commands.command(name='rmv', help='Move all the members in one room to another.')
    @is_owner_predicator()
    async def room_move_all_text_command(
        self,
        ctx: commands.Context,
        source_id: Optional[int] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        target_id: Optional[int] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        source_channel = ctx.guild.get_channel(source_id)
        target_channel = ctx.guild.get_channel(target_id)
        if not source_channel or not target_channel:
            await self.handler.send_message(ctx, content='\U0001F6AB One or both channel IDs are invalid.')
            return
        if not isinstance(source_channel, discord.VoiceChannel) or not isinstance(target_channel, discord.VoiceChannel):
            raise ValueError('\U0001F6AB Both source and target must be voice channels.')
        for member in source_channel.members:
            try:
                await member.move_to(target_channel)
            except discord.Forbidden:
                logger.warning(f'\U0001F6AB Missing permissions to move {member}.')
            except discord.HTTPException:
                logger.warning(f'\U0001F6AB Failed to move {member} due to a network error.')
        await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Moved all members from {source_channel.mention} to {target_channel.mention}.')
        
    @app_commands.command(name='roleid', description='Get the ID of a role by name in this server.')
    @is_owner_developer_app_predicator()
    @app_commands.describe(role_name='The name of the role to look up')
    async def get_role_id_app_command(
        self,
        interaction: discord.Interaction,
        role_name: str
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in a server.')
        role = discord.utils.get(interaction.guild.roles, name=role_name)
        if role:
            await self.send(interaction, content=f'{self.get_random_emoji()} Role `{role.name}` has ID `{role.id}`.')
        else:
            await self.send(interaction, content=f'\U0001F6AB No role named "{role_name}" found in this server.')

    @commands.command(name='roleid', help='Get the ID of a role by name in this server.')
    @is_owner_developer_predicator()
    async def get_role_id(self, ctx: commands.Context, *, role_name: str):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in a server.')
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Role `{role.name}` has ID `{role.id}`.')
        else:
            await self.handler.send_message(ctx, content=f'\U0001F6AB No role named "{role_name}" found in this server.')
    
    @app_commands.command(name='rmute', description='Mutes all members in a VC (except yourself).')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_app_predicator()
    async def room_mute_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in a server.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        muted_members, skipped_members, failed_members = [], [], []
        async with self.bot.db_pool.acquire() as conn:
            for member in channel_obj.members:
                if member.id == interaction.user.id:
                    continue
                try:
                    row = await conn.fetchrow('''
                        SELECT 1 FROM active_voice_mutes
                        WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target='user'
                        AND (expires_at IS NULL OR expires_at > NOW())
                    ''', interaction.guild.id, member.id, channel_obj.id)
                    if row:
                        skipped_members.append(member)
                        continue
                    if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                        await member.edit(mute=True)
                    await conn.execute('''
                        INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, target, room_name, expires_at, reason)
                        VALUES ($1, $2, $3, 'user', '', NULL, $4)
                    ''', interaction.guild.id, member.id, channel_obj.id, 'Muted via room_mute')
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'mute', member.id, interaction.user.id, interaction.guild.id, channel_obj.id, 'Muted via room_mute')
                    muted_members.append(member)
                except Exception as e:
                    logger.warning(f'Mute failed for {member.mention}: {e}')
                    failed_members.append(member)
        summary = f'{self.get_random_emoji()} Muted {len(muted_members)} member(s) in {channel_obj.mention}.'
        if skipped_members:
            summary += f'\n\U000026A0\U0000FE0F Skipped {len(skipped_members)}.'
        if failed_members:
            summary += f'\n\U0001F6AB Failed to mute {len(failed_members)}.'
        return await self.send(interaction, content=summary, allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='rmute', help='Mutes all members in a VC (except yourself).')
    @is_owner_predicator()
    async def room_mute_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        muted_members, skipped_members, failed_members = [], [], []
        async with self.bot.db_pool.acquire() as conn:
            for member in channel_obj.members:
                if member.id == ctx.author.id:
                    continue
                try:
                    row = await conn.fetchrow('''
                        SELECT 1 FROM active_voice_mutes
                        WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target='user'
                        AND (expires_at IS NULL OR expires_at > NOW())
                    ''', ctx.guild.id, member.id, channel_obj.id)
                    if row:
                        skipped_members.append(member)
                        continue
                    if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                        await member.edit(mute=True)
                    await conn.execute('''
                        INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, target, room_name, expires_at, reason)
                        VALUES ($1, $2, $3, 'user', '', NULL, $4)
                    ''', ctx.guild.id, member.id, channel_obj.id, 'Muted via room_mute')
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'mute', member.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Muted via room_mute')
                    muted_members.append(member)
                except Exception as e:
                    logger.warning(f'Mute failed for {member.mention}: {e}')
                    failed_members.append(member)
        summary = f'{self.get_random_emoji()} Muted {len(muted_members)} member(s) in {channel_obj.mention}.'
        if skipped_members:
            summary += f'\n\U000026A0\U0000FE0F Skipped {len(skipped_members)}.'
        if failed_members:
            summary += f'\n\U0001F6AB Failed to mute {len(failed_members)}.'
        return await self.handler.send_message(ctx, content=summary, allowed_mentions=discord.AllowedMentions.none())

    @app_commands.command(name='smute', description='Mutes a member throughout the entire guild.')
    @app_commands.describe(member='Tag a member or include their snowflake ID', reason='Optional reason (required for 7 days or more)')
    @is_server_muter_app_predicator()
    async def server_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None,
        reason: Optional[str] = ''
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        member_obj = await self.resolve_member_app(interaction, member)
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from target: `{member}`.')
        if member_obj.bot:
            return await self.send(interaction, content='\U0001F6AB You cannot server mute the bot.')
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
            ''', member_obj.id, interaction.guild.id)
            await conn.execute('''
                INSERT INTO active_server_voice_mutes (guild_id, discord_snowflake, reason)
                VALUES ($1, $2, $3)
                ON CONFLICT (guild_id, discord_snowflake) DO UPDATE
                SET reason = EXCLUDED.reason
            ''', interaction.guild.id, member_obj.id, reason or 'No reason provided')
        if member_obj.voice and member_obj.voice.channel:
            await member_obj.edit(mute=True)
        return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has been server muted for reason: {reason}', allowed_mentions=discord.AllowedMentions.none())
            
    @commands.command(name='smute', help='Mutes a member throughout the entire guild.')
    @is_server_muter_predicator()
    async def server_mute_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        *,
        reason: Optional[str] = commands.parameter(default='', description='Optional reason (required for 7 days or more)')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
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
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been server muted for reason: {reason}', allowed_mentions=discord.AllowedMentions.none())

    @app_commands.command(name='stages', description='Lists stage mute statistics.')
    @app_commands.describe(target='"all", channel name/ID/mention')
    @is_owner_developer_coordinator_app_predicator(None)
    async def list_stages_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in interaction.user.roles)
        channel_obj=await self.resolve_channel_app(interaction, target)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        is_owner_or_dev,is_mod_or_coord = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        async with self.bot.db_pool.acquire() as conn:
            if target.lower() == 'all':
                if not is_owner_or_dev: return await self.send(interaction, content='\U0001F6AB Only owners/devs can view all stages.')
                stages = await conn.fetch('''
                    SELECT s.channel_id, s.room_name, s.initiator_id, s.expires_at, COUNT(v.discord_snowflake) AS active_mutes
                    FROM active_stages s LEFT JOIN active_voice_mutes v ON s.guild_id=v.guild_id AND s.channel_id=v.channel_id AND v.target='room'
                    WHERE s.guild_id=$1 GROUP BY s.channel_id, s.room_name, s.initiator_id, s.expires_at ORDER BY s.channel_id
                ''', interaction.guild.id)
                if not stages:
                    return await self.send(interaction, content=f'\U0001F6AB No active stages in {interaction.guild.name}.')
                pages, chunk_size = [], 8
                for i in range(0, len(stages), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F399 Active Stages in {interaction.guild.name}',
                        color=discord.Color.purple()
                    )
                    for s in stages[i:i+chunk_size]:
                        ch = interaction.guild.get_channel(s['channel_id'])
                        ch_name = ch.mention if ch else f'Unknown Channel ({s["channel_id"]})'
                        embed.add_field(name=ch_name, value=f'Active stage mutes: {s["active_mutes"]}', inline=False)
                    pages.append(embed)
                paginator = UserPaginator(self.bot, interaction, pages)
                return await paginator.start()
            if not (is_owner_or_dev or is_mod_or_coord or is_team_member):
                return await self.send(interaction, content=f'\U0001F6AB You do not have permission to view stages in {channel_obj.mention}.')
            stage = await conn.fetchrow('''
                SELECT initiator_id, expires_at FROM active_stages WHERE guild_id=$1 AND channel_id=$2 AND room_name = $3
            ''', interaction.guild.id, channel_obj.id, room_name)
            if not stage:
                return await self.send(interaction, content=f'\U0001F6AB No active stage in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            mutes = await conn.fetch('''
                SELECT discord_snowflake, expires_at, reason FROM active_voice_mutes WHERE guild_id=$1 AND channel_id=$2 AND target='room' AND room_name = $3
            ''', interaction.guild.id, channel_obj.id, room_name)
            coordinators = await conn.fetch('''
                SELECT discord_snowflake FROM stage_coordinators WHERE guild_id=$1 AND channel_id=$2
            ''', interaction.guild.id, channel_obj.id)
            coordinator_mentions = [member.mention for c in coordinators if (member:=interaction.guild.get_member(c['discord_snowflake']))]
            coordinator_str = ', '.join(coordinator_mentions) if coordinator_mentions else 'No coordinators'
            initiator = interaction.guild.get_member(stage['initiator_id'])
            initiator_name = initiator.mention if initiator else f'`{stage["initiator_id"]}`'
            expires = self.fmt_duration(stage['expires_at']) if stage['expires_at'] else 'No expiration'
            lines = []
            for m in mutes:
                user = interaction.guild.get_member(m['discord_snowflake'])
                duration_str = self.fmt_duration(m['expires_at']) if m['expires_at'] else 'No expiration'
                reason = m['reason'] or 'No reason provided'
                lines.append(f'• {user.mention} — {reason} — {duration_str}')
            description = f'Initiator: {initiator_name}\nStage Expires: {expires}\nCoordinators: {coordinator_str}\nActive stage mutes: {len(lines)}\n\n'+'\n'.join(lines)
            pages, chunk_size = [], 18
            for i in range(0, len(description.splitlines()), chunk_size):
                embed=discord.Embed(
                    title=f'\U0001F399 Stage info for {channel_obj.mention}',
                    description='\n'.join(description.splitlines()[i:i+chunk_size]),
                    color=discord.Color.purple()
                )
                pages.append(embed)
            paginator = UserPaginator(self.bot, interaction, pages)
            return await paginator.start()
            
    @commands.command(name='stages', help='Lists stage mute statistics.')
    @is_owner_developer_coordinator_predicator(None)
    async def list_stages_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in ctx.author.roles)
        channel_obj = await self.resolve_channel(ctx, target)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall(ctx)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to view stages in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if target and target.lower() == 'all':
                if not is_owner_or_dev:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Only owners/devs can view all stages.')
                stages = await conn.fetch('''
                    SELECT s.channel_id, s.initiator_id, s.expires_at, COUNT(v.discord_snowflake) AS active_mutes
                    FROM active_stages s
                    LEFT JOIN active_voice_mutes v
                        ON s.guild_id=v.guild_id AND s.channel_id=v.channel_id AND v.target='room'
                    WHERE s.guild_id=$1 AND s.room_name = $2
                    GROUP BY s.channel_id, s.initiator_id, s.expires_at
                    ORDER BY s.channel_id
                ''', ctx.guild.id, room_name)
                if not stages:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No active stages in {ctx.guild.name}.')
                pages, chunk_size = [], 8
                for i in range(0, len(stages), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F399 Active Stages in {ctx.guild.name}',
                        color=discord.Color.purple()
                    )
                    for s in stages[i:i+chunk_size]:
                        ch = ctx.guild.get_channel(s['channel_id'])
                        ch_name = ch.mention if ch else f'Unknown Channel ({s["channel_id"]})'
                        embed.add_field(name=ch_name, value=f'Active stage mutes: {s["active_mutes"]}', inline=False)
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            stage = await conn.fetchrow('''
                SELECT initiator_id, expires_at
                FROM active_stages
                WHERE guild_id=$1 AND channel_id=$2 AND room_name = $3
            ''', ctx.guild.id, channel_obj.id, room_name)
            if not stage:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No active stage in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            mutes = await conn.fetch('''
                SELECT discord_snowflake, expires_at, reason
                FROM active_voice_mutes
                WHERE guild_id=$1 AND channel_id=$2 AND target='room' and room_name = $3
            ''', ctx.guild.id, channel_obj.id, room_name)
            coordinators = await conn.fetch('''
                SELECT discord_snowflake
                FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2
            ''', ctx.guild.id, channel_obj.id)
            coordinator_mentions = []
            for c in coordinators:
                member = ctx.guild.get_member(c['discord_snowflake'])
                if member:
                    coordinator_mentions.append(member.mention)
            coordinator_str = ', '.join(coordinator_mentions) if coordinator_mentions else 'No coordinators'
            initiator = ctx.guild.get_member(stage['initiator_id'])
            initiator_name = initiator.mention if initiator else f'`{stage["initiator_id"]}`'
            expires = self.fmt_duration(stage['expires_at']) if stage['expires_at'] else 'No expiration'
            lines = []
            for m in mutes:
                user = ctx.guild.get_member(m['discord_snowflake'])
                duration_str = self.fmt_duration(m['expires_at']) if m['expires_at'] else 'No expiration'
                reason = m['reason'] or 'No reason provided'
                lines.append(f'• {user.mention} — {reason} — {duration_str}')
            description = (
                f'Initiator: {initiator_name}\n'
                f'Stage Expires: {expires}\n'
                f'Coordinators: {coordinator_str}\n'
                f'Active stage mutes: {len(lines)}\n\n' + '\n'.join(lines)
            )
            pages, chunk_size = [], 18
            for i in range(0, len(description.splitlines()), chunk_size):
                embed = discord.Embed(
                    title=f'\U0001F399 Stage info for {channel_obj.mention}',
                    description='\n'.join(description.splitlines()[i:i+chunk_size]),
                    color=discord.Color.purple()
                )
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
            
    @app_commands.command(name='survey', description='Survey moderators, developers, owners, and coordinators in the current or specified channel.')
    @app_commands.describe(channel='Tag a voice/stage channel')
    @is_owner_developer_app_predicator()
    async def stage_survey_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if not isinstance(channel_obj, discord.VoiceChannel):
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        owners, developers, moderators, coordinators = [], [], [], []
        for member in channel_obj.members:
            if await is_owner_member(member, self.bot): owners.append(member)
            elif await is_developer_member(member, self.bot): developers.append(member)
            elif await is_coordinator_via_objects(member, channel_obj): coordinators.append(member)
            elif await is_moderator_via_objects(member, channel_obj): moderators.append(member)
        def fmt(users): return ', '.join(u.mention for u in users) if users else '*None*'
        msg = (
            f'\U0001F50D **Survey results for {channel_obj.mention}:**\n'
            f'\n**Owners:** {fmt(owners)}'
            f'\n**Developers:** {fmt(developers)}'
            f'\n**Moderators:** {fmt(moderators)}'
            f'\n**Coordinators:** {fmt(coordinators)}'
            f'\n\nTotal surveyed: {len(channel_obj.members)}'
        )
        await self.send(interaction, content=msg, allowed_mentions=discord.AllowedMentions.none())
    
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
        channel_obj = await self.resolve_channel(ctx, channel)
        if not isinstance(channel_obj, (discord.VoiceChannel, discord.StageChannel)):
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid voice or stage channel.')
        owners, developers, moderators, coordinators = [], [], [], []
        for member in channel_obj.members:
            if await is_owner_member(member, ctx.bot): owners.append(member)
            elif await is_developer_member(member, ctx.bot): developers.append(member)
            elif await is_coordinator_via_objects(member, channel_obj): coordinators.append(member)
            elif await is_moderator_via_objects(member, channel_obj): moderators.append(member)
        def fmt(users):
            return ', '.join(u.mention for u in users) if users else '*None*'
        msg = (
            f'\U0001F50D **Survey results for {channel_obj.mention}:**\n'
            f'\n**Owners:** {fmt(owners)}'
            f'\n**Developers:** {fmt(developers)}'
            f'\n**Moderators:** {fmt(moderators)}'
            f'\n**Coordinators:** {fmt(coordinators)}'
            f'\n\nTotal surveyed: {len(channel_obj.members)}'
        )
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())

    @app_commands.command(name='temp', description='Mark a channel as a temporary room and assign an owner.')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID', owner='Tag a member or include their snowflake ID')
    @is_owner_developer_predicator()
    async def create_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None,
        owner: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if not isinstance(channel_obj, discord.VoiceChannel):
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        member_obj = await self.resolve_member_app(interaction, owner)
        temp_channel = TempChannel(channel_obj, channel_obj.name)
        self.temp_rooms.setdefault(interaction.guild.id, {})[channel_obj.name] = temp_channel
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (discord_snowflake, coordinator_room_names, coordinator_channel_ids, moderator_room_names, moderator_channel_ids)
                VALUES ($1, ARRAY[]::TEXT[], ARRAY[]::BIGINT[], ARRAY[]::TEXT[], ARRAY[]::BIGINT[])
                ON CONFLICT (discord_snowflake) DO NOTHING
            ''', member_obj.id)
            old_owner_row = await conn.fetchrow(
                'SELECT owner_snowflake, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_name=$2',
                interaction.guild.id, channel_obj.name
            )
            old_owner = old_owner_row['owner_snowflake'] if old_owner_row else None
            old_room_snowflake = old_owner_row['room_snowflake'] if old_owner_row else None
            if old_owner and old_owner != member_obj.id:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, coordinator_room_names, coordinator_channel_ids, moderator_room_names, moderator_channel_ids)
                    VALUES ($1, ARRAY[]::TEXT[], ARRAY[]::BIGINT[], ARRAY[]::TEXT[], ARRAY[]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO NOTHING
                ''', old_owner)
                channel_ids_to_remove = [ch for ch in (old_room_snowflake, channel_obj.id) if ch]
                for ch_id in channel_ids_to_remove:
                    await conn.execute('''
                        UPDATE users
                        SET coordinator_room_names = array_remove(coordinator_room_names, $1),
                            coordinator_channel_ids = array_remove(coordinator_channel_ids, $2),
                            moderator_room_names = array_remove(moderator_room_names, $1),
                            moderator_channel_ids = array_remove(moderator_channel_ids, $2),
                            updated_at = NOW()
                        WHERE discord_snowflake = $3
                    ''', channel_obj.name, ch_id, old_owner)
            await conn.execute('''
                UPDATE users
                SET moderator_room_names = array_remove(moderator_room_names, $1),
                    moderator_channel_ids = array_remove(moderator_channel_ids, $2),
                    updated_at = NOW()
                WHERE discord_snowflake = $3
            ''', channel_obj.name, channel_obj.id, member_obj.id)
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, room_name, owner_snowflake, room_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_snowflake, room_name) DO UPDATE SET owner_snowflake=$3, room_snowflake=$4
            ''', interaction.guild.id, channel_obj.name, member_obj.id, channel_obj.id)
            await conn.execute('''
                UPDATE users
                SET coordinator_room_names = CASE WHEN NOT $1=ANY(coordinator_room_names) THEN array_append(coordinator_room_names,$1) ELSE coordinator_room_names END,
                    coordinator_channel_ids = CASE WHEN NOT $2=ANY(coordinator_channel_ids) THEN array_append(coordinator_channel_ids,$2) ELSE coordinator_channel_ids END,
                    updated_at = NOW()
                WHERE discord_snowflake = $3
            ''', channel_obj.name, channel_obj.id, member_obj.id)
            await self.send(interaction, content=f'{self.get_random_emoji()} Channel {channel_obj.mention} is now owned by {member_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='temp', help='Mark a channel as a temporary room and assign an owner.')
    @is_owner_developer_predicator()
    async def create_temp_room_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        owner: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        if not channel:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please provide a valid channel or ID.')
        if not owner:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please provide the owner\'s Discord ID.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if not isinstance(channel_obj, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        member_obj = await self.resolve_member(ctx, owner)
        temp_channel = TempChannel(channel_obj, channel_obj.name)
        self.temp_rooms.setdefault(ctx.guild.id, {})[channel_obj.name] = temp_channel
        async with ctx.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (discord_snowflake, coordinator_room_names, coordinator_channel_ids, moderator_room_names, moderator_channel_ids)
                VALUES ($1, ARRAY[]::TEXT[], ARRAY[]::BIGINT[], ARRAY[]::TEXT[], ARRAY[]::BIGINT[])
                ON CONFLICT (discord_snowflake) DO NOTHING
            ''', member_obj.id)
            old_owner_row = await conn.fetchrow(
                'SELECT owner_snowflake, room_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_name=$2',
                ctx.guild.id, channel_obj.name
            )
            old_owner = old_owner_row['owner_snowflake'] if old_owner_row else None
            old_room_snowflake = old_owner_row['room_snowflake'] if old_owner_row else None
            if old_owner and old_owner != member_obj.id:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, coordinator_room_names, coordinator_channel_ids, moderator_room_names, moderator_channel_ids)
                    VALUES ($1, ARRAY[]::TEXT[], ARRAY[]::BIGINT[], ARRAY[]::TEXT[], ARRAY[]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO NOTHING
                ''', old_owner)
                channel_ids_to_remove = [ch for ch in (old_room_snowflake, channel_obj.id) if ch]
                for ch_id in channel_ids_to_remove:
                    await conn.execute('''
                        UPDATE users
                        SET coordinator_room_names = array_remove(coordinator_room_names, $1),
                            coordinator_channel_ids = array_remove(coordinator_channel_ids, $2),
                            moderator_room_names = array_remove(moderator_room_names, $1),
                            moderator_channel_ids = array_remove(moderator_channel_ids, $2),
                            updated_at = NOW()
                        WHERE discord_snowflake = $3
                    ''', channel_obj.name, ch_id, old_owner)
            await conn.execute('''
                UPDATE users
                SET moderator_room_names = array_remove(moderator_room_names, $1),
                    moderator_channel_ids = array_remove(moderator_channel_ids, $2),
                    updated_at = NOW()
                WHERE discord_snowflake = $3
            ''', channel_obj.name, channel_obj.id, member_obj.id)
            await conn.execute('''
                INSERT INTO temporary_rooms (guild_snowflake, room_name, owner_snowflake, room_snowflake)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_snowflake, room_name) DO UPDATE SET owner_snowflake=$3, room_snowflake=$4
            ''', ctx.guild.id, channel_obj.name, member_obj.id, channel_obj.id)
            await conn.execute('''
                UPDATE users
                SET coordinator_room_names = CASE WHEN NOT $1=ANY(coordinator_room_names) THEN array_append(coordinator_room_names,$1) ELSE coordinator_room_names END,
                    coordinator_channel_ids = CASE WHEN NOT $2=ANY(coordinator_channel_ids) THEN array_append(coordinator_channel_ids,$2) ELSE coordinator_channel_ids END,
                    updated_at = NOW()
                WHERE discord_snowflake = $3
            ''', channel_obj.name, channel_obj.id, member_obj.id)
            await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Channel {channel_obj.mention} is now owned by {member_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())

    @app_commands.command(name='temps', description='List temporary rooms with matching command aliases.')
    @is_owner_developer_predicator()
    async def check_temp_rooms_app_command(
        self,
        interaction: discord.Interaction
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in servers.')
        aliases = self.bot.command_aliases.get(interaction.guild.id, {})
        rooms = self.temp_rooms.get(interaction.guild.id, {})
        temp_aliases = aliases.get('temp_room_aliases', {})
        if not rooms:
            return await self.send(interaction, content='No temporary rooms found.')
        listings = []
        for room_name, room in rooms.items():
            room_id = room.channel.id
            listings.append(f"{room_name} ({room_id})")
            for alias_type, group in temp_aliases.items():
                for alias_name,meta in group.items():
                    a_id = meta.get('channel_id')
                    a_room = meta.get('room_name')
                    if a_id == room_id and a_room == room_name:
                        listings.append(f"  ↳ {alias_name} ({alias_type})")
        output = '\n'.join(listings)
        return await self.send(interaction, content=f'{self.get_random_emoji()}\n{output}')
        
    @commands.command(name='temps', help='List temporary rooms with matching command aliases.')
    @is_owner_developer_predicator()
    async def check_temp_rooms_text_command(
        self,
        ctx: commands.Context
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        aliases = self.bot.command_aliases.get(ctx.guild.id, {})
        rooms = self.temp_rooms.get(ctx.guild.id,{})
        temp_aliases = aliases.get('temp_room_aliases', {})
        if not rooms:
            return await self.handler.send_message(ctx, content='\U0001F6AB No temporary rooms found.')
        listings = []
        for room_name, room in rooms.items():
            room_id = room.channel.id
            listings.append(f"{room_name} ({room_id})")
            for alias_type, group in temp_aliases.items():
                for alias_name,meta in group.items():
                    a_id = meta.get('channel_id')
                    a_room = meta.get('room_name')
                    if a_id == room_id and a_room == room_name:
                        listings.append(f"  ↳ {alias_name} ({alias_type})")
        output = '\n'.join(listings)
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()}\n{output}')
        
    @app_commands.command(name='tmutes', description='Lists text-mute statistics.')
    @app_commands.describe(target='"all", channel name/ID/mention, or user mention/ID')
    @is_owner_developer_coordinator_moderator_role_app_predicator(None)
    async def list_text_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in interaction.user.roles)
        member_obj = await self.resolve_member(interaction, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(interaction, target)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        if member_obj:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall_app(interaction)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use this command (`tmutes`) for {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if target and target.lower()=='all' and is_owner_or_dev:
                records = await conn.fetch('''SELECT discord_snowflake, channel_id, room_name, reason, expires_at FROM active_text_mutes WHERE guild_id = $1 ORDER BY room_name, channel_id, discord_snowflake''', interaction.guild.id)
                if not records:
                    return await self.send(interaction, content=f'\U0001F6AB No users are currently text-muted in {interaction.guild.name}.')
                grouped = defaultdict(list)
                for r in records:
                    grouped[r['channel_id']].append(r)
                pages, chunk_size = [], 18
                for ch_id, entries in sorted(grouped.items()):
                    ch = interaction.guild.get_channel(ch_id)
                    ch_name = ch.mention
                    for i in range(0, len(entries), chunk_size):
                        embed= discord.Embed(
                            title=f'\U0001F4DA Text-mute records for {ch_name}',
                            color=discord.Color.orange()
                        )
                        for e in entries[i:i+chunk_size]:
                            user = interaction.guild.get_member(e['discord_snowflake'])
                            mention = user.name if user else f'`{e["discord_snowflake"]}`'
                            reason = e['reason'] or 'No reason provided'
                            duration_str = self.fmt_duration(e['expires_at'])
                            embed.add_field(name=mention, value=f'Reason: {reason}\nDuration: {duration_str}', inline=False)
                        pages.append(embed)
                paginator = UserPaginator(self.bot, interaction, pages)
                return await paginator.start()
            elif member_obj:
                records = await conn.fetch('''SELECT channel_id, room_name, reason, expires_at FROM active_text_mutes WHERE discord_snowflake = $1 AND guild_id = $2 AND room_name = $3''', member_obj.id, interaction.guild.id, room_name)
                if not records: return await self.send(interaction, content=f'\U0001F6AB {member_obj.mention} is not text-muted in any channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in records:
                    ch = interaction.guild.get_channel(r['channel_id'])
                    duration_str = self.fmt_duration(r['expires_at'])
                    lines.append(f'• {ch.mention} — {r["reason"]} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F4DA Text-mute records for {member_obj.display_name}',
                        description='\n'.join(lines[i:i+chunk_size]),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = UserPaginator(self.bot, interaction, pages)
                return await paginator.start()
            elif channel_obj:
                records = await conn.fetch('''SELECT discord_snowflake, room_name, reason, expires_at FROM active_text_mutes WHERE channel_id = $1 AND guild_id = $2 AND room_name = $3''', channel_obj.id, interaction.guild.id, room_name)
                if not records:
                    return await self.send(interaction, content=f'\U0001F6AB No users are currently text-muted in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in records:
                    user = interaction.guild.get_member(r['discord_snowflake'])
                    if not user:
                        continue
                    duration_str = self.fmt_duration(r['expires_at'])
                    lines.append(f'• {user.mention} — {r["reason"]} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F4DA Text-mute records for {channel_obj.name}',
                        description='\n'.join(lines[i:i+chunk_size]),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = UserPaginator(self.bot, interaction, pages)
                return await paginator.start()
        return await self.send(interaction, content='\U0001F6AB You must specify "all", a member, or a text channel.')
    
    @commands.command(name='tmutes', help='Lists text-mute statistics.')
    @is_owner_developer_coordinator_moderator_role_predicator(None)
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        valid_role_ids = [r['role_id'] for r in rows]
        is_team_member = any(r.id in valid_role_ids for r in ctx.author.roles)
        member_obj = await self.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.resolve_channel(ctx, target)
        room_name = ''
        temp_room_obj = None
        if channel_obj:
            for guild_temp_rooms in self.temp_rooms.values():
                for temp_channel in guild_temp_rooms.values():
                    if temp_channel.channel.id == channel_obj.id:
                        temp_room_obj = temp_channel
                        room_name = temp_channel.room_name
                        break
                if temp_room_obj:
                    break
        if member_obj:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_overall(ctx)
        else:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord and not is_team_member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`tmutes`) in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if target and target.lower() == 'all' and is_owner_or_dev:
                records = await conn.fetch('''
                    SELECT discord_snowflake, channel_id, room_name, reason, expires_at
                    FROM active_text_mutes
                    WHERE guild_id = $1
                    ORDER BY room_name, channel_id, discord_snowflake
                ''', ctx.guild.id)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are currently text-muted in {ctx.guild.name}.')
                grouped = defaultdict(list)
                for r in records:
                    key = (r['channel_id'], r['room_name'])
                    grouped[key].append(r)
                pages, chunk_size = [], 18
                for (ch_id, rname), entries in sorted(grouped.items()):
                    ch = ctx.guild.get_channel(ch_id)
                    ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                    if rname:
                        ch_name = f"{ch_name}"
                    for i in range(0, len(entries), chunk_size):
                        embed = discord.Embed(
                            title=f'\U0001F4DA Text-mute records for {ch_name}',
                            color=discord.Color.orange()
                        )
                        for e in entries[i:i + chunk_size]:
                            user = ctx.guild.get_member(e['discord_snowflake'])
                            mention = user.name if user else f'<@{e["discord_snowflake"]}>'
                            reason = e['reason'] or 'No reason provided'
                            duration_str = self.fmt_duration(e['expires_at'])
                            embed.add_field(name=mention, value=f'Reason: {reason}\nDuration: {duration_str}', inline=False)
                        pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            elif member_obj:
                records = await conn.fetch('''
                    SELECT channel_id, room_name, reason, expires_at
                    FROM active_text_mutes
                    WHERE discord_snowflake = $1 AND guild_id = $2
                ''', member_obj.id, ctx.guild.id)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not text-muted in any channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in records:
                    ch = ctx.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'Channel ID `{r["channel_id"]}`'
                    if r['room_name']:
                        ch_name = f"{ch_name}"
                    duration_str = self.fmt_duration(r['expires_at'])
                    lines.append(f'• {ch_name} — {r["reason"]} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F4DA Text-mute records for {member_obj.display_name}',
                        description='\n'.join(lines[i:i+chunk_size]),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            elif channel_obj:
                records = await conn.fetch('''
                    SELECT discord_snowflake, reason, expires_at
                    FROM active_text_mutes
                    WHERE channel_id = $1 AND guild_id = $2 AND room_name = $3
                ''', channel_obj.id, ctx.guild.id, room_name)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are currently text-muted in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in records:
                    user = ctx.guild.get_member(r['discord_snowflake'])
                    mention = user.mention if user else f'<@{r["discord_snowflake"]}>'
                    duration_str = self.fmt_duration(r['expires_at'])
                    lines.append(f'• {mention} — {r["reason"]} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F4DA Text-mute records for {channel_obj.name}',
                        description='\n'.join(lines[i:i+chunk_size]),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F6AB You must specify "all", a member, or a text channel.')
        
    @app_commands.command(name='trole', description='Adds or enables a team role in the permission table.')
    @is_owner_predicator()
    @app_commands.describe(role_ref='Tag a role or provide its ID')
    async def team_role_interaction(
        self,
        interaction: discord.Interaction,
        *,
        role_ref: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        if not role_ref:
            return await self.send(interaction, content='\U0001F6AB You must mention a role or provide its ID.')
        role = None
        if interaction.data.get('resolved', {}).get('roles'):
            roles_resolved = interaction.data['resolved']['roles']
            role_id = next(iter(roles_resolved))
            role = interaction.guild.get_role(int(role_id))
        else:
            try:
                role_id = int(role_ref.replace('<@&','').replace('>',''))
                role = interaction.guild.get_role(role_id)
            except ValueError:
                pass
        if not role:
            return await self.send(interaction, content='\U0001F6AB Could not find that role. Please mention it or use its numeric ID.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO role_permissions (role_id, is_team_member)
                VALUES ($1, TRUE)
                ON CONFLICT (role_id)
                DO UPDATE SET is_team_member = TRUE;
            ''', role.id)
        await self.send(interaction, content=f'{self.get_random_emoji()} Role `{role.name}` (`{role.id}`) is now authorized for moderator-level checks.')

    @commands.command(name='trole', help='Adds or enables a team role in the permission table.')
    @is_owner_predicator()
    async def team_role_text_command(
        self,
        ctx: commands.Context,
        *,
        role_ref: Optional[str] = commands.parameter(description='Tag a role or provide its ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        if not role_ref:
            return await self.handler.send_message(ctx, content='\U0001F6AB You must mention a role or provide its ID.')
        role = None
        if ctx.message.role_mentions:
            role = ctx.message.role_mentions[0]
        else:
            try:
                role_id = int(role_ref)
                role = ctx.guild.get_role(role_id)
            except ValueError:
                pass
        if not role:
            return await self.handler.send_message(ctx, content='\U0001F6AB Could not find that role. Please mention it or use its numeric ID.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO role_permissions (role_id, is_team_member)
                VALUES ($1, TRUE)
                ON CONFLICT (role_id)
                DO UPDATE SET is_team_member = TRUE;
            ''', role.id)
        await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Role `{role.name}` (`{role.id}`) is now authorized for moderator-level checks.')
            
    @app_commands.command(name='xadmin', description='Revokes server mute privileges from a user.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_app_predicator()
    async def demote_server_muter_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        member_obj = await self.resolve_member_app(interaction, member)
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await self.send(interaction, content='\U0001F6AB You cannot revoke server mute for the bot.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET server_muter_guild_ids = (SELECT ARRAY(
                    SELECT unnest(server_muter_guild_ids)
                     EXCEPT SELECT $2
                ))
                WHERE discord_snowflake = $1
            ''', member_obj.id, interaction.guild.id)
        await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} no longer has server mute privileges.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='xadmin', help='Revokes server mute privileges from a user.')
    @is_owner_predicator()
    async def demote_server_muter_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj or not member:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
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
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} no longer has server mute privileges.', allowed_mentions=discord.AllowedMentions.none())
        
    @app_commands.command(name='xalias', description='Deletes an alias.')
    @is_owner_developer_coordinator_app_predicator(None)
    @app_commands.describe(alias_name='Include an alias name')
    async def delete_alias_app_command(
        self,
        interaction: discord.Interaction,
        alias_name: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        guild_aliases = self.bot.command_aliases.get(interaction.guild.id, self.bot.command_aliases.default_factory())
        alias_dict = None
        alias_entry = None
        alias_room = None
        alias_type = None
        for c in ('cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute'):
            if alias_name in guild_aliases.get('channel_aliases', {}).get(c, {}):
                alias_dict = 'channel_aliases'
                alias_type = c
                alias_entry = guild_aliases['channel_aliases'][c][alias_name]
                break
        if not alias_type:
            for c in ('role', 'unrole'):
                if alias_name in guild_aliases.get('role_aliases', {}).get(c, {}):
                    alias_dict = 'role_aliases'
                    alias_type = c
                    alias_entry = guild_aliases['role_aliases'][c][alias_name]
                    break
        if not alias_type:
            for c in ('mute', 'unmute', 'ban', 'unban', 'tmute', 'untmute'):
                t = guild_aliases.get('temp_room_aliases', {}).get(c, {})
                if alias_name in t:
                    alias_dict = 'temp_room_aliases'
                    alias_type = c
                    alias_entry = t[alias_name]
                    alias_room = alias_entry.get('room_name')
                    break
        if not alias_type:
            return await self.send(interaction, content=f'\U0001F6AB Alias `{alias_name}` not found.')
        channel_obj = None
        if alias_dict in ('channel_aliases', 'temp_room_aliases'):
            cid = alias_entry.get('channel_id') if isinstance(alias_entry, dict) else alias_entry
            channel_obj = interaction.guild.get_channel(cid) if cid else None
            if channel_obj:
                is_owner_or_dev, _ = await check_owner_dev_coord_mod_app(interaction, channel_obj)
                if not is_owner_or_dev:
                    async with self.bot.db_pool.acquire() as conn:
                        row = await conn.fetchrow('SELECT coordinator_channel_ids FROM users WHERE discord_snowflake=$1', interaction.user.id)
                    allowed = row.get('coordinator_channel_ids') or []
                    if cid not in allowed:
                        return await self.send(interaction, content=f'\U0001F6AB You do not have permission to delete alias `{alias_name}` in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if alias_dict == 'temp_room_aliases' and alias_room:
                await conn.execute('DELETE FROM command_aliases WHERE guild_id=$1 AND alias_type=$2 AND alias_name=$3 AND room_name=$4', interaction.guild.id, alias_type, alias_name, alias_room)
            else:
                await conn.execute('DELETE FROM command_aliases WHERE guild_id=$1 AND alias_type=$2 AND alias_name=$3', interaction.guild.id, alias_type, alias_name)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1,$2,$3,$4,$5,$6)
                ''', 'delete_alias', None, interaction.user.id, interaction.guild.id, cid, f'Deleted alias {alias_name}')
        if self.bot.get_command(alias_name):
            self.bot.remove_command(alias_name)
        if alias_dict == 'temp_room_aliases':
            guild_aliases['temp_room_aliases'][alias_type].pop(alias_name, None)
        elif alias_dict == 'channel_aliases':
            guild_aliases['channel_aliases'][alias_type].pop(alias_name, None)
        else:
            guild_aliases['role_aliases'][alias_type].pop(alias_name, None)
        return await self.send(interaction, content=f'{self.get_random_emoji()} Deleted alias `{alias_name}` from `{alias_type}`.')

    @commands.command(name='xalias', help='Deletes an alias.')
    @is_owner_developer_coordinator_predicator(None)
    async def delete_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: Optional[str] = commands.parameter(default=None, description='Include an alias name')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        if not alias_name or not alias_name.strip():
            return await self.handler.send_message(ctx, content='\U0001F6AB `alias_name` cannot be empty.')
        guild_aliases = self.bot.command_aliases.get(ctx.guild.id, {})
        alias_dict = None
        alias_entry = None
        alias_room = None
        alias_type = None
        for c in ('cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute'):
            if alias_name in guild_aliases.get('channel_aliases', {}).get(c, {}):
                alias_dict = 'channel_aliases'
                alias_type = c
                alias_entry = guild_aliases['channel_aliases'][c][alias_name]
                break
        if not alias_type:
            for c in ('role', 'unrole'):
                if alias_name in guild_aliases.get('role_aliases', {}).get(c, {}):
                    alias_dict = 'role_aliases'
                    alias_type = c
                    alias_entry = guild_aliases['role_aliases'][c][alias_name]
                    break
        if not alias_type:
            for c in ('mute', 'unmute', 'ban', 'unban', 'tmute', 'untmute'):
                t = guild_aliases.get('temp_room_aliases', {}).get(c, {})
                if alias_name in t:
                    alias_dict = 'temp_room_aliases'
                    alias_type = c
                    alias_entry = t[alias_name]
                    alias_room = alias_entry.get('room_name')
                    break
        if not alias_type:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Alias `{alias_name}` not found.')
        channel_obj = None
        cid = alias_entry.get('channel_id') if isinstance(alias_entry, dict) else alias_entry
        if alias_dict in ('channel_aliases', 'temp_room_aliases'):
            channel_obj = ctx.guild.get_channel(cid) if cid else None
            if channel_obj:
                is_owner_or_dev, _ = await check_owner_dev_coord_mod(ctx, channel_obj)
                if not is_owner_or_dev:
                    async with ctx.bot.db_pool.acquire() as conn:
                        row = await conn.fetchrow('SELECT coordinator_channel_ids FROM users WHERE discord_snowflake=$1', ctx.author.id)
                    allowed = row.get('coordinator_channel_ids') or []
                    if cid not in allowed:
                        return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to delete alias `{alias_name}` in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            if alias_dict == 'temp_room_aliases' and alias_room:
                await conn.execute('DELETE FROM command_aliases WHERE guild_id=$1 AND alias_type=$2 AND alias_name=$3 AND room_name=$4', ctx.guild.id, alias_type, alias_name, alias_room)
            else:
                await conn.execute('DELETE FROM command_aliases WHERE guild_id=$1 AND alias_type=$2 AND alias_name=$3', ctx.guild.id, alias_type, alias_name)
            await conn.execute('''
                INSERT INTO moderation_logs(action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES($1,$2,$3,$4,$5,$6)
            ''', 'delete_alias', None, ctx.author.id, ctx.guild.id, cid, f'Deleted alias {alias_name}')
        if self.bot.get_command(alias_name):
            self.bot.remove_command(alias_name)
        if alias_dict == 'temp_room_aliases':
            guild_aliases['temp_room_aliases'][alias_type].pop(alias_name, None)
        elif alias_dict == 'channel_aliases':
            guild_aliases['channel_aliases'][alias_type].pop(alias_name, None)
        else:
            guild_aliases['role_aliases'][alias_type].pop(alias_name, None)
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Deleted alias `{alias_name}` from `{alias_type}`.')
    
    @app_commands.command(name='xcap', description='Delete a duration cap for bans, mutes and text mutes.')
    @is_owner_developer_predicator()
    @app_commands.describe(channel='Tag a channel or include its snowflake ID', moderation_type='One of: `mute`, `ban`, `tmute`')
    async def undo_cap_interaction(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None,
        moderation_type: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        is_temp_room = False
        room_name = ''
        temp_rooms = self.temp_rooms.get(interaction.guild.id, {})
        for temp_channel in temp_rooms.values():
            if temp_channel.room_name == channel_obj.name:
                channel_obj = temp_channel
                is_temp_room = True
                room_name = temp_channel.room_name
                break
        valid_types = {'mute','ban','tmute'}
        if moderation_type not in valid_types:
            return await self.send(interaction, content=f'\U0001F6AB Invalid moderation type. Must be one of: {", ".join(valid_types)}')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4',
                interaction.guild.id, channel_obj.id, moderation_type, room_name
            )
            if is_temp_room:
                room_display = room_name
            else:
                room_display = channel_obj.mention
            if row:
                original_duration = row['duration']
                msg = f'{self.get_random_emoji()} Cap deleted on {room_display} for {moderation_type} of duration {original_duration}.'
            else:
                return await self.send(interaction, content=f'\U0001F6AB No caps exist in {room_display} for {moderation_type}.')
            await conn.execute('''
                DELETE FROM active_caps
                WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4
            ''', interaction.guild.id, channel_obj.id, moderation_type, room_name)
        await self.send(interaction, content=msg, allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='xcap', help='Delete a duration cap for bans, mutes and text mutes.')
    @is_owner_developer_predicator()
    async def undo_cap_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        moderation_type: Optional[str] = commands.parameter(default=None, description='One of: `mute`, `ban`, `tmute`')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        is_temp_room = False
        room_name = ''
        temp_rooms = self.temp_rooms.get(ctx.guild.id, {})
        for temp_channel in temp_rooms.values():
            if temp_channel.room_name == channel_obj.name:
                channel_obj = temp_channel
                is_temp_room = True
                room_name = temp_channel.room_name
                break
        valid_types = {'mute', 'ban', 'tmute'}
        if moderation_type not in valid_types:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid moderation type. Must be one of: {", ".join(valid_types)}')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4',
                ctx.guild.id, channel_obj.id, moderation_type, room_name
            )
            if is_temp_room:
                room_display = room_name
            else:
                room_display = channel_obj.mention
            if row:
                original_duration = row['duration']
                msg = f'{self.get_random_emoji()} Cap deleted on {room_display} for {moderation_type} of duration {original_duration}.'
            else:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No caps exist in {room_display} for {moderation_type}.')
            await conn.execute('''
                DELETE FROM active_caps
                WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4
            ''', ctx.guild.id, channel_obj.id, moderation_type, room_name)
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())
        
    @app_commands.command(name='xcoord', description='Revokes coordinator access from a user in a specific voice channel or temporary room.')
    @is_owner_developer_app_predicator()
    @app_commands.describe(member='Tag a member or include their snowflake ID', channel='Tag a channel or include its snowflake ID')
    async def demote_coordinator_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        member_obj = await self.resolve_member_app(interaction, member)
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT coordinator_channel_ids, coordinator_room_names FROM users WHERE discord_snowflake = $1', member_obj.id)
            if not row:
                return await self.send(interaction, content=f'\U0001F6AB {member_obj.mention} is not found in the coordinator database.', allowed_mentions=discord.AllowedMentions.none())
            current_channel_ids = row.get('coordinator_channel_ids', []) or []
            current_rooms = row.get('coordinator_room_names', []) or []
            if channel_obj.id in current_channel_ids:
                await conn.execute('''
                    UPDATE users
                    SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $2),
                        updated_at = NOW()
                    WHERE discord_snowflake = $1
                ''', member_obj.id, channel_obj.id)
            temp_rooms_to_remove = []
            for room in current_rooms:
                temp_alias = self.bot.command_aliases.get(interaction.guild.id, self.bot.command_aliases.default_factory()) \
                    .get('temp_room_aliases', {}) \
                    .get('mute', {}) \
                    .get(room)
                if temp_alias and temp_alias.get('channel_id') == (channel_obj.id if channel_obj else -1):
                    temp_rooms_to_remove.append(room)
            if temp_rooms_to_remove:
                for room_name in temp_rooms_to_remove:
                    await conn.execute('''
                        UPDATE users
                        SET coordinator_room_names = array_remove(coordinator_room_names, $2),
                            updated_at = NOW()
                        WHERE discord_snowflake = $1
                    ''', member_obj.id, room_name)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1, $2, $3, $4, $5, $6)
            ''', 'remove_coordinator', member_obj.id, interaction.user.id, interaction.guild.id, channel_obj.id if channel_obj else None, 'Removed a coordinator from a voice channel/temp room')
            updated_row = await conn.fetchrow('SELECT coordinator_channel_ids, coordinator_room_names FROM users WHERE discord_snowflake = $1', member_obj.id)
            remaining_channels = updated_row.get('coordinator_channel_ids', []) if updated_row else []
            remaining_rooms = updated_row.get('coordinator_room_names', []) if updated_row else []
            guild_voice_channel_ids = [vc.id for vc in interaction.guild.voice_channels]
            has_remaining_guild_channels = any(ch_id in guild_voice_channel_ids for ch_id in remaining_channels)
            has_remaining_rooms = bool(remaining_rooms)
            if not has_remaining_guild_channels and not has_remaining_rooms:
                return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention}\'s coordinator access has been completely revoked from {channel_obj.mention} and in {interaction.guild.name} (no remaining channels or temp rooms).', allowed_mentions=discord.AllowedMentions.none())
            else:
                return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention}\'s coordinator access has been revoked from {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='xcoord', help='Revokes coordinator access from a user in a specific voice channel.')
    @is_owner_developer_predicator()
    async def demote_coordinator_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT coordinator_channel_ids, coordinator_room_names FROM users WHERE discord_snowflake=$1', member_obj.id)
            if not row:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not found in the coordinator database.', allowed_mentions=discord.AllowedMentions.none())
            current_channel_ids = row.get('coordinator_channel_ids', []) or []
            current_rooms = row.get('coordinator_room_names', []) or []
            if channel_obj.id in current_channel_ids:
                await conn.execute('UPDATE users SET coordinator_channel_ids=array_remove(coordinator_channel_ids,$2), updated_at=NOW() WHERE discord_snowflake=$1', member_obj.id, channel_obj.id)
            temp_rooms = await conn.fetch('SELECT room_name FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2', ctx.guild.id, channel_obj.id)
            for tr in temp_rooms:
                room_name = tr['room_name']
                if room_name in current_rooms:
                    await conn.execute('UPDATE users SET coordinator_room_names=array_remove(coordinator_room_names,$2), updated_at=NOW() WHERE discord_snowflake=$1', member_obj.id, room_name)
            await conn.execute('INSERT INTO moderation_logs(action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES($1,$2,$3,$4,$5,$6)', 'remove_coordinator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Removed a coordinator from a voice channel/temp room')
            updated_row = await conn.fetchrow('SELECT coordinator_channel_ids, coordinator_room_names FROM users WHERE discord_snowflake=$1', member_obj.id)
            remaining_channels = updated_row.get('coordinator_channel_ids', []) if updated_row else []
            remaining_rooms = updated_row.get('coordinator_room_names', []) if updated_row else []
            guild_voice_channel_ids = [vc.id for vc in ctx.guild.voice_channels]
            has_remaining_guild_channels = any(ch_id in guild_voice_channel_ids for ch_id in remaining_channels)
            has_remaining_rooms = bool(remaining_rooms)
            if not has_remaining_guild_channels and not has_remaining_rooms:
                return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention}\'s coordinator access has been completely revoked from {channel_obj.mention} and in {ctx.guild.name} (no remaining channels or temp rooms).', allowed_mentions=discord.AllowedMentions.none())
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention}\'s coordinator access has been revoked from {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                
    @app_commands.command(name='xdev', description='Removes a developer.')
    @is_owner_app_predicator()
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    async def delete_developer_app(
        self,
        interaction: discord.Interaction,
        member: str = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        member_obj = await self.resolve_member_app(interaction, member)
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from input: {member}')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET developer_guild_ids = array_remove(developer_guild_ids, $2),
                    updated_at = NOW()
                WHERE discord_snowflake = $1
            ''', member_obj.id, interaction.guild.id)
        return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention}\'s developer access has been revoked in {interaction.guild.name}.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='xdev', help='Removes a developer.')
    @is_owner_predicator()
    async def delete_developer_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
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
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention}\'s developer access has been revoked in {ctx.guild.name}.', allowed_mentions=discord.AllowedMentions.none())
    
    @app_commands.command(name='xmod', description="Revokes a member's VC moderator role for a given channel.")
    @is_owner_developer_coordinator_app_predicator(None)
    @app_commands.describe(member='Tag a member or include their snowflake ID', channel='Tag a channel or include its snowflake ID')
    async def delete_moderator_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        member_obj = await self.resolve_member_app(interaction, member)
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod_app(interaction, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.send(interaction, content=f'\U0001F6AB You do not have permission to use this command (`xmod`) in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT moderator_channel_ids, moderator_room_names FROM users WHERE discord_snowflake = $1', member_obj.id)
            if not row:
                return await self.send(interaction, content=f'\U0001F6AB {member_obj.mention} is not found in the moderator database.', allowed_mentions=discord.AllowedMentions.none())
            current_channel_ids = row.get('moderator_channel_ids', []) or []
            current_rooms = row.get('moderator_room_names', []) or []
            if channel_obj.id in current_channel_ids:
                await conn.execute('UPDATE users SET moderator_channel_ids = array_remove(moderator_channel_ids, $2), updated_at = NOW() WHERE discord_snowflake = $1', member_obj.id, channel_obj.id)
            temp_rooms_to_remove = []
            for room in current_rooms:
                temp_alias = self.bot.command_aliases.get(guild.id, self.bot.command_aliases.default_factory()).get('temp_room_aliases', {}).get('mute', {}).get(room)
                if temp_alias and temp_alias.get('channel_id') == channel_obj.id:
                    temp_rooms_to_remove.append(room)
            for room_name in temp_rooms_to_remove:
                await conn.execute('UPDATE users SET moderator_room_names = array_remove(moderator_room_names, $2), updated_at = NOW() WHERE discord_snowflake = $1', member_obj.id, room_name)
            await conn.execute('INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1, $2, $3, $4, $5, $6)', 'remove_moderator', member_obj.id, interaction.user.id, interaction.guild.id, channel_obj.id, 'Removed a moderator from the channel/temp room')
            updated_row = await conn.fetchrow('SELECT moderator_channel_ids, moderator_room_names FROM users WHERE discord_snowflake = $1', member_obj.id)
            remaining_channels = updated_row.get('moderator_channel_ids', []) if updated_row else []
            remaining_rooms = updated_row.get('moderator_room_names', []) if updated_row else []
            guild_voice_channel_ids = [vc.id for vc in interaction.guild.voice_channels]
            has_remaining_guild_channels = any(ch_id in guild_voice_channel_ids for ch_id in remaining_channels)
            has_remaining_rooms = bool(remaining_rooms)
            if not has_remaining_guild_channels and not has_remaining_rooms:
                return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has had all moderator access revoked in {channel_obj.mention} and {interaction.guild.name} (no remaining channels or temp rooms).', allowed_mentions=discord.AllowedMentions.none())
            else:
                return await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has been revoked moderator access in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                
    @commands.command(name='xmod', help='Revokes a member\'s VC moderator role for a given channel.')
    @is_owner_developer_coordinator_predicator(None)
    async def delete_moderator_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`xmod`) in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT moderator_channel_ids, moderator_room_names FROM users WHERE discord_snowflake=$1', member_obj.id)
            if not row:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not found in the moderator database.', allowed_mentions=discord.AllowedMentions.none())
            current_channel_ids = row.get('moderator_channel_ids', []) or []
            current_rooms = row.get('moderator_room_names', []) or []
            if channel_obj.id in current_channel_ids:
                await conn.execute('UPDATE users SET moderator_channel_ids=array_remove(moderator_channel_ids,$2), updated_at=NOW() WHERE discord_snowflake=$1', member_obj.id, channel_obj.id)
            temp_rooms = await conn.fetch('SELECT room_name FROM temporary_rooms WHERE guild_snowflake=$1 AND room_snowflake=$2', ctx.guild.id, channel_obj.id)
            for tr in temp_rooms:
                room_name = tr['room_name']
                if room_name in current_rooms:
                    await conn.execute('UPDATE users SET moderator_room_names=array_remove(moderator_room_names,$2), updated_at=NOW() WHERE discord_snowflake=$1', member_obj.id, room_name)
            await conn.execute('INSERT INTO moderation_logs(action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES($1,$2,$3,$4,$5,$6)', 'remove_moderator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Removed a moderator from the channel/temp room')
            updated_row = await conn.fetchrow('SELECT moderator_channel_ids, moderator_room_names FROM users WHERE discord_snowflake=$1', member_obj.id)
            remaining_channels = updated_row.get('moderator_channel_ids', []) if updated_row else []
            remaining_rooms = updated_row.get('moderator_room_names', []) if updated_row else []
            guild_voice_channel_ids = [vc.id for vc in ctx.guild.voice_channels]
            has_remaining_guild_channels = any(ch_id in guild_voice_channel_ids for ch_id in remaining_channels)
            has_remaining_rooms = bool(remaining_rooms)
            if not has_remaining_guild_channels and not has_remaining_rooms:
                return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has had all moderator access revoked from {channel_obj.mention} and {ctx.guild.name} (no remaining channels or temp rooms).', allowed_mentions=discord.AllowedMentions.none())
            return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been revoked moderator access in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
    
    @app_commands.command(name='xrmute', description='Unmutes all members in a VC (except yourself).')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_app_predicator()
    async def room_unmute_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content='\U0001F6AB Please specify a valid target.')
        unmuted_members, skipped_members, failed_members = [], [], []
        async with self.bot.db_pool.acquire() as conn:
            for member in channel_obj.members:
                try:
                    row = await conn.fetchrow('''
                        SELECT 1 FROM active_voice_mutes
                        WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target='user'
                        AND (expires_at IS NULL OR expires_at > NOW())
                    ''', interaction.guild.id, member.id, channel_obj.id)
                    if not row:
                        skipped_members.append(member)
                        continue
                    if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                        await member.edit(mute=False)
                    await conn.execute('''
                        DELETE FROM active_voice_mutes
                        WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target='user'
                    ''', interaction.guild.id, member.id, channel_obj.id)
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'unmute', member.id, interaction.user.id, interaction.guild.id, channel_obj.id, 'Unmuted via room_unmute')
                    unmuted_members.append(member)
                except Exception as e:
                    logger.warning(f'Unmute failed for {member.mention}: {e}')
                    failed_members.append(member)
        summary = f'{self.get_random_emoji()} Unmuted {len(unmuted_members)} member(s) in {channel_obj.mention}.'
        if skipped_members:
            summary += f'\n\U000026A0\U0000FE0F Skipped {len(skipped_members)}.'
        if failed_members:
            summary += f'\n\U0001F6AB Failed to unmute {len(failed_members)}.'
        return await self.send(interaction, content=summary, allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='xrmute', help='Unmutes all members in a VC (except yourself).')
    @is_owner_predicator()
    async def room_unmute_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        unmuted_members, skipped_members, failed_members = [], [], []
        async with self.bot.db_pool.acquire() as conn:
            for member in channel_obj.members:
                if member.id == ctx.author.id:
                    continue
                try:
                    row = await conn.fetchrow('''
                        SELECT 1 FROM active_voice_mutes
                        WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target='user'
                        AND (expires_at IS NULL OR expires_at > NOW())
                    ''', ctx.guild.id, member.id, channel_obj.id)
                    if not row:
                        skipped_members.append(member)
                        continue
                    if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                        await member.edit(mute=False)
                    await conn.execute('''
                        DELETE FROM active_voice_mutes
                        WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target='user'
                    ''', ctx.guild.id, member.id, channel_obj.id)
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'unmute', member.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Unmuted via room_unmute')
                    unmuted_members.append(member)
                except Exception as e:
                    logger.warning(f'Unmute failed for {member.mention}: {e}')
                    failed_members.append(member)
        summary = f'{self.get_random_emoji()} Unmuted {len(unmuted_members)} member(s) in {channel_obj.mention}.'
        if skipped_members:
            summary += f'\n\U000026A0\U0000FE0F Skipped {len(skipped_members)}.'
        if failed_members:
            summary += f'\n\U0001F6AB Failed to unmute {len(failed_members)}.'
        return await self.handler.send_message(ctx, content=summary, allowed_mentions=discord.AllowedMentions.none())
    
    @app_commands.command(name='xsmute', description='Unmutes a member throughout the entire guild.')
    @is_server_muter_app_predicator()
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    async def undo_server_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
            return await self.send(interaction, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET server_mute_guild_ids = (SELECT ARRAY(
                    SELECT unnest(server_mute_guild_ids)
                    EXCEPT SELECT $2
                ))
                WHERE discord_snowflake = $1
            ''', member_obj.id, guild_id)
            await conn.execute('''
                DELETE FROM active_server_voice_mutes
                WHERE discord_snowflake = $1 AND guild_id = $2
            ''', member_obj.id, guild_id)
        if member_obj.voice and member_obj.voice.channel:
            await member_obj.edit(mute=False)
        await self.send(interaction, content=f'{self.get_random_emoji()} {member_obj.mention} has been server unmuted.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='xsmute', help='Unmutes a member throughout the entire guild.')
    @is_server_muter_predicator()
    async def undo_server_mute_text_command(
            self,
            ctx: commands.Context,
            member: Optional[str] = commands.parameter(description='Tag a member or include their snowflake ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.resolve_member(ctx, member)
        if not member_obj:
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
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} {member_obj.mention} has been server unmuted.', allowed_mentions=discord.AllowedMentions.none())
    
    @app_commands.command(name='xstage', description='Destroy the stage in the current channel.')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_coordinator_app_predicator()
    async def stage_quit_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.send(interaction, content=f'\U0001F6AB {channel_obj.mention} is not a valid target.')
        temp_rooms = self.temp_rooms.get(interaction.guild.id, {})
        room_name = ''
        is_temp_room = False
        for temp_channel in temp_rooms.values():
            if temp_channel.id == getattr(channel_obj, 'id', None) or temp_channel.name == getattr(channel_obj, 'name', None):
                room_name = temp_channel.room_name
                is_temp_room = True
                channel_obj = temp_channel
                break
        async with self.bot.db_pool.acquire() as conn:
            stage = await conn.fetchrow('SELECT initiator_id FROM active_stages WHERE guild_id=$1 AND channel_id=$2', interaction.guild.id, channel_obj.id)
            if not stage:
                return await self.send(interaction, content='\U000026A0\U0000FE0F No active stage found.')
            is_owner_or_dev, _ = await check_owner_dev_coord_app(interaction, channel_obj)
            if interaction.user.id != stage['initiator_id'] and not is_owner_or_dev:
                return await self.send(interaction, content='\U0001F6AB Only the stage initiator can end this stage.')
            await conn.execute('''
                DELETE FROM active_stages
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', interaction.guild.id, channel_obj.id, room_name)
            await conn.execute('''
                DELETE FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', interaction.guild.id, channel_obj.id, room_name)
            await conn.execute('DELETE FROM active_voice_mutes WHERE guild_id=$1 AND channel_id=$2 AND target=\'room\'', interaction.guild.id, channel_obj.id)
            for member in channel_obj.members:
                user_mute = await conn.fetchrow('''
                    SELECT 1 FROM active_voice_mutes
                    WHERE guild_id=$1 AND channel_id=$2 AND discord_snowflake=$3 AND target='room' AND room_name=$4
                ''', interaction.guild.id, channel_obj.id, member.id, room_name)
                if not user_mute and member.voice and member.voice.mute:
                    try:
                        await member.edit(mute=False, reason='Stage ended — no user-specific mute found')
                    except discord.Forbidden:
                        logger.warning(f'Failed to mute member: {member.mention}.')
        return await self.send(interaction, content=f'{self.get_random_emoji()} Stage in {channel_obj.mention} has ended.')
        
    @commands.command(name='xstage', help='Destroy the stage in the current channel.')
    @is_owner_developer_coordinator_predicator()
    async def stage_quit_text_command(
        self,
        ctx: commands.Context,
        *,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB {channel_obj.mention} is not a valid target.')
        temp_rooms = self.temp_rooms.get(ctx.guild.id, {})
        room_name = ''
        is_temp_room = False
        for temp_channel in temp_rooms.values():
            if temp_channel.name == getattr(channel_obj, 'name', None):
                room_name = temp_channel.room_name
                is_temp_room = True
                channel_obj = temp_channel
                break
        async with self.bot.db_pool.acquire() as conn:
            stage = await conn.fetchrow('SELECT initiator_id FROM active_stages WHERE guild_id=$1 AND channel_id=$2', ctx.guild.id, channel_obj.id)
            if not stage:
                return await self.handler.send_message(ctx, content='\U000026A0\U0000FE0F No active stage found.')
            is_owner_or_dev, _ = await check_owner_dev_coord(ctx, channel_obj)
            if ctx.author.id != stage['initiator_id'] and not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F6AB Only the stage initiator can end this stage.')
            await conn.execute('''
                DELETE FROM active_stages
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', ctx.guild.id, channel_obj.id, room_name)
            await conn.execute('''
                DELETE FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', ctx.guild.id, channel_obj.id, room_name)
            await conn.execute('DELETE FROM active_voice_mutes WHERE guild_id=$1 AND channel_id=$2 AND target=\'room\'', ctx.guild.id, channel_obj.id)
            for member in channel_obj.members:
                user_mute = await conn.fetchrow('''
                    SELECT 1 FROM active_voice_mutes
                    WHERE guild_id=$1 AND channel_id=$2 AND discord_snowflake=$3 AND target='room' AND room_name=$4
                ''', ctx.guild.id, channel_obj.id, member.id, room_name)
                if not user_mute and member.voice and member.voice.mute:
                    try:
                        await member.edit(mute=False, reason='Stage ended — no user-specific mute found')
                    except discord.Forbidden:
                        logger.warning(f'Failed to unmute {member.mention}')
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Stage in {channel_obj.mention} has ended.')

    @app_commands.command(name='xtemp', description='Delete a temporary room.')
    @is_owner_developer_app_predicator()
    async def delete_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str]
    ):
        if not interaction.guild:
            return await self.send(interaction, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel_app(interaction, channel)
        rooms = self.temp_rooms.get(interaction.guild.id,{})
        temp_room = rooms.get(channel_obj.name)
        if not temp_room:
            await self.send(interaction, content='\U0001F6AB This channel is not registered as a temporary room.')
            return
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT owner_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_name=$2', interaction.guild.id, channel_obj.name)
            owner_id = row['owner_snowflake'] if row else None
            if owner_id:
                await conn.execute('UPDATE users SET coordinator_room_names=array_remove(coordinator_room_names,$1), coordinator_channel_ids=array_remove(coordinator_channel_ids,$2), moderator_room_names=array_remove(moderator_room_names,$1), moderator_channel_ids=array_remove(moderator_channel_ids,$2), updated_at=NOW() WHERE discord_snowflake=$3', channel_obj.name, channel_obj.id, owner_id)
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake=$1 AND room_name=$2', interaction.guild.id, channel_obj.name)
        del rooms[channel_obj.name]
        await self.send(interaction, content=f'{self.get_random_emoji()} Temporary room {channel_obj.mention} has been removed.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='xtemp', help='Delete a temporary room.')
    @is_owner_developer_predicator()
    async def delete_temp_room_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.resolve_channel(ctx, channel)
        rooms = self.temp_rooms.get(ctx.guild.id,{})
        temp_room = rooms.get(channel_obj.name)
        if not temp_room:
            await self.handler.send_message(ctx, content='\U0001F6AB This channel is not registered as a temporary room.')
            return
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT owner_snowflake FROM temporary_rooms WHERE guild_snowflake=$1 AND room_name=$2', ctx.guild.id, channel_obj.name)
            owner_id = row['owner_snowflake'] if row else None
            if owner_id:
                await conn.execute('UPDATE users SET coordinator_room_names=array_remove(coordinator_room_names,$1), coordinator_channel_ids=array_remove(coordinator_channel_ids,$2), moderator_room_names=array_remove(moderator_room_names,$1), moderator_channel_ids=array_remove(moderator_channel_ids,$2), updated_at=NOW() WHERE discord_snowflake=$3', channel_obj.name, channel_obj.id, owner_id)
            await conn.execute('DELETE FROM temporary_rooms WHERE guild_snowflake=$1 AND room_name=$2', ctx.guild.id, channel_obj.name)
        del rooms[channel_obj.name]
        await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Temporary room {channel_obj.mention} has been removed.', allowed_mentions=discord.AllowedMentions.none())
        
    @app_commands.command(name='xtrole', description='Removes a team role from the permission table.')
    @app_commands.describe(role='Select a role or provide its numeric ID')
    @is_owner_app_predicator()
    async def undo_team_role_app_command(
        self,
        interaction: discord.Interaction,
        role: Optional[discord.Role] = None
    ):
        if not interaction.guild:
            return await self.send(interaction, content='This command must be used in a server.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE role_permissions
                SET is_team_member = FALSE
                WHERE role_id = $1;
            ''', role.id)
        await self.send(interaction, content=f'{self.get_random_emoji()} Role `{role.name}` (`{role.id}`) is no longer authorized for logging commands.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='xtrole', help='Removes a team role from the permission table.')
    @is_owner_predicator()
    async def undo_team_role_text_command(
        self,
        ctx: commands.Context,
        role_ref: Optional[str] = commands.parameter(description='Tag a role or provide its ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        if not role_ref:
            return await self.handler.send_message(ctx, content='\U0001F6AB You must mention a role or provide its ID.')
        role = None
        if ctx.message.role_mentions:
            role = ctx.message.role_mentions[0]
        else:
            try:
                role_id = int(role_ref)
                role = ctx.guild.get_role(role_id)
            except ValueError:
                pass
        if not role:
            return await self.handler.send_message(ctx, content='\U0001F6AB Could not find that role. Please mention it or use its snowflake ID.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE role_permissions
                SET is_team_member = FALSE
                WHERE role_id = $1;
            ''', role.id)
        await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Role `{role.name}` (`{role.id}`) is no longer authorized for logging commands.')
        
    def create_ban_log_pages(self, ctx: commands.Context, member: discord.Member, channel: discord.VoiceChannel, duration_display: Optional[str], reason: Optional[str], executor: discord.Member, expires_at: Optional[datetime], command_used: Optional[str], was_in_channel: bool = False, is_modification: bool = False, guild: discord.Guild = None, highest_role: Optional[str] = ''):
        if expires_at is None: color, ban_type, duration_emoji = 0xDC143C, '🔒 Permanent', '♾️'
        elif (expires_at - datetime.now(timezone.utc)).days >= 7: color, ban_type, duration_emoji = 0xFF6B35, '⏰ Extended', '📅'
        else: color, ban_type, duration_emoji = 0xFF8C00, '⏱️ Temporary', '⏰'
        title = '🔄 Ban Modified' if is_modification else '🔨 User Banned'
        guild = guild or ctx.guild
        embed_user = discord.Embed(title=f"{title} - User Identity", color=color, timestamp=datetime.now(timezone.utc))
        embed_user.description = f"**Target:** {member.mention} banned from {channel.mention}"
        embed_user.set_thumbnail(url=executor.display_avatar.url)
        user_priority = f"**Display Name:** {member.display_name}\n**Username:** @{member.name}\n**User ID:** `{member.id}`\n**Account Age:** <t:{int(member.created_at.timestamp())}:R>"
        if member.joined_at: user_priority += f"\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
        embed_user.add_field(name='👤 Target User (Priority Info)', value=user_priority, inline=False)
        exec_priority = f"**Moderator:** {executor.display_name} (@{executor.name})\n**Mod ID:** `{executor.id}`\n**Top Role:** {highest_role or executor.top_role.mention}"
        embed_user.add_field(name='👮‍♂️ Executed By', value=exec_priority, inline=True)
        ctx_info = f"**Original Message ID:** `{ctx.message.id}`\n**Message Link:** [Jump to Message]({ctx.message.jump_url})\n**Command Channel:** {ctx.channel.mention}\n**Command Used:** `{command_used}`"
        embed_user.add_field(name='📱 Command Context', value=ctx_info, inline=True)
        if expires_at:
            time_left = expires_at - datetime.now(timezone.utc)
            hours_left = round(time_left.total_seconds() / 3600, 1)
            days_left = time_left.days
            duration_info = f'**Type:** {ban_type}\n**Duration:** {duration_display}\n**Expires:** <t:{int(expires_at.timestamp())}:F>\n**Time Left:** '
            duration_info += f'{days_left}d {hours_left % 24:.1f}h' if days_left > 0 else f'{hours_left}h'
        else:
            duration_info = f'**Type:** {ban_type}\n**Duration:** {duration_display}\n**Status:** Permanent Ban'
        embed_user.add_field(name=f'{duration_emoji} Ban Duration', value=duration_info, inline=False)
        embed_user.add_field(name='📝 Reason', value=f'```{reason if reason else "No reason provided"}```', inline=False)
        embed_user.set_footer(text=f"Ban Ref: {member.id}-{channel.id} | Msg: {ctx.message.id}", icon_url=guild.icon.url if guild and guild.icon else None)
        embed_duration = discord.Embed(title=f"{title} - Duration Info", color=color, timestamp=datetime.now(timezone.utc))
        embed_duration.add_field(name=f'{duration_emoji} Ban Duration', value=duration_info, inline=False)
        action_details = f"**Was in Channel:** {'✅ Yes' if was_in_channel else '\U0001F6AB No'}\n**Action Type:** {'Modification' if is_modification else 'New Ban'}\n**Server:** {guild.name} (`{guild.id}`)"
        embed_duration.add_field(name='⚙️ Action Details', value=action_details, inline=True)
        channel_basic = f"**Channel:** {channel.mention} (`{channel.id}`)\n**Category:** {channel.category.name if channel.category else 'None'}"
        embed_duration.add_field(name='📍 Channel Info', value=channel_basic, inline=True)
        embeds = [embed_user, embed_duration]
        if reason:
            reason_chunks = [reason[i:i+1000] for i in range(0, len(reason), 1000)]
            if len(reason_chunks) > 1:
                for i, chunk in enumerate(reason_chunks):
                    reason_embed = discord.Embed(title=f"{title} - Reason (cont.)", color=color, timestamp=datetime.now(timezone.utc))
                    reason_embed.add_field(name=f'📝 Reason (Part {i+1})', value=f'```{chunk}```', inline=False)
                    embeds.append(reason_embed)
        return embeds
        
    def create_text_mute_log_pages(self, ctx: commands.Context, member: discord.Member, channel: discord.VoiceChannel, duration_display: Optional[str], reason: Optional[str], executor: discord.Member, expires_at: Optional[datetime], command_used: Optional[str], was_in_channel: bool = False, is_modification: bool = False, guild: discord.Guild = None, highest_role: Optional[str] = ''):
        if expires_at is None: color, mute_type, duration_emoji = 0xDC143C, '🔒 Permanent', '♾️'
        elif (expires_at - datetime.now(timezone.utc)).days >= 7: color, mute_type, duration_emoji = 0xFF6B35, '⏰ Extended', '📅'
        else: color, mute_type, duration_emoji = 0xFF8C00, '⏱️ Temporary', '⏰'
        title = '🔄 Text Mute Modified' if is_modification else '🔇 User Text Muted'
        guild = guild or ctx.guild
        embed_user = discord.Embed(title=f"{title} - User Identity", color=color, timestamp=datetime.now(timezone.utc))
        embed_user.description = f"**Target:** {member.mention} text muted in {channel.mention}"
        embed_user.set_thumbnail(url=executor.display_avatar.url)
        user_priority = f"**Display Name:** {member.display_name}\n**Username:** @{member.name}\n**User ID:** `{member.id}`\n**Account Age:** <t:{int(member.created_at.timestamp())}:R>"
        if member.joined_at: user_priority += f"\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
        embed_user.add_field(name='👤 Target User (Priority Info)', value=user_priority, inline=False)
        exec_priority = f"**Moderator:** {executor.display_name} (@{executor.name})\n**Mod ID:** `{executor.id}`\n**Top Role:** {highest_role or executor.top_role.mention}"
        embed_user.add_field(name='👮‍♂️ Executed By', value=exec_priority, inline=True)
        ctx_info = f"**Original Message ID:** `{ctx.message.id}`\n**Message Link:** [Jump to Message]({ctx.message.jump_url})\n**Command Channel:** {ctx.channel.mention}\n**Command Used:** `{command_used}`"
        embed_user.add_field(name='📱 Command Context', value=ctx_info, inline=True)
        if expires_at:
            time_left = expires_at - datetime.now(timezone.utc)
            hours_left = round(time_left.total_seconds() / 3600, 1)
            days_left = time_left.days
            duration_info = f'**Type:** {mute_type}\n**Duration:** {duration_display}\n**Expires:** <t:{int(expires_at.timestamp())}:F>\n**Time Left:** '
            duration_info += f'{days_left}d {hours_left % 24:.1f}h' if days_left > 0 else f'{hours_left}h'
        else:
            duration_info = f'**Type:** {mute_type}\n**Duration:** {duration_display}\n**Status:** Permanent Text Mute'
        embed_user.add_field(name=f'{duration_emoji} Text Mute Duration', value=duration_info, inline=False)
        embed_user.add_field(name='📝 Reason', value=f'```{reason if reason else "No reason provided"}```', inline=False)
        embed_user.set_footer(text=f"Text Mute Ref: {member.id}-{channel.id} | Msg: {ctx.message.id}", icon_url=guild.icon.url if guild and guild.icon else None)
        embed_duration = discord.Embed(title=f"{title} - Duration Info", color=color, timestamp=datetime.now(timezone.utc))
        embed_duration.add_field(name=f'{duration_emoji} Text Mute Duration', value=duration_info, inline=False)
        action_details = f"**Was in Channel:** {'✅ Yes' if was_in_channel else '\U0001F6AB No'}\n**Action Type:** {'Modification' if is_modification else 'New Text Mute'}\n**Server:** {guild.name} (`{guild.id}`)"
        embed_duration.add_field(name='⚙️ Action Details', value=action_details, inline=True)
        channel_basic = f"**Channel:** {channel.mention} (`{channel.id}`)\n**Category:** {channel.category.name if channel.category else 'None'}"
        embed_duration.add_field(name='📍 Channel Info', value=channel_basic, inline=True)
        embeds = [embed_user, embed_duration]
        if reason:
            reason_chunks = [reason[i:i+1000] for i in range(0, len(reason), 1000)]
            if len(reason_chunks) > 1:
                for i, chunk in enumerate(reason_chunks):
                    reason_embed = discord.Embed(title=f"{title} - Reason (cont.)", color=color, timestamp=datetime.now(timezone.utc))
                    reason_embed.add_field(name=f'📝 Reason (Part {i+1})', value=f'```{chunk}```', inline=False)
                    embeds.append(reason_embed)
        return embeds
    
    def create_voice_mute_log_pages(self, ctx: commands.Context, member: discord.Member, channel: discord.VoiceChannel, duration_display: Optional[str], reason: Optional[str], executor: discord.Member, expires_at: Optional[datetime], command_used: Optional[str], was_in_channel: bool = False, is_modification: bool = False, guild: discord.Guild = None, highest_role: Optional[str] = ''):
        if expires_at is None:
            color, mute_type, duration_emoji = 0xDC143C, '🔒 Permanent', '♾️'
        elif (expires_at - datetime.now(timezone.utc)).days >= 7:
            color, mute_type, duration_emoji = 0xFF6B35, '⏰ Extended', '📅'
        else:
            color, mute_type, duration_emoji = 0xFF8C00, '⏱️ Temporary', '⏰'
        title = '🔄 Voice Mute Modified' if is_modification else '🎙️ User Voice Muted'
        guild = guild or ctx.guild
        embed_user = discord.Embed(title=f"{title} - User Identity", color=color, timestamp=datetime.now(timezone.utc))
        embed_user.description = f"**Target:** {member.mention} voice muted in {channel.mention}"
        embed_user.set_thumbnail(url=executor.display_avatar.url)
        user_priority = f"**Display Name:** {member.display_name}\n**Username:** @{member.name}\n**User ID:** `{member.id}`\n**Account Age:** <t:{int(member.created_at.timestamp())}:R>"
        if member.joined_at:
            user_priority += f"\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
        embed_user.add_field(name='👤 Target User (Priority Info)', value=user_priority, inline=False)
        exec_priority = f"**Moderator:** {executor.display_name} (@{executor.name})\n**Mod ID:** `{executor.id}`\n**Top Role:** {highest_role or executor.top_role.mention}"
        embed_user.add_field(name='👮‍♂️ Executed By', value=exec_priority, inline=True)
        ctx_info = f"**Original Message ID:** `{ctx.message.id}`\n**Message Link:** [Jump to Message]({ctx.message.jump_url})\n**Command Channel:** {ctx.channel.mention}\n**Command Used:** `{command_used}`"
        embed_user.add_field(name='📱 Command Context', value=ctx_info, inline=True)
        if expires_at:
            time_left = expires_at - datetime.now(timezone.utc)
            hours_left, days_left = round(time_left.total_seconds() / 3600, 1), time_left.days
            duration_info = f'**Type:** {mute_type}\n**Duration:** {duration_display}\n**Expires:** <t:{int(expires_at.timestamp())}:F>\n**Time Left:** '
            duration_info += f'{days_left}d {hours_left % 24:.1f}h' if days_left > 0 else f'{hours_left}h'
        else:
            duration_info = f'**Type:** {mute_type}\n**Duration:** {duration_display}\n**Status:** Permanent Voice Mute'
        embed_user.add_field(name=f'{duration_emoji} Voice Mute Duration', value=duration_info, inline=False)
        embed_user.add_field(name='📝 Reason', value=f'```{reason if reason else "No reason provided"}```', inline=False)
        embed_user.set_footer(text=f"Voice Mute Ref: {member.id}-{channel.id} | Msg: {ctx.message.id}", icon_url=guild.icon.url if guild and guild.icon else None)
        embed_duration = discord.Embed(title=f"{title} - Duration Info", color=color, timestamp=datetime.now(timezone.utc))
        embed_duration.add_field(name=f'{duration_emoji} Voice Mute Duration', value=duration_info, inline=False)
        action_details = f"**Was in Channel:** {'✅ Yes' if was_in_channel else '\U0001F6AB No'}\n**Action Type:** {'Modification' if is_modification else 'New Voice Mute'}\n**Server:** {guild.name} (`{guild.id}`)"
        embed_duration.add_field(name='⚙️ Action Details', value=action_details, inline=True)
        channel_basic = f"**Channel:** {channel.mention} (`{channel.id}`)\n**Category:** {channel.category.name if channel.category else 'None'}"
        embed_duration.add_field(name='📍 Channel Info', value=channel_basic, inline=True)
        embeds = [embed_user, embed_duration]
        if reason:
            reason_chunks = [reason[i:i+1000] for i in range(0, len(reason), 1000)]
            if len(reason_chunks) > 1:
                for i, chunk in enumerate(reason_chunks):
                    reason_embed = discord.Embed(title=f"{title} - Reason (cont.)", color=color, timestamp=datetime.now(timezone.utc))
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
        if days > 0:
            return f"{days} day(s)"
        elif hours > 0:
            return f"{hours} hour(s)"
        elif minutes > 0:
            return f"{minutes} minute(s)"
        else:
            return f"less than a minute"
    
    async def get_caps_for_channel(self, guild_id: int, channel_id: int) -> list[tuple[str, str]]:
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
    
    async def resolve_channel(self, ctx: commands.Context, value: Optional[Union[int, str, discord.TextChannel, discord.VoiceChannel]]) -> Optional[Union[discord.TextChannel, discord.VoiceChannel]]:
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

    async def resolve_member_app(self, interaction: discord.Interaction, value: Optional[Union[int, str, discord.Member]]) -> Optional[discord.Member]:
        try:
            if isinstance(value, discord.Member):
                logger.debug(f"Direct member: {value.id}")
                return value
            if isinstance(value, int):
                m = interaction.guild.get_member(value)
                if not m:
                    try: m = await interaction.guild.fetch_member(value)
                    except discord.NotFound: m = None
                if m:
                    logger.debug(f"Resolved member by int ID: {m.id}")
                    return m
            if isinstance(value, str):
                if value.isdigit():
                    mid = int(value)
                    m = interaction.guild.get_member(mid)
                    if not m:
                        try: m = await interaction.guild.fetch_member(mid)
                        except discord.NotFound: m = None
                    if m:
                        logger.debug(f"Resolved member by str ID: {m.id}")
                        return m
                if value.startswith('<@') and value.endswith('>'):
                    mid = int(value[2:-1].replace('!', ''))
                    m = interaction.guild.get_member(mid)
                    if not m:
                        try: m = await interaction.guild.fetch_member(mid)
                        except discord.NotFound: m = None
                    if m:
                        logger.debug(f"Member mention resolved: {m.id}")
                        return m
                    return m
        except Exception as e:
            logger.warning(f"Member resolution error: {e}")
        return None
    
    async def resolve_channel_app(self, interaction: discord.Interaction, value: Optional[Union[int, str, discord.TextChannel, discord.VoiceChannel]]) -> Optional[Union[discord.TextChannel, discord.VoiceChannel]]:
        try:
            if isinstance(value, (discord.TextChannel, discord.VoiceChannel)):
                logger.debug(f"Direct channel: {value.id}")
                return value
            if isinstance(value, int):
                c = interaction.guild.get_channel(value)
                if isinstance(c, (discord.TextChannel, discord.VoiceChannel)):
                    logger.debug(f"Resolved channel by int ID: {c.id}")
                    return c
            if isinstance(value, str):
                if value.isdigit():
                    cid = int(value)
                    c = interaction.guild.get_channel(cid)
                    if isinstance(c, (discord.TextChannel, discord.VoiceChannel)):
                        logger.debug(f"Resolved channel by str ID: {c.id}")
                        return c
                if value.startswith('<#') and value.endswith('>'):
                    cid = int(value[2:-1])
                    c = interaction.guild.get_channel(cid)
                    if c:
                        logger.debug(f"Channel mention resolved: {c.id}")
                        return c
        except Exception as e:
            logger.warning(f"Channel resolution error: {e}")
        return interaction.channel

    def get_random_emoji(self):
        return random.choice(VEGAN_EMOJIS)
                
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
            return target, f'{value} day(s)' if sign > 0 else f'{value} day(s)'
        if duration.endswith(('h', 'hr', 'hrs', 'hour', 'hours')):
            value = int(duration.rstrip('hrshours'))
            delta = timedelta(hours=value * sign)
            target = (base if (is_relative and base) else datetime.now(timezone.utc)) + delta
            return target, f'{value} hour(s)' if sign > 0 else f'{value} hour(s)'
        if duration.endswith(('m', 'min', 'mins', 'minute', 'minutes')):
            value = int(duration.rstrip('minsmutes'))
            delta = timedelta(minutes=value * sign)
            target = (base if (is_relative and base) else datetime.now(timezone.utc)) + delta
            return target, f'{value} minute(s)' if sign > 0 else f'{value} minute(s)'
        value = int(duration)
        delta = timedelta(hours=value * sign)
        target = (base if (is_relative and base) else datetime.now(timezone.utc)) + delta
        return target, f'{value} hour(s)' if sign > 0 else f'{value} hour(s)'

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
    
    async def send(self, interaction, **kw):
        await interaction.response.send_message(**kw, ephemeral=True)
            
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

    def setup_backup_directory(self, backup_dir: Optional[str]) -> str:
        os.makedirs(backup_dir, exist_ok=True)
        return backup_dir
        
async def setup(bot: DiscordBot):
    cog = Hybrid(bot)
    await bot.add_cog(cog)
    cog.list_bans_app_command.callback._team_command = True
    cog.list_mutes_text_command.callback._team_command = True
    cog.list_mutes_app_command.callback._team_command = True
    cog.list_bans_text_command.callback._team_command = True
    cog.list_text_mutes_text_command.callback._team_command = True
    cog.list_text_mutes_app_command.callback._team_command = True

class TempChannel:
    def __init__(self, channel: discord.abc.GuildChannel, room_name: str):
        self.channel = channel
        self.is_temp_room = True
        self.room_name = room_name

    def __getattr__(self, name):
        return getattr(self.channel, name)
