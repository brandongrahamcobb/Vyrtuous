''' event_listeners.py

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
import discord
from datetime import datetime, timezone
import os
import time
from discord.ext import commands
from vyrtuous.inc.helpers import *
from vyrtuous.utils.setup_logging import logger
from vyrtuous.service.check_service import *
from vyrtuous.service.discord_message_service import DiscordMessageService, ChannelPaginator
from vyrtuous.bot.discord_bot import DiscordBot
import asyncio
from typing import defaultdict

class EventListeners(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.db_pool)
        self.join_log = defaultdict(list)
        self._ready_done = False
    
    async def fetch_active_bans(self, conn, guild_id: int, user_id: int, temp_room_name: str = None):
        if temp_room_name:
            return await conn.fetch('''
                SELECT channel_id, expires_at
                FROM active_bans
                WHERE guild_id = $1
                  AND discord_snowflake = $2
                  AND room_name = $3
                  AND (expires_at IS NULL OR expires_at > NOW())
            ''', guild_id, user_id, temp_room_name)
        else:
            return await conn.fetch('''
                SELECT channel_id, expires_at
                FROM active_bans
                WHERE guild_id = $1
                  AND discord_snowflake = $2
                  AND (expires_at IS NULL OR expires_at > NOW())
            ''', guild_id, user_id)
    
    async def fetch_active_text_mutes(self, conn, guild_id: int, user_id: int, temp_room_name: str = None):
        if temp_room_name:
            return await conn.fetch('''
                SELECT channel_id
                FROM active_text_mutes
                WHERE guild_id = $1
                  AND discord_snowflake = $2
                  AND room_name = $3
                  AND (expires_at IS NULL OR expires_at > NOW())
            ''', guild_id, user_id, temp_room_name)
        else:
            return await conn.fetch('''
                SELECT channel_id
                FROM active_text_mutes
                WHERE guild_id = $1
                  AND discord_snowflake = $2
                  AND (expires_at IS NULL OR expires_at > NOW())
            ''', guild_id, user_id)
    
    @commands.Cog.listener()
    async def on_guild_channel_create(self, channel: discord.abc.GuildChannel):
        guild = channel.guild
        name = channel.name.lower()
        for c in guild.channels:
            if c.id != channel.id and c.name.lower() == name:
                return
        async with self.bot.db_pool.acquire() as conn:
            temp_room = await conn.fetchrow(
                '''
                SELECT room_name
                FROM temporary_rooms
                WHERE guild_snowflake=$1
                  AND room_name=$2
                ''',
                guild.id, name
            )
            if not temp_room:
                print(f"DEBUG: Not a temporary room, exiting")
                return
            await conn.execute(
                '''
                UPDATE temporary_rooms
                SET room_snowflake=$3
                WHERE guild_snowflake=$1
                  AND room_name=$2
                ''',
                guild.id, name, channel.id
            )
            await conn.execute(
                '''
                UPDATE command_aliases
                SET channel_id=$3
                WHERE guild_id=$1
                  AND room_name=$2
                ''',
                guild.id, name, channel.id
            )
        guild_aliases = self.bot.command_aliases.setdefault(guild.id, {})
        temp_aliases = guild_aliases.get('temp_room_aliases', {})
        for alias_type, aliases in temp_aliases.items():
            for alias_name, data in aliases.items():
                if data.get("room_name") == name:
                    data["channel_id"] = channel.id
                    
    # Done
    @commands.Cog.listener()
    async def on_voice_state_update(self, member: discord.Member, before: discord.VoiceState, after: discord.VoiceState) -> None:
        if before.channel == after.channel and before.mute == after.mute and before.self_mute == after.self_mute:
            return
        if member.bot:
            return
        after_channel = after.channel
        before_channel = before.channel
        user_id = member.id
        should_be_muted = False
        just_manual_unmute = False
        async with self.db_pool.acquire() as conn:
            user_data = await conn.fetchrow('SELECT server_mute_guild_ids FROM users WHERE discord_snowflake = $1', user_id)
            server_mute_guild_ids = user_data['server_mute_guild_ids'] or [] if user_data else []
            if member.guild.id in server_mute_guild_ids:
                return
            active_stage = None
            coordinator_ids = set()
            if after_channel:
                temp_room = await conn.fetchrow('''
                    SELECT room_name FROM temporary_rooms
                    WHERE guild_snowflake = $1 AND room_snowflake = $2
                ''', member.guild.id, after_channel.id)
                room_name = temp_room['room_name'] if temp_room else ''
                try:
                    if temp_room:
                        room_name = temp_room['room_name']
                        active_stage = await conn.fetchrow('''
                            SELECT expires_at FROM active_stages
                            WHERE guild_id = $1 AND channel_id = $2 AND room_name = $3
                        ''', member.guild.id, after_channel.id, room_name)
                        coordinators = await conn.fetch('''
                            SELECT discord_snowflake FROM stage_coordinators
                            WHERE guild_id = $1 AND channel_id = $2 AND room_name = $3
                        ''', member.guild.id, after_channel.id, room_name)
                    else:
                        active_stage = await conn.fetchrow('''
                            SELECT expires_at FROM active_stages
                            WHERE guild_id = $1 AND channel_id = $2
                        ''', member.guild.id, after_channel.id)
                        coordinators = await conn.fetch('''
                            SELECT discord_snowflake FROM stage_coordinators
                            WHERE guild_id = $1 AND channel_id = $2
                        ''', member.guild.id, after_channel.id)
                    if after_channel is not None and after_channel != before.channel:
                        if active_stage:
                            now = time.time()
                            self.join_log[member.id] = [t for t in self.join_log[member.id] if now - t < 300]
                            if len(self.join_log[member.id]) < 1:
                                self.join_log[member.id].append(now)
                                expires = active_stage['expires_at']
                                embed = discord.Embed(
                                    title=f"\U0001F399 {after_channel.name} — Stage Mode",
                                    description=f"Ends <t:{int(expires.timestamp())}:R>",
                                    color=discord.Color.green()
                                )
                                embed.add_field(name="\u200b", value="**Ask to speak!**", inline=False)
                                await after_channel.send(embed=embed)
                        rows = await conn.fetch('''
                            SELECT channel_id, discord_snowflake, reason
                            FROM active_flags
                            WHERE guild_id = $1 AND discord_snowflake = $2
                        ''', member.guild.id, user_id)
                        if rows:
                            grouped = {}
                            for row in rows:
                                grouped.setdefault(row['channel_id'], []).append(row)
                            context_records = grouped.get(after_channel.id)
                            if after_channel.id == 1222056499959042108 and context_records:
                                if context_records and after_channel.id:
                                    embeds = []
                                    embed = discord.Embed(
                                        title=f'\u26A0\uFE0F {member.display_name} is flagged',
                                        color=discord.Color.red()
                                    )
                                    embed.set_thumbnail(url=member.display_avatar.url)
                                    for record in context_records:
                                        reason = record['reason'] or 'No reason provided'
                                        embed.add_field(name=f'Channel: {after_channel.mention}', value=f'Reason: {reason}', inline=False)
                                    other_channels = [ch_id for ch_id in grouped.keys() if ch_id != after_channel.id]
                                    if other_channels:
                                        ch_mentions = []
                                        for ch_id in other_channels:
                                            ch = member.guild.get_channel(ch_id)
                                            if not ch:
                                                ch = await member.guild.fetch_channel(ch_id)
                                            ch_mentions.append(ch.mention if ch else f'Channel ID `{ch_id}`')
                                        embed.add_field(name='Other flagged channels', value='\n'.join(ch_mentions), inline=False)
                                    embeds.append(embed)
                                    for ch_id in other_channels:
                                        records = grouped[ch_id]
                                        ch = member.guild.get_channel(ch_id)
                                        ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                                        embed = discord.Embed(
                                            title=f'\u26A0\uFE0F {member.display_name} is flagged in {ch_name}',
                                            color=discord.Color.red()
                                        )
                                        embed.set_thumbnail(url=member.display_avatar.url)
                                        for record in records:
                                            reason = record['reason'] or 'No reason provided'
                                            embed.add_field(name='Channel', value=f'{ch_name}\nReason: {reason}', inline=False)
                                        embeds.append(embed)
                                    now = time.time()
                                    self.join_log[member.id] = [t for t in self.join_log[member.id] if now - t < 300]
                                    if len(self.join_log[member.id]) < 1:
                                        self.join_log[member.id].append(now)
                                        if len(embeds) == 1:
                                            await after_channel.send(embed=embeds[0])
                                        else:
                                            paginator = ChannelPaginator(self.bot, after_channel, embeds)
                                            await paginator.start()
                    coordinator_ids = {r['discord_snowflake'] for r in coordinators}
                    is_owner_or_dev = await is_owner_developer_via_objects(member, self.bot)
                    if before.mute != after.mute:
                        if before.mute and not after.mute and before_channel:
                            if temp_room:
                                room_name = temp_room['room_name']
                                result = await conn.execute('''
                                    DELETE FROM active_voice_mutes
                                    WHERE guild_id = $1
                                      AND discord_snowflake = $2
                                      AND channel_id = $3
                                      AND room_name = $4
                                      AND target = 'user'
                                      AND expires_at IS NOT NULL
                                ''', member.guild.id, user_id, before_channel.id, room_name)
                            else:
                                result = await conn.execute('''
                                    DELETE FROM active_voice_mutes
                                    WHERE guild_id = $1
                                      AND discord_snowflake = $2
                                      AND channel_id = $3
                                      AND target = 'user'
                                      AND expires_at IS NOT NULL
                                ''', member.guild.id, user_id, before_channel.id)
                            just_manual_unmute = True
                        elif active_stage and before.mute and not after.mute and before_channel:
                            try:
                                if temp_room:
                                    room_name = temp_room['room_name']
                                    await conn.execute('''
                                        DELETE FROM active_voice_mutes
                                        WHERE guild_id = $1
                                          AND discord_snowflake = $2
                                          AND channel_id = $3
                                          AND room_name = $4
                                          AND target = 'room'
                                    ''', member.guild.id, user_id, after_channel.id, room_name)
                                else:
                                    await conn.execute('''
                                        DELETE FROM active_voice_mutes
                                        WHERE guild_id = $1
                                          AND discord_snowflake = $2
                                          AND channel_id = $3
                                          AND target = 'room'
                                    ''', member.guild.id, user_id, after_channel.id)
                                just_manual_unmute = True
                                await member.edit(mute=False, reason=f'Removing stage mute in {before_channel.name}')
                                return
                            except discord.Forbidden:
                                logger.debug(f'No permission to unmute {member.display_name}')
                                raise
                            except discord.HTTPException as e:
                                logger.debug(f'Failed to unmute {member.display_name}: {e}')
                                raise
                    if active_stage and (user_id not in coordinator_ids and not is_owner_or_dev) and not after.mute and before_channel != after_channel and after_channel:
                         expires_at = active_stage['expires_at']
                         try:
                             if temp_room:
                                 room_name = temp_room['room_name']
                                 await conn.execute('''
                                     INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, target, room_name)
                                     VALUES ($1, $2, $3, $4, 'room', $5)
                                     ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
                                     DO UPDATE SET expires_at = EXCLUDED.expires_at
                                 ''', member.guild.id, user_id, after_channel.id, expires_at, room_name)
                             else:
                                 await conn.execute('''
                                     INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, target)
                                     VALUES ($1, $2, $3, $4, 'room')
                                     ON CONFLICT (guild_id, discord_snowflake, channel_id, target)
                                     DO UPDATE SET expires_at = EXCLUDED.expires_at
                                 ''', member.guild.id, user_id, after_channel.id, expires_at)
                             await member.edit(mute=True, reason=f'Enforcing stage mute in {after_channel.name}')
                             return
                         except discord.Forbidden:
                             logger.debug(f'No permission to mute {member.display_name}')
                         except discord.HTTPException as e:
                             logger.debug(f'Failed to mute {member.display_name}: {e}')
                    existing_mute_row = await conn.fetchrow('''
                        SELECT expires_at
                        FROM active_voice_mutes
                        WHERE guild_id = $1
                          AND discord_snowflake = $2
                          AND channel_id = $3
                          AND target = 'user'
                    ''', member.guild.id, user_id, after_channel.id)
                    if existing_mute_row:
                        if existing_mute_row['expires_at'] and existing_mute_row['expires_at'] <= datetime.now(timezone.utc):
                            if temp_room:
                                room_name = temp_room['room_name']
                                existing_mute_row = await conn.fetchrow('''
                                    SELECT expires_at
                                    FROM active_voice_mutes
                                    WHERE guild_id = $1
                                      AND discord_snowflake = $2
                                      AND channel_id = $3
                                      AND room_name = $4
                                      AND target = 'user'
                                ''', member.guild.id, user_id, after_channel.id, room_name)
                            else:
                                existing_mute_row = await conn.fetchrow('''
                                    SELECT expires_at
                                    FROM active_voice_mutes
                                    WHERE guild_id = $1
                                      AND discord_snowflake = $2
                                      AND channel_id = $3
                                      AND target = 'user'
                                ''', member.guild.id, user_id, after_channel.id)
                            existing_mute_row = None
                            should_be_muted = False
                    if not before.mute and after.mute and after_channel:
                        if not existing_mute_row:
                            other_cog = self.bot.get_cog("Hybrid")
                            if other_cog is not None:
                                if member.id not in other_cog.super["members"]:
                                    if temp_room:
                                        room_name = temp_room['room_name']
                                        await conn.execute('''
                                            INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, target, room_name)
                                            VALUES ($1, $2, $3, NOW() + interval '1 hour', 'user', $4)
                                            ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target) DO UPDATE
                                            SET expires_at = EXCLUDED.expires_at
                                        ''', member.guild.id, user_id, after_channel.id, room_name)
                                    else:
                                        await conn.execute('''
                                            INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, target)
                                            VALUES ($1, $2, $3, NOW() + interval '1 hour', 'user')
                                            ON CONFLICT (guild_id, discord_snowflake, channel_id, target) DO UPDATE
                                            SET expires_at = EXCLUDED.expires_at
                                        ''', member.guild.id, user_id, after_channel.id)
                                else:
                                    embed = discord.Embed(
                                        title=f'\u1F4AB {member.display_name} is a hero!',
                                        description=f'{member.display_name} cannot be muted.',
                                        color=discord.Color.gold()
                                    )
                                    embed.set_thumbnail(url=member.display_avatar.url)
                                    await after_channel.send(embed=embed)
                    existing_mute_row = await conn.fetchrow('''
                        SELECT expires_at
                        FROM active_voice_mutes
                        WHERE guild_id = $1
                          AND discord_snowflake = $2
                          AND channel_id = $3
                          AND target = 'user'
                    ''', member.guild.id, user_id, after_channel.id)
                    should_be_muted = False
                    if existing_mute_row:
                        if not existing_mute_row['expires_at'] or existing_mute_row['expires_at'] > datetime.now(timezone.utc):
                            should_be_muted = True
                    if just_manual_unmute and existing_mute_row:
                        if temp_room:
                            room_name = temp_room['room_name']
                            existing_mute_row = await conn.fetchrow('''
                                SELECT expires_at
                                FROM active_voice_mutes
                                WHERE guild_id = $1
                                  AND discord_snowflake = $2
                                  AND channel_id = $3
                                  AND room_name = $4
                                  AND target = 'user'
                            ''', member.guild.id, user_id, after_channel.id, room_name)
                        else:
                            existing_mute_row = await conn.fetchrow('''
                                SELECT expires_at
                                FROM active_voice_mutes
                                WHERE guild_id = $1
                                  AND discord_snowflake = $2
                                  AND channel_id = $3
                                  AND target = 'user'
                            ''', member.guild.id, user_id, after_channel.id)
                        records = [r for r in records if member.guild.get_channel(r['channel_id'])]
                        if not records:
                            return
                        description_lines = []
                        perm_lines = []
                        for record in records:
                            channel_obj = member.guild.get_channel(record['channel_id'])
                            channel_mention = channel_obj.mention if channel_obj else f'`{record['channel_id']}`'
                            reason = record['reason']
                            if record['expires_at'] is None:
                                perm_lines.append(f'• {channel_mention} — {reason}')
                        if perm_lines:
                            perm_embed = discord.Embed(title=f'\U0001F507 Attempted unmute: {member.display_name}', description='\n'.join(perm_lines)+'\n\nThis user must be unmuted manually via the bot.',     color=discord.Color.red())
                            await after_channel.send(embed=perm_embed)
                    if should_be_muted and not after.mute:
                        try:
                            await member.edit(mute=True, reason=f'Enforcing mute in {after_channel.name} (found in arrays)')
                        except discord.Forbidden:
                            logger.debug(f'No permission to mute {member.display_name}')
                        except discord.HTTPException as e:
                            logger.debug(f'Failed to mute {member.display_name}: {e}')
                    elif not should_be_muted and after.mute:
                        try:
                            await member.edit(mute=False, reason=f'Auto-unmuting in {after_channel.name} (no mute record)')
                        except discord.Forbidden:
                            logger.debug(f'No permission to unmute {member.display_name}')
                        except discord.HTTPException as e:
                            logger.debug(f'Failed to unmute {member.display_name}: {e}')
                except Exception as e:
                    logger.exception('Error in on_voice_state_update', exc_info=e)
            
    @commands.Cog.listener()
    async def on_member_join(self, member: discord.Member) -> None:
        user_id = member.id
        guild = member.guild
        async with self.db_pool.acquire() as conn:
            temp_room_name = None
            temp_room = await conn.fetchrow('''
                SELECT room_name FROM temporary_rooms
                WHERE guild_snowflake = $1
                  AND owner_snowflake = $2
            ''', guild.id, user_id)
            if temp_room:
                temp_room_name = temp_room['room_name']
            bans = await self.fetch_active_bans(conn, guild.id, user_id, temp_room_name)
            text_mutes = await self.fetch_active_text_mutes(conn, guild.id, user_id, temp_room_name)
            for row in bans:
                channel = guild.get_channel(row['channel_id'])
                if not channel or not isinstance(channel, (discord.TextChannel, discord.VoiceChannel)):
                    continue
                if row['expires_at'] and row['expires_at'] < datetime.now(timezone.utc):
                    continue
                try:
                    overwrite = channel.overwrites_for(member)
                    overwrite.view_channel = False
                    await channel.set_permissions(member, overwrite=overwrite, reason='Reinstating active channel ban')
                except discord.Forbidden:
                    print(f'Missing permissions to ban in channel {channel.id}')
                except discord.HTTPException as e:
                    print(f'Failed to apply ban for {member} in {channel.id}: {e}')
            for row in text_mutes:
                channel = guild.get_channel(row['channel_id'])
                if not channel or not isinstance(channel, discord.TextChannel):
                    continue
                try:
                    overwrite = channel.overwrites_for(member)
                    overwrite.send_messages = False
                    await channel.set_permissions(member, overwrite=overwrite, reason='Reinstating text mute')
                except discord.Forbidden:
                    print(f'Missing permissions to text mute in channel {channel.id}')
                except discord.HTTPException as e:
                    print(f'Failed to apply text mute for {member} in {channel.id}: {e}')

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if after.author.bot:
            return
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)
                
    @commands.Cog.listener()
    async def on_command_error(self, ctx, error):
        if isinstance(error, commands.MissingRequiredArgument):
            missing = error.param.name
            await ctx.reply(f'\U0001F6AB Missing required argument: `{missing}`')
            return
        if isinstance(error, commands.CheckFailure):
            return await send_check_failure_embed(ctx, error)
            
#    @commands.Cog.listener()
#    async def on_command(self, ctx):
#        await ctx.send("Bot is currently down. Changes will not be saved permanently.")

    @commands.Cog.listener()
    async def on_ready(self):
        if getattr(self, "_ready_done", False):
            return
        self._ready_done = True
        hybrid_cog = self.bot.get_cog("Hybrid")
        if hybrid_cog:
            await hybrid_cog.load_temp_rooms()
        async with self.bot.db_pool.acquire() as conn:
            bans = await conn.fetch('SELECT discord_snowflake, channel_id, room_name FROM active_bans')
            texts = await conn.fetch('SELECT discord_snowflake, channel_id, room_name FROM active_text_mutes')
            voices = await conn.fetch('SELECT discord_snowflake, channel_id, room_name FROM active_voice_mutes')
            ban_set = {(r['discord_snowflake'], r['channel_id'], r['room_name'] or '') for r in bans}
            text_set = {(r['discord_snowflake'], r['channel_id'], r['room_name'] or '') for r in texts}
            voice_set = {(r['discord_snowflake'], r['channel_id'], r['room_name'] or '') for r in voices}
            for guild in self.bot.guilds:
                for channel in guild.channels:
                    temp_room = await conn.fetchrow('''
                        SELECT room_name FROM temporary_rooms
                        WHERE guild_snowflake = $1 AND room_snowflake = $2
                    ''', guild.id, channel.id)
                    room_name = temp_room['room_name'] if temp_room else ''
                    if isinstance(channel, discord.TextChannel):
                        for overwrite_obj, overwrite in channel.overwrites.items():
                            uid = overwrite_obj.id
                            if overwrite.send_messages is False and (uid, channel.id, room_name) not in text_set:
                                await channel.set_permissions(overwrite_obj, send_messages=None)
                            if overwrite.view_channel is False and (uid, channel.id, room_name) not in ban_set:
                                await channel.set_permissions(overwrite_obj, overwrite=None)
                    if isinstance(channel, discord.VoiceChannel):
                        for member in channel.members:
                            if member.voice and member.voice.mute and (member.id, channel.id, room_name) not in voice_set:
                                await member.edit(mute=False)
    
async def setup(bot: DiscordBot):
    await bot.add_cog(EventListeners(bot))
