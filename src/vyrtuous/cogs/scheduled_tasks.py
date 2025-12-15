''' scheduled_tasks.py A discord.py cog containing scheduled tasks for the Vyrtuous bot.

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
from datetime import datetime, timezone
from discord.ext import commands, tasks
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.utils.database import Database
from vyrtuous.utils.setup_logging import logger

import discord

class ScheduledTasks(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        
    async def cog_load(self):
        if not self.backup_database.is_running():
            self.backup_database.start()
        if not self.check_expired_bans.is_running():
            self.check_expired_bans.start()
        if not self.check_expired_voice_mutes.is_running():
            self.check_expired_voice_mutes.start()
        if not self.check_expired_text_mutes.is_running():
            self.check_expired_text_mutes.start()
        if not self.check_expired_stages.is_running():
            self.check_expired_stages.start()
        if not self.check_active_bans.is_running():
            self.check_active_bans.start()

    @tasks.loop(minutes=5)
    async def check_active_bans(self):
        async with self.bot.db_pool.acquire() as conn:
            for guild in self.bot.guilds:
                for channel in guild.channels:
                    if not isinstance(channel, (discord.TextChannel, discord.VoiceChannel)):
                        continue
                    for member in channel.members:
                        perms = channel.permissions_for(member)
                        if not perms.view_channel:
                            break
                        ban = await conn.fetchrow('''
                            SELECT expires_at
                            FROM active_bans
                            WHERE guild_id = $1
                              AND discord_snowflake = $2
                              AND channel_id = $3
                              AND room_name = $4
                        ''', guild.id, member.id, channel.id, channel.name)
                        if ban:
                            expires_at = ban['expires_at']
                            now = datetime.now(timezone.utc)
                            if expires_at is None or expires_at > now:
                                await channel.set_permissions(member, overwrite=discord.PermissionOverwrite(view_channel=False))
                                
    @tasks.loop(minutes=5)
    async def check_expired_bans(self):
        try:
            now = datetime.now(timezone.utc)
            async with self.bot.db_pool.acquire() as conn:
                expired = await conn.fetch('''
                    SELECT guild_id, discord_snowflake, channel_id, room_name
                    FROM active_bans
                    WHERE expires_at <= $1
                ''', now)
                for record in expired:
                    try:
                        user_id = record['discord_snowflake']
                        guild_id = record['guild_id']
                        channel_id = record['channel_id']
                        room_name = record['room_name']
                        guild = self.bot.get_guild(guild_id)
                        channel = self.bot.get_channel(channel_id)
                        if not channel and guild:
                            try:
                                channel = await self.bot.fetch_channel(channel_id)
                            except discord.NotFound:
                                channel = None
                        if channel is None:
                            logger.info(f'Channel {channel_id} not found, cleaning up expired ban')
                            await conn.execute('''
                                DELETE FROM active_bans
                                WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND room_name=$4
                            ''', guild_id, user_id, channel_id, room_name)
                            continue
                        member = guild.get_member(user_id)
                        if member is None:
                            try:
                                member = await guild.fetch_member(user_id)
                            except discord.NotFound:
                                logger.info(f'Member {user_id} not found in guild {guild_id}, cleaning up expired ban')
                                await conn.execute('''
                                    DELETE FROM active_bans
                                    WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND room_name=$4
                                ''', guild_id, user_id, channel_id, room_name)
                                continue
                        await conn.execute('''
                            DELETE FROM active_bans
                            WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND room_name=$4
                        ''', guild_id, user_id, channel_id, room_name)
                        try:
                            await channel.set_permissions(member, overwrite=None)
                            logger.info(f'Removed ban override for user {user_id} in channel {channel_id}')
                        except discord.Forbidden:
                            logger.warning(f'No permission to remove ban override for user {user_id} in channel {channel_id}')
                        except discord.HTTPException as e:
                            logger.error(f'Failed to remove permission override for user {user_id} in channel {channel_id}: {e}')
                    except Exception as e:
                        logger.error(f'Error processing expired ban for user {user_id} in guild {guild_id}: {e}', exc_info=True)
                        continue
        except Exception as e:
            logger.error(f'Error in check_expired_bans task: {e}', exc_info=True)
    
    @tasks.loop(seconds=15)
    async def check_expired_voice_mutes(self):
        try:
            now = datetime.now(timezone.utc)
            async with self.bot.db_pool.acquire() as conn:
                expired = await conn.fetch('''
                    SELECT guild_id, discord_snowflake, channel_id, target, room_name
                    FROM active_voice_mutes
                    WHERE expires_at IS NOT NULL AND expires_at <= $1
                ''', now)
                for record in expired:
                    try:
                        guild_id = record['guild_id']
                        user_id = record['discord_snowflake']
                        channel_id = record['channel_id']
                        target = record['target']
                        room_name = record['room_name']
                        guild = self.bot.get_guild(guild_id)
                        channel = self.bot.get_channel(channel_id)
                        if not channel and guild:
                            try:
                                channel = await guild.fetch_channel(channel_id)
                            except discord.NotFound:
                                channel = None
                        if channel is None:
                            logger.info(f'Guild {guild_id} or channel {channel_id} not found, cleaning up expired voice mute')
                            await conn.execute('''
                                DELETE FROM active_voice_mutes
                                WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target=$4 AND room_name=$5
                            ''', guild_id, user_id, channel_id, target, room_name)
                            continue
                        member = guild.get_member(user_id)
                        if member is None:
                            try:
                                member = await guild.fetch_member(user_id)
                            except discord.NotFound:
                                logger.info(f'Member {user_id} not found in guild {guild_id}, cleaning up expired voice mute')
                                await conn.execute('''
                                    DELETE FROM active_voice_mutes
                                    WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target=$4 AND room_name=$5
                                ''', guild_id, user_id, channel_id, target, room_name)
                                continue
                        await conn.execute('''
                            DELETE FROM active_voice_mutes
                            WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target=$4 AND room_name=$5
                        ''', guild_id, user_id, channel_id, target, room_name)
                        if member.voice and member.voice.channel and member.voice.channel.id == channel_id:
                            try:
                                await member.edit(mute=False)
                                logger.info(f'Unmuted user {user_id} in channel {channel_id}')
                            except discord.Forbidden:
                                logger.warning(f'No permission to unmute user {user_id} in channel {channel_id}')
                            except discord.HTTPException as e:
                                logger.error(f'Failed to unmute user {user_id} in channel {channel_id}: {e}')
                        else:
                            logger.info(f'User {user_id} not in voice channel {channel_id}, skipping unmute')
                    except Exception as e:
                        logger.error(f'Error processing expired voice mute for user {user_id} in guild {guild_id}: {e}', exc_info=True)
                        continue
        except Exception as e:
            logger.error(f'Error in check_expired_voice_mutes task: {e}', exc_info=True)

    
    @tasks.loop(minutes=1)
    async def check_expired_stages(self):
        try:
            now = datetime.now(timezone.utc)
            async with self.bot.db_pool.acquire() as conn:
                expired = await conn.fetch('''
                    SELECT guild_id, channel_id, room_name
                    FROM active_stages
                    WHERE expires_at IS NOT NULL
                      AND expires_at <= $1
                ''', now)
                for record in expired:
                    try:
                        guild_id = record['guild_id']
                        channel_id = record['channel_id']
                        room_name = record['room_name']
                        guild = self.bot.get_guild(guild_id)
                        muted_members = await conn.fetch('''
                            SELECT discord_snowflake
                            FROM active_voice_mutes
                            WHERE guild_id = $1 AND channel_id = $2 AND target = $3 AND room_name = $4
                        ''', guild_id, channel_id, 'room', room_name)
                        await conn.execute('''
                            DELETE FROM active_voice_mutes
                            WHERE guild_id = $1 AND channel_id = $2 AND target = 'room' AND room_name = $3
                        ''', guild_id, channel_id, room_name)
                        await conn.execute('''
                            DELETE FROM stage_coordinators
                            WHERE guild_id = $1 AND channel_id = $2 AND room_name = $3
                        ''', guild_id, channel_id, room_name)
                        await conn.execute('''
                            DELETE FROM active_stages
                            WHERE guild_id = $1 AND channel_id = $2 AND room_name = $3
                        ''', guild_id, channel_id, room_name)
                        if guild:
                            for member_record in muted_members:
                                user_id = member_record['discord_snowflake']
                                member = guild.get_member(user_id)
                                if member is None:
                                    try:
                                        member = await guild.fetch_member(user_id)
                                    except discord.NotFound:
                                        logger.error(f'Failed to retrieve member {user_id} from expired stage.')
                                        continue
                                if member and member.voice and member.voice.mute:
                                    try:
                                        await member.edit(mute=False, reason='Stage room closed or unmuted automatically')
                                        logger.info(f'Unmuted user {user_id} after stage {channel_id} expired')
                                    except discord.Forbidden:
                                        logger.warning(f'No permission to unmute user {user_id} in expired stage {channel_id}')
                                    except discord.HTTPException as e:
                                        logger.error(f'Failed to unmute user {user_id} in expired stage {channel_id}: {e}')
                        else:
                            logger.info(f'Guild {guild_id} not found when processing expired stage {channel_id}')
                        logger.info(f'Cleaned up expired stage for channel {channel_id} in guild {guild_id}')
                    except Exception as e:
                        logger.error(f'Error processing expired stage for channel {channel_id} in guild {guild_id}: {e}', exc_info=True)
                        continue
        except Exception as e:
            logger.error(f'Error in check_expired_stages task: {e}', exc_info=True)

    
    @tasks.loop(minutes=1)
    async def check_expired_text_mutes(self):
        try:
            now = datetime.now(timezone.utc)
            async with self.bot.db_pool.acquire() as conn:
                expired = await conn.fetch('''
                    SELECT guild_id, discord_snowflake, channel_id, room_name
                    FROM active_text_mutes
                    WHERE expires_at IS NOT NULL AND expires_at <= $1
                ''', now)
                if not expired:
                    return
                for record in expired:
                    try:
                        user_id = record['discord_snowflake']
                        guild_id = record['guild_id']
                        channel_id = record['channel_id']
                        room_name = record['room_name']
                        guild = self.bot.get_guild(guild_id)
                        channel = self.bot.get_channel(channel_id)
                        if not channel and guild:
                            try:
                                channel = await guild.fetch_channel(channel_id)
                            except discord.NotFound:
                                channel = None
                        if guild is None or channel is None:
                            logger.info(f'Guild {guild_id} or channel {channel_id} not found, cleaning up expired text mute')
                            await conn.execute('''
                                DELETE FROM active_text_mutes
                                WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND room_name=$4
                            ''', guild_id, user_id, channel_id, room_name)
                            continue
                        member = guild.get_member(user_id)
                        if member is None:
                            try:
                                member = await guild.fetch_member(user_id)
                            except discord.NotFound:
                                logger.info(f'Member {user_id} not found in guild {guild_id}, cleaning up expired text mute')
                                await conn.execute('''
                                    DELETE FROM active_text_mutes
                                    WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND room_name=$4
                                ''', guild_id, user_id, channel_id, room_name)
                                continue
                        await conn.execute('''
                            DELETE FROM active_text_mutes
                            WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND room_name=$4
                        ''', guild_id, user_id, channel_id, room_name)
                        try:
                            await channel.set_permissions(member, send_messages=None)
                            logger.info(f'Removed text mute override for user {user_id} in channel {channel_id}')
                        except discord.Forbidden:
                            logger.warning(f'No permission to remove mute override for user {user_id} in channel {channel_id}')
                        except discord.HTTPException as e:
                            logger.error(f'Failed to remove permission override for user {user_id} in channel {channel_id}: {e}')
                    except Exception as e:
                        logger.error(f'Error processing expired text mute for user {user_id} in guild {guild_id}: {e}', exc_info=True)
        except Exception as e:
            logger.error(f'Error in check_expired_text_mutes task: {e}', exc_info=True)
                
    @tasks.loop(hours=24)
    async def backup_database(self) -> None:
        try:
            backup = Database()
            backup.create_backup_directory()
            backup.execute_backup()
            logger.info(f'Backup completed successfully.')
        except Exception as e:
            logger.error(f'Error during database backup: {e}')

    @backup_database.before_loop
    async def before_backup(self):
        await self.bot.wait_until_ready()
    
    @check_expired_bans.before_loop
    async def before_check_expired_bans(self):
        await self.bot.wait_until_ready()
    
    @check_expired_voice_mutes.before_loop
    async def before_check_expired_voice_mutes(self):
        await self.bot.wait_until_ready()
    
    @check_expired_text_mutes.before_loop
    async def before_check_expired_text_mutes(self):
        await self.bot.wait_until_ready()
    
    @check_expired_stages.before_loop
    async def before_check_expired_stages(self):
        await self.bot.wait_until_ready()

async def setup(bot: DiscordBot):
    await bot.add_cog(ScheduledTasks(bot))
