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
from vyrtuous.utils.ban import Ban
from vyrtuous.utils.database import Database
from vyrtuous.utils.setup_logging import logger
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.stage import Stage
from vyrtuous.utils.video_room import VideoRoom
from vyrtuous.utils.voice_mute import VoiceMute

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
        if not self.update_video_room_status.is_running():
            self.update_video_room_status.start()
                                
    @tasks.loop(minutes=5)
    async def check_expired_bans(self):
        try:
            now = datetime.now(timezone.utc)
            async with self.bot.db_pool.acquire() as conn:
                expired_bans = await Ban.fetch_by_expired(now=now)
                for expired_ban in expired_bans:
                    try:
                        user_id = expired_ban.member_snowflake
                        guild_id = expired_ban.guild_snowflake
                        channel_id = expired_ban.channel_snowflake
                        guild = self.bot.get_guild(guild_id)
                        channel = self.bot.get_channel(channel_id)
                        if not channel and guild:
                            try:
                                channel = await self.bot.fetch_channel(channel_id)
                            except discord.NotFound:
                                channel = None
                        if channel is None:
                            logger.info(f'Channel {channel_id} not found, cleaning up expired ban')
                            await Ban.delete_by_channel_guild_and_member(channel_snowflake=channel_id, guild_snowflake=guild_id, member_snowflake=user_id)
                            continue
                        member = guild.get_member(user_id)
                        if member is None:
                            try:
                                member = await guild.fetch_member(user_id)
                            except discord.NotFound:
                                logger.info(f'Member {user_id} not found in guild {guild_id}, cleaning up expired ban')
                                await Ban.delete_by_channel_guild_and_member(channel_snowflake=channel_id, guild_snowflake=guild_id, member_snowflake=user_id)
                                continue
                        await Ban.delete_by_channel_guild_and_member(channel_snowflake=channel_id, guild_snowflake=guild_id, member_snowflake=user_id)
                        try:
                            await channel.set_permissions(member, overwrite=None)
                            logger.info(f'Removed ban override for user {user_id} in channel {channel_id}')
                        except discord.Forbidden:
                            logger.warning(f'No permission to remove ban override for user {user_id} in channel {channel_id}')
                        except discord.HTTPException as e:
                            logger.error(f'Failed to remove permission override for user {user_id} in channel {channel_id}: {e}.')
                    except Exception as e:
                        logger.error(f'Error processing expired ban for user {user_id} in guild {guild_id}: {e}.', exc_info=True)
                        continue
        except Exception as e:
            logger.error(f'Error in check_expired_bans task: {e}.', exc_info=True)
    
    @tasks.loop(seconds=15)
    async def check_expired_voice_mutes(self):
        try:
            now = datetime.now(timezone.utc)
            async with self.bot.db_pool.acquire() as conn:
                expired_voice_mutes = await VoiceMute.fetch_by_expired(now=now)
                for expired_voice_mute in expired_voice_mutes:
                    try:
                        guild_id = expired_voice_mute.guild_snowflake
                        user_id = expired_voice_mute.member_snowflake
                        channel_id = expired_voice_mute.channel_snowflake
                        target = expired_voice_mute.target
                        guild = self.bot.get_guild(guild_id)
                        channel = self.bot.get_channel(channel_id)
                        if not channel and guild:
                            try:
                                channel = await guild.fetch_channel(channel_id)
                            except discord.NotFound:
                                channel = None
                        if channel is None:
                            logger.info(f'Guild {guild_id} or channel {channel_id} not found, cleaning up expired voice mute')
                            await VoiceMute.delete_by_channel_guild_and_target(channel_snowflake=channel_id, guild_snowflake=guild_id, target=target)
                            continue
                        member = guild.get_member(user_id)
                        if member is None:
                            try:
                                member = await guild.fetch_member(user_id)
                            except discord.NotFound:
                                logger.info(f'Member {user_id} not found in guild {guild_id}, cleaning up expired voice mute')
                                await VoiceMute.delete_by_channel_guild_member_and_target(channel_snowflake=channel_id, guild_snowflake=guild_id, member_snowflake=member.id, target="user")
                                continue
                        await VoiceMute.delete_by_channel_guild_member_and_target(channel_snowflake=channel_id, guild_snowflake=guild_id, member_snowflake=member.id, target="user")
                        if member.voice and member.voice.channel and member.voice.channel.id == channel_id:
                            try:
                                await member.edit(mute=False)
                                logger.info(f'Unmuted user {user_id} in channel {channel_id}')
                            except discord.Forbidden:
                                logger.warning(f'No permission to unmute user {user_id} in channel {channel_id}')
                            except discord.HTTPException as e:
                                logger.error(f'Failed to unmute user {user_id} in channel {channel_id}: {e}.')
                        else:
                            logger.info(f'User {user_id} not in voice channel {channel_id}, skipping unmute')
                    except Exception as e:
                        logger.error(f'Error processing expired voice mute for user {user_id} in guild {guild_id}: {e}.', exc_info=True)
                        continue
        except Exception as e:
            logger.error(f'Error in check_expired_voice_mutes task: {e}.', exc_info=True)

    
    @tasks.loop(minutes=1)
    async def check_expired_stages(self):
        try:
            now = datetime.now(timezone.utc)
            async with self.bot.db_pool.acquire() as conn:
                expired_stages = await Stage.fetch_by_expired(now)
                for expired_stage in expired_stages:
                    try:
                        guild_id = expired_stage.guild_snowflake
                        channel_id = expired_stage.channel_snowflake
                        guild = self.bot.get_guild(guild_id)
                        voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_id, guild_snowflake=guild_id, target="room")
                        await VoiceMute.delete_by_channel_guild_and_target(channel_snowflake=channel_id, guild_snowflake=guild_id, target="room")
                        await Stage.delete_by_channel_and_guild(channel_snowflake=channel_id, guild_snowflake=guild_id)
                        if guild:
                            for voice_mute in voice_mutes:
                                user_id = voice_mute.member_snowflake
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
                                        logger.error(f'Failed to unmute user {user_id} in expired stage {channel_id}: {e}.')
                        else:
                            logger.info(f'Guild {guild_id} not found when processing expired stage {channel_id}')
                        logger.info(f'Cleaned up expired stage for channel {channel_id} in guild {guild_id}')
                    except Exception as e:
                        logger.error(f'Error processing expired stage for channel {channel_id} in guild {guild_id}: {e}.', exc_info=True)
                        continue
        except Exception as e:
            logger.error(f'Error in check_expired_stages task: {e}.', exc_info=True)

    @tasks.loop(minutes=5)
    async def update_video_room_status(self):
        video_rooms = await VideoRoom.fetch_all()
        for video_room in video_rooms:
            channel = await self.bot.fetch_channel(video_room.channel_snowflake)
            try:
                await channel.edit(status="Video-Only Room", reason="Reset video-only room status")
            except discord.Forbidden:
                logger.info("Failed to enforce video room status")
    
    @tasks.loop(minutes=1)
    async def check_expired_text_mutes(self):
        try:
            now = datetime.now(timezone.utc)
            async with self.bot.db_pool.acquire() as conn:
                expired_text_mutes = await TextMute.fetch_by_expired(now)
                if not expired_text_mutes:
                    return
                for expired_text_mute in expired_text_mutes:
                    try:
                        user_id = expired_text_mute.member_snowflake
                        guild_id = expired_text_mute.guild_snowflake
                        channel_id = expired_text_mute.channel_snowflake
                        guild = self.bot.get_guild(guild_id)
                        channel = self.bot.get_channel(channel_id)
                        if not channel and guild:
                            try:
                                channel = await guild.fetch_channel(channel_id)
                            except discord.NotFound:
                                channel = None
                        if guild is None or channel is None:
                            logger.info(f'Guild {guild_id} or channel {channel_id} not found, cleaning up expired text mute')
                            await TextMute.delete_by_channel_guild_and_member(channel_snowflake=channel_id, guild_snowflake=guild_id, member_snowflake=user_id)
                            continue
                        member = guild.get_member(user_id)
                        if member is None:
                            try:
                                member = await guild.fetch_member(user_id)
                            except discord.NotFound:
                                logger.info(f'Member {user_id} not found in guild {guild_id}, cleaning up expired text mute')
                                await TextMute.delete_by_channel_guild_and_member(channel_snowflake=channel_id, guild_snowflake=guild_id, member_snowflake=user_id)
                                continue
                        await TextMute.delete_by_channel_guild_and_member(channel_snowflake=channel_id, guild_snowflake=guild_id, member_snowflake=user_id)
                        try:
                            await channel.set_permissions(member, send_messages=None)
                            logger.info(f'Removed text mute override for user {user_id} in channel {channel_id}')
                        except discord.Forbidden:
                            logger.warning(f'No permission to remove mute override for user {user_id} in channel {channel_id}')
                        except discord.HTTPException as e:
                            logger.error(f'Failed to remove permission override for user {user_id} in channel {channel_id}: {e}.')
                    except Exception as e:
                        logger.error(f'Error processing expired text mute for user {user_id} in guild {guild_id}: {e}.', exc_info=True)
        except Exception as e:
            logger.error(f'Error in check_expired_text_mutes task: {e}.', exc_info=True)
                
    @tasks.loop(hours=24)
    async def backup_database(self):
        try:
            db = Database()
            db.create_backup_directory()
            db.execute_backup()
            logger.info(f'Backup completed successfully.')
        except Exception as e:
            logger.error(f'Error during database backup: {e}.')

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
    

    @update_video_room_status.before_loop
    async def before_update_video_room_status(self):
        await self.bot.wait_until_ready()

async def setup(bot: DiscordBot):
    await bot.add_cog(ScheduledTasks(bot))
