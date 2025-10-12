''' scheduled_tasks.py

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
from datetime import datetime, timezone
import os
import subprocess

from discord.ext import commands, tasks
from vyrtuous.inc.helpers import *
from vyrtuous.utils.setup_logging import logger
from vyrtuous.bot.discord_bot import DiscordBot

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
    
    @staticmethod
    def perform_backup(db_user: str, db_name: str, db_host: str, db_password: str, backup_dir: str) -> str:
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
    
    @staticmethod
    def setup_backup_directory(backup_dir: str) -> str:
        os.makedirs(backup_dir, exist_ok=True)
        return backup_dir
        
    @commands.after_invoke
    async def after_invoke(self, ctx: commands.Context) -> None:
        if hasattr(self.bot, 'db_pool'):
            await self.bot.db_pool.close()

    @tasks.loop(minutes=5)
    async def check_expired_bans(self):
        now = datetime.now(timezone.utc)
        async with self.bot.db_pool.acquire() as conn:
            expired = await conn.fetch('''
                SELECT guild_id, discord_snowflake, channel_id
                FROM active_bans
                WHERE expires_at <= $1
            ''', now)
        for record in expired:
            user_id = record['discord_snowflake']
            guild = self.bot.get_guild(record['guild_id'])
            if guild is None:
                continue
            channel_id = record['channel_id']
            channel = self.bot.get_channel(channel_id)
            if not channel:
                try:
                    channel = await guild.fetch_channel(channel_id)
                except discord.NotFound:
                    continue
            member = guild.get_member(user_id)
            if member is None:
                try:
                    member = await guild.fetch_member(user_id)
                except discord.NotFound:
                    continue
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute('''
                    DELETE FROM active_bans
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                ''', guild.id, user_id, channel_id)
            try:
                await channel.set_permissions(member, overwrite=None)
            except discord.Forbidden:
                logger.warning(f'No permission to remove ban override for user {user_id} in channel {channel_id}.')
            except discord.HTTPException as e:
                logger.error(f'Failed to remove permission override: {e}')
    
    @tasks.loop(seconds=15)
    async def check_expired_voice_mutes(self):
        now = datetime.now(timezone.utc)
        async with self.bot.db_pool.acquire() as conn:
            expired = await conn.fetch('''
                SELECT guild_id, discord_snowflake, channel_id, target
                FROM active_voice_mutes
                WHERE expires_at IS NOT NULL
                  AND expires_at <= $1
            ''', now)
            for record in expired:
                guild = self.bot.get_guild(record['guild_id'])
                if guild is None:
                    continue
                channel = self.bot.get_channel(record['channel_id'])
                if not channel:
                    try:
                        channel = await guild.fetch_channel(record['channel_id'])
                    except discord.NotFound:
                        continue
                if not isinstance(channel, discord.VoiceChannel):
                    continue
                member = guild.get_member(record['discord_snowflake'])
                if member is None:
                    try:
                        member = await guild.fetch_member(record['discord_snowflake'])
                    except discord.NotFound:
                        continue
                if member.voice and member.voice.channel and member.voice.channel.id == record['channel_id']:
                    try:
                        await member.edit(mute=False)
                    except discord.HTTPException:
                        logger.warning(f'No permission to unmute user {member.id} in channel {channel.id}.')
                await conn.execute('''
                    DELETE FROM active_voice_mutes
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3 AND target = $4
                ''', record['guild_id'], record['discord_snowflake'], record['channel_id'], record['target'])
    
    @tasks.loop(minutes=1)
    async def check_expired_stages(self):
        now = datetime.now(timezone.utc)
        async with self.bot.db_pool.acquire() as conn:
            expired = await conn.fetch('''
                SELECT guild_id, channel_id
                FROM active_stages
                WHERE expires_at IS NOT NULL
                  AND expires_at <= $1
            ''', now)
            for record in expired:
                guild_id, channel_id = record['guild_id'], record['channel_id']
                records = await conn.fetch('SELECT discord_snowflake FROM active_voice_mutes WHERE guild_id = $1 AND channel_id = $2 AND target = $3', guild_id, channel_id, 'room')
                await conn.execute('''
                    DELETE FROM active_voice_mutes
                    WHERE guild_id = $1 AND channel_id = $2 AND target = 'room'
                ''', guild_id, channel_id)
                guild = self.bot.get_guild(guild_id)
                if guild:
                    for record in records:
                        member = guild.get_member(record['user_id'])
                        if member and member.voice and member.voice.mute:
                            await member.edit(mute=False, reason='Stage room closed or unmuted automatically')
                await conn.execute('''
                    DELETE FROM stage_coordinators
                    WHERE guild_id = $1 AND channel_id = $2
                ''', guild_id, channel_id)
                await conn.execute('''
                    DELETE FROM active_stages
                    WHERE guild_id = $1 AND channel_id = $2
                ''', guild_id, channel_id)
    
    @tasks.loop(minutes=1)
    async def check_expired_text_mutes(self):
        now = datetime.now(timezone.utc)
        async with self.bot.db_pool.acquire() as conn:
            expired = await conn.fetch('''
                SELECT guild_id, discord_snowflake, channel_id
                FROM active_text_mutes
                WHERE expires_at IS NOT NULL
                  AND expires_at <= $1
            ''', now)
            if not expired:
                return
            for record in expired:
                user_id = record['discord_snowflake']
                guild = self.bot.get_guild(record['guild_id'])
                if guild is None:
                    continue
                channel_id = record['channel_id']
                channel = self.bot.get_channel(channel_id)
                if not channel:
                    try:
                        channel = await guild.fetch_channel(channel_id)
                    except discord.NotFound:
                        continue
                if not isinstance(channel, discord.VoiceChannel):
                    continue
                member = guild.get_member(user_id)
                if member is None:
                    try:
                        member = await guild.fetch_member(user_id)
                    except discord.NotFound:
                        continue
                await conn.execute('''
                    DELETE FROM active_text_mutes
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                ''', guild.id, member.id, channel_id)
                try:
                    await channel.set_permissions(member, send_messages=None)
                except discord.Forbidden:
                    logger.warning(f'No permission to remove mute override for user {user_id} in channel {channel_id}.')
                except discord.HTTPException as e:
                    logger.error(f'Failed to remove permission override: {e}')

    @tasks.loop(hours=24)
    async def backup_database(self) -> None:
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
        except Exception as e:
            logger.error(f'Error during database backup: {e}')

    @backup_database.before_loop
    async def before_backup(self):
        await self.bot.wait_until_ready()

async def setup(bot: DiscordBot):
    await bot.add_cog(ScheduledTasks(bot))
