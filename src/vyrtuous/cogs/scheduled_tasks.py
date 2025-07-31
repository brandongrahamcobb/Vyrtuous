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
import datetime
import os
import subprocess

from discord.ext import commands, tasks
from vyrtuous.inc.helpers import *
from vyrtuous.utils.setup_logging import logger


class ScheduledTasks(commands.Cog):

    def __init__(self, bot):
        self.backup_database.start()
        self.bot = bot
        self.config = bot.config
        self.check_expired_bans.start()
        self.check_expired_voice_mutes.start()
    
    def perform_backup(self, db_user: str, db_name: str, db_host: str, db_password: str, backup_dir: str) -> str:
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        backup_file = os.path.join(backup_dir, f'backup_{timestamp}.sql')
        dump_command = [
            'pg_dump',
            '-U', db_user,
            '-h', db_host,
            '-d', db_name,
            '-F', 'p',
            '-f', backup_file,
        ]
        # Pass PGPASSWORD via environment for security
        env = os.environ.copy()
        env["PGPASSWORD"] = db_password
        result = subprocess.run(
            dump_command,
            capture_output=True,
            text=True,
            env=env,
        )
        if result.returncode != 0:
            raise RuntimeError(f'Backup failed: {result.stderr}')
        return backup_file
    
    def setup_backup_directory(self, backup_dir: str) -> str:
        os.makedirs(backup_dir, exist_ok=True)
        return backup_dir
        
    @commands.after_invoke
    async def after_invoke(self, ctx) -> None:
        if hasattr(self.bot, 'db_pool'):
            await self.bot.db_pool.close()
            
    # @tasks.loop(minutes=5)
    # async def check_expired_bans(self):
    #     now = datetime.datetime.utcnow()
    #     async with self.bot.db_pool.acquire() as conn:
    #         rows = await conn.fetch(
    #             '''
    #             SELECT user_id, channel_id
    #             FROM ban_expirations
    #             WHERE expires_at <= $1
    #             ''',
    #             now
    #         )
    #         for row in rows:
    #             user_id = row['user_id']
    #             channel_id = row['channel_id']
    #             channel = self.bot.get_channel(channel_id)
    #             if not isinstance(channel, discord.VoiceChannel):
    #                 continue
    #             guild = channel.guild
    #             member = guild.get_member(user_id)
    #             if member and channel:
    #                 await self.remove_ban_role(member, channel)
    #             await conn.execute(
    #                 '''
    #                 DELETE FROM ban_expirations
    #                 WHERE user_id = $1 AND channel_id = $2
    #                 ''',
    #                 user_id, channel_id
    #             )
    #
    # async def remove_ban_role(self, member: discord.Member, channel: discord.VoiceChannel):
    #     role_id = self.bot.command_aliases.get(member.guild.id, {}).get('role', {}).get(str(channel.id))
    #     if not role_id:
    #         return
    #     role = member.guild.get_role(role_id)
    #     if role and role in member.roles:
    #         await member.remove_roles(role, reason='Ban expired')

    async def remove_ban_role(self, member: discord.Member, channel: discord.VoiceChannel):
        guild_id = member.guild.id
        channel_id = channel.id

        async with self.bot.db_pool.acquire() as conn:
            role_id = await conn.fetchval(
                '''
                SELECT role_id
                FROM ban_roles
                WHERE guild_id = $1
                  AND channel_id = $2
                ''',
                guild_id, channel_id
            )
        if not role_id:
            return

        role = member.guild.get_role(role_id)
        if role and role in member.roles:
            try:
                await member.remove_roles(role, reason='Ban expired')
            except discord.Forbidden:
                logger.warning(f"❌ Missing permissions to remove role {role_id} from {member}.")

    @tasks.loop(minutes=5)
    async def check_expired_bans(self):
        now = datetime.datetime.utcnow()
        async with self.bot.db_pool.acquire() as conn:
            expired = await conn.fetch(
                '''
                SELECT user_id, channel_id
                FROM ban_expirations
                WHERE expires_at <= $1
                ''',
                now
            )
        for record in expired:
            user_id = record['user_id']
            channel_id = record['channel_id']
            channel = self.bot.get_channel(channel_id)
            if not isinstance(channel, discord.VoiceChannel):
                continue
            guild = channel.guild
            member = guild.get_member(user_id)
            if not member:
                continue
            await self.remove_ban_role(member, channel)
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute(
                    'DELETE FROM ban_expirations WHERE user_id = $1 AND channel_id = $2',
                    user_id, channel_id
                )
                await conn.execute(
                    'DELETE FROM active_bans WHERE user_id = $1 AND channel_id = $2',
                    user_id, channel_id
                )
                await conn.execute(
                    '''
                    UPDATE users
                    SET ban_channel_ids = array_remove(ban_channel_ids, $2),
                        updated_at      = NOW()
                    WHERE user_id = $1
                    ''',
                    user_id, channel_id
                )

    @tasks.loop(minutes=5)
    async def check_expired_voice_mutes(self):
        now = datetime.datetime.utcnow()
        async with self.bot.db_pool.acquire() as conn:
            expired = await conn.fetch(
                '''
                SELECT user_id, channel_id
                FROM active_mutes
                WHERE expires_at IS NOT NULL
                  AND expires_at <= $1
                ''',
                now
            )
        for record in expired:
            user_id = record['user_id']
            channel_id = record['channel_id']
            channel = self.bot.get_channel(channel_id)
            if not isinstance(channel, discord.VoiceChannel):
                continue
            guild = channel.guild
            member = guild.get_member(user_id)
            if not member:
                continue
            if member.voice and member.voice.channel and member.voice.channel.id == channel_id:
                try:
                    await member.edit(mute=False)
                except discord.HTTPException:
                    pass
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute(
                    'DELETE FROM active_mutes WHERE user_id = $1 AND channel_id = $2',
                    user_id, channel_id
                )
                await conn.execute(
                    '''
                    UPDATE users
                    SET mute_channel_ids = array_remove(mute_channel_ids, $2),
                        updated_at       = NOW()
                    WHERE user_id = $1
                    ''',
                    user_id, channel_id
                )
                await conn.execute(
                    '''
                    UPDATE users
                    SET server_mute_channel_ids = array_remove(server_mute_channel_ids, $2),
                        updated_at              = NOW()
                    WHERE user_id = $1
                    ''',
                    user_id, channel_id
                )
    @tasks.loop(hours=24)
    async def backup_database(self) -> None:
        try:
            backup_dir = self.setup_backup_directory('./backups')
            backup_file = self.perform_backup(
                db_user=os.getenv("POSTGRES_USER"),
                db_name=os.getenv("POSTGRES_DATABASE"),
                db_host=os.getenv("POSTGRES_HOST"),
                db_password=os.getenv("POSTGRES_PASSWORD"),
                backup_dir=backup_dir
            )
            logger.info(f'Backup completed successfully: {backup_file}')
        except Exception as e:
            logger.error(f'Error during database backup: {e}')

    @backup_database.before_loop
    async def before_backup(self):
        await self.bot.wait_until_ready()
        
    @check_expired_bans.before_loop
    async def before_backup(self):
        await self.bot.wait_until_ready()

    @check_expired_voice_mutes.before_loop
    async def before_backup(self):
        await self.bot.wait_until_ready()

async def setup(bot: commands.Bot):
    await bot.add_cog(ScheduledTasks(bot))
