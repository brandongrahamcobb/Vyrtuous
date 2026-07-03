"""!/bin/python3
scheduled_tasks.py A discord.py cog containing scheduled tasks for the Vyrtuous bot.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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
"""

import discord
from discord.ext import commands, tasks

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database import Database
from vyrtuous.utils.moderation import (ban_service, text_mute_service,
                                       voice_mute_service)
from vyrtuous.utils.rooms import automute_room_service
from vyrtuous.utils.statistics import system_monitoring_service
from vyrtuous.utils.tracking import bug_service
from vyrtuous.utils.users import active_member_service


class ScheduledTasks(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.__bot = bot

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
        if not self.check_expired_bugs.is_running():
            self.check_expired_bugs.start()
        if not self.temporarily_cleanup_overwrites.is_running():
            self.temporarily_cleanup_overwrites.start()
        if not self.save_active_members.is_running():
            self.save_active_members.start()
        if not self.system_monitoring.is_running():
            self.system_monitoring.start()
        if not self.cleanup_stale_overwrites.is_running():
            self.cleanup_stale_overwrites.start()

    @tasks.loop(hours=1)
    async def remove_inactive_members(self):
        for guild in self.__bot.guilds:
            count = await active_member_service.remove_inactive_members(guild=guild)
            self.__bot.logger.info(
                f"Removed {count} inactive members from {guild.name}."
            )

    @tasks.loop(minutes=5)
    async def system_monitoring(self):
        await system_monitoring_service.log_cpu_seconds()
        await system_monitoring_service.log_rx_bytes()
        await system_monitoring_service.log_tx_bytes()

    @tasks.loop(hours=72)
    async def cleanup_stale_overwrites(self):
        for guild in self.__bot.guilds:
            for channel in guild.channels:
                if isinstance(channel, discord.VoiceChannel):
                    for target, overwrite in channel.overwrites.items():
                        if overwrite.is_empty():
                            try:
                                await channel.set_permissions(target, overwrite=None)
                                self.__bot.logger.info(
                                    f"Cleaned up stale overwrite for {target.mention} in {channel.mention}."
                                )
                            except discord.HTTPException:
                                self.__bot.logger.info(
                                    f"Failed to cleaned up stale overwrite for {target.mention} in {channel.mention}."
                                )
                            except discord.Forbidden:
                                self.__bot.logger.info(
                                    f"Failed to cleaned up stale overwrite for {target.mention} in {channel.mention}."
                                )

    @tasks.loop(minutes=1)
    async def save_active_members(self):
        await active_member_service.save_active_members()
        self.__bot.logger.info("Saved active members.")

    @tasks.loop(minutes=5)
    async def check_expired_bans(self):
        await ban_service.clean_expired()
        self.__bot.logger.info("Cleaned up expired bans.")

    @tasks.loop(seconds=15)
    async def check_expired_voice_mutes(self):
        await voice_mute_service.clean_expired()
        self.__bot.logger.info("Cleaned up expired voice-mutes.")

    @tasks.loop(minutes=1)
    async def check_expired_stages(self):
        await automute_room_service.clean_expired()
        self.__bot.logger.info("Cleaned up expired stages.")

    @tasks.loop(hours=8)
    async def check_expired_bugs(self):
        await bug_service.clean_expired()
        self.__bot.logger.info("Sent developer log to developers.")

    @tasks.loop(minutes=1)
    async def check_expired_text_mutes(self):
        await text_mute_service.clean_expired()
        self.__bot.logger.info("Cleaned up expired text-mutes.")

    @tasks.loop(hours=1)
    async def temporarily_cleanup_overwrites(self):
        await ban_service.clean_overwrites()
        await text_mute_service.clean_overwrites()
        self.__bot.logger.info("Reset ban and text-mute overwrites.")

    @tasks.loop(hours=24)
    async def backup_database(self):
        try:
            db = Database(config=self.__bot.config)
            db.create_backup_directory()
            db.execute_backup()
            self.__bot.logger.info("Backup completed successfully.")
        except Exception as e:
            self.__bot.logger.error(
                f"Error during database backup: {str(e).capitalize()}"
            )


async def setup(bot: DiscordBot):
    await bot.add_cog(ScheduledTasks(bot))
