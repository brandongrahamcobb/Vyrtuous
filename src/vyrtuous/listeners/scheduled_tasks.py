"""!/bin/python3
scheduled_tasks.py A discord.py cog containing scheduled tasks for the Vyrtuous bot.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import time

import discord
from discord.ext import commands, tasks

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import ChannelState
from vyrtuous.clean import (
    clean_active_member_service,
    clean_automute_service,
    clean_ban_service,
    clean_text_mute_service,
    clean_video_only_channels_service,
    clean_voice_mute_service,
)
from vyrtuous.db.database import Database
from vyrtuous.utils.moderation import flag_service
from vyrtuous.utils.statistics import system_monitoring_service
from vyrtuous.utils.users import active_member_service


class ScheduledTasks(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def cog_load(self) -> None:
        if not self.backup_database.is_running():
            self.backup_database.start()
        if not self.check_expired_bans.is_running():
            self.check_expired_bans.start()
        if not self.check_expired_voice_mutes.is_running():
            self.check_expired_voice_mutes.start()
        if not self.check_expired_text_mutes.is_running():
            self.check_expired_text_mutes.start()
        if not self.check_expired_automutes.is_running():
            self.check_expired_automutes.start()
        if not self.temporarily_cleanup_overwrites.is_running():
            self.temporarily_cleanup_overwrites.start()
        if not self.save_active_members.is_running():
            self.save_active_members.start()
        if not self.system_monitoring.is_running():
            self.system_monitoring.start()
        if not self.cleanup_stale_overwrites.is_running():
            self.cleanup_stale_overwrites.start()
        if not self.notify_loop.is_running():
            self.notify_loop.start()
        if not self.check_expired_video_only_channels.is_running():
            self.check_expired_video_only_channels.start()
        if not self.clean_inactive_members.is_running():
            self.clean_inactive_members.start()

    @tasks.loop(hours=1)
    async def clean_inactive_members(self) -> None:
        count = await clean_active_member_service.clean_inactive_members()
        self.__bot.logger.debug(f"Removed {count} inactive members.")

    @tasks.loop(minutes=5)
    async def system_monitoring(self) -> None:
        await system_monitoring_service.log_cpu_seconds()
        await system_monitoring_service.log_rx_bytes()
        await system_monitoring_service.log_tx_bytes()

    @tasks.loop(hours=72)
    async def cleanup_stale_overwrites(self) -> None:
        for guild in self.__bot.guilds:
            for channel in guild.channels:
                if isinstance(channel, discord.VoiceChannel):
                    for target, overwrite in channel.overwrites.items():
                        if overwrite.is_empty() and isinstance(target, discord.Member):
                            try:
                                await channel.set_permissions(target, overwrite=None)
                                self.__bot.logger.debug(
                                    f"Cleaned up stale overwrite for {target.mention} in {channel.mention}."
                                )
                            except discord.HTTPException:
                                self.__bot.logger.debug(
                                    f"Failed to cleaned up stale overwrite for {target.mention} in {channel.mention}."
                                )

    @tasks.loop(seconds=30)
    async def notify_loop(self) -> None:
        bot: DiscordBot = DiscordBot.get_instance()
        joined_at = bot.registry.get(ChannelState).joined_at
        now = time.time()
        for (
            channel_snowflake,
            guild_snowflake,
            member_snowflake,
        ), joined in joined_at.items():
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None or not isinstance(
                channel, discord.channel.VocalGuildChannel
            ):
                continue
            member = guild.get_member(member_snowflake)
            if member is None:
                continue
            if member in channel.members:
                elapsed = now - joined
                if elapsed >= 900 and int(elapsed // 900) > int((elapsed - 60) // 900):
                    await flag_service.warn(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        member_snowflake=member_snowflake,
                    )

    @tasks.loop(minutes=1)
    async def save_active_members(self) -> None:
        await active_member_service.save_and_update_active_members()
        self.__bot.logger.debug("Saved active members.")

    @tasks.loop(minutes=5)
    async def check_expired_bans(self) -> None:
        await clean_ban_service.clean_expired_bans()
        self.__bot.logger.debug("Cleaned up expired bans.")

    @tasks.loop(seconds=15)
    async def check_expired_voice_mutes(self) -> None:
        await clean_voice_mute_service.clean_expired_voice_mutes()
        self.__bot.logger.debug("Cleaned up expired voice-mutes.")

    @tasks.loop(minutes=1)
    async def check_expired_automutes(self) -> None:
        await clean_automute_service.clean_expired_automutes()
        self.__bot.logger.debug("Cleaned up expired automute channels.")

    @tasks.loop(minutes=1)
    async def check_expired_video_only_channels(self) -> None:
        await clean_video_only_channels_service.clean_expired_video_only_channels()
        self.__bot.logger.debug("Cleaned up expired video-only channels.")

    @tasks.loop(minutes=1)
    async def check_expired_text_mutes(self) -> None:
        await clean_text_mute_service.clean_expired_text_mutes()
        self.__bot.logger.debug("Cleaned up expired text-mutes.")

    @tasks.loop(hours=1)
    async def temporarily_cleanup_overwrites(self) -> None:
        await clean_ban_service.clean_ban_overwrites()
        await clean_text_mute_service.clean_text_mute_overwrites()
        self.__bot.logger.debug("Reset ban and text-mute overwrites.")

    @tasks.loop(hours=24)
    async def backup_database(self) -> None:
        try:
            db = Database(config=self.__bot.config)
            db.create_backup_directory()
            db.execute_backup()
            self.__bot.logger.debug("Backup completed successfully.")
        except Exception as e:
            self.__bot.logger.error(
                f"Error during database backup: {str(e).capitalize()}"
            )

async def setup(bot: DiscordBot):
    await bot.add_cog(ScheduledTasks(bot))
