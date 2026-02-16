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

from datetime import datetime, timedelta, timezone

import discord
from discord.ext import commands, tasks

from vyrtuous.ban.ban_service import BanService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.database import Database
from vyrtuous.developer.developer import Developer
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_service import DurationService

# from vyrtuous.video_room.video_room import VideoRoom
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.stage_room.stage import Stage
from vyrtuous.stage_room.stage_service import StageService
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.text_mute.text_mute import TextMute
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.logger import logger
from vyrtuous.voice_mute.voice_mute import VoiceMute
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService


class ScheduledTasks(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.__bot = bot
        self.config = bot.config
        self.__author_service = AuthorService()
        self.__bot = bot
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__duration_service = DurationService()
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__stream_service = StreamService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
        )
        self.__ban_service = BanService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__text_mute_service = TextMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__voice_mute_service = VoiceMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__stage_service = StageService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__guild_owner_service = GuildOwnerService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__developer_service = DeveloperService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__sysadmin_service = SysadminService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__bug_service = BugService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )

    async def cog_load(self):
        if not self.backup_database.is_running():
            self.backup_database.start()
        if not self.check_expired_bans.is_running():
            self.check_expired_bans.start()
        if not self.check_expired_voice_mutes.is_running():
            self.check_expired_voice_mutes.start()
        if not self.check_expired_text_mutes.is_running():
            self.check_expired_text_mutes.start()
        # if not self.check_expired_stages.is_running():
        #     self.check_expired_stages.start()
        if not self.check_expired_bugs.is_running():
            self.check_expired_bugs.start()
        if not self.check_guild_owners.is_running():
            self.check_guild_owners.start()
        if not self.check_sysadmin.is_running():
            self.check_sysadmin.start()
        if not self.temporarily_cleanup_overwrites.is_running():
            self.temporarily_cleanup_overwrites.start()

    #        if not self.update_video_room_status.is_running():
    #            self.update_video_room_status.start()

    @tasks.loop(minutes=5)
    async def check_expired_bans(self):
        await self.__ban_service.clean_expired()
        logger.info("Cleaned up expired bans.")

    @tasks.loop(seconds=15)
    async def check_expired_voice_mutes(self):
        await self.__voice_mute_service.clean_expired()
        logger.info("Cleaned up expired voice-mutes.")

    # @tasks.loop(minutes=1)
    # async def check_expired_stages(self):
    #     self.__database_factory.model = Stage
    #     expired_stages = await self.__database_factory.select(expired=True)
    #     if expired_stages:
    #         for expired_stage in expired_stages:
    #             channel_snowflake = int(expired_stage.channel_snowflake)
    #             guild_snowflake = int(expired_stage.guild_snowflake)
    #             guild = self.__bot.get_guild(guild_snowflake)
    #             if guild is None:
    #                 await Stage.delete(
    #                     channel_snowflake=channel_snowflake,
    #                     guild_snowflake=guild_snowflake,
    #                 )
    #                 logger.info(
    #                     f"Unable to locate guild {guild_snowflake}, cleaning up expired stage."
    #                 )
    #                 continue
    #             channel = guild.get_channel(channel_snowflake)
    #             if channel is None:
    #                 await Stage.delete(
    #                     channel_snowflake=channel_snowflake,
    #                     guild_snowflake=guild_snowflake,
    #                 )
    #                 logger.info(
    #                     f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired voice-mute."
    #                 )
    #                 continue
    #             voice_mutes = await VoiceMute.select(
    #                 channel_snowflake=channel.id,
    #                 guild_snowflake=guild.id,
    #                 target="room",
    #             )
    #             await Stage.delete(
    #                 channel_snowflake=channel.id, guild_snowflake=guild.id
    #             )
    #             for voice_mute in voice_mutes:
    #                 member_snowflake = voice_mute.member_snowflake
    #                 member = guild.get_member(member_snowflake)
    #                 if member is None:
    #                     await VoiceMute.delete(
    #                         channel_snowflake=channel.id,
    #                         member_snowflake=int(member_snowflake),
    #                         guild_snowflake=guild.id,
    #                         target="room",
    #                     )
    #                     logger.info(
    #                         f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}) from expired stage."
    #                     )
    #                     continue
    #                 await VoiceMute.delete(
    #                     channel_snowflake=channel.id,
    #                     member_snowflake=member.id,
    #                     guild_snowflake=guild.id,
    #                     target="room",
    #                 )
    #                 if (
    #                     member.voice
    #                     and member.voice.channel
    #                     and member.voice.mute
    #                     and member.voice_channel.id == channel.id
    #                 ):
    #                     try:
    #                         await member.edit(
    #                             mute=False, reason="Stage room closed automatically."
    #                         )
    #                         logger.info(
    #                             f"Undone voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in in guild {guild.name} ({guild_snowflake}) after stage expired."
    #                         )
    #                     except discord.Forbidden as e:
    #                         logger.warning(
    #                             f"Unable to undo voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}). {str(e).capitalize()}"
    #                         )
    #                 else:
    #                     logger.info(
    #                         f"Member {member.display_name} ({member.id}) is not in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), skipping undo voice-mute."
    #                     )
    #         logger.info("Cleaned up expired stages.")

    #    @tasks.loop(minutes=5)
    #    async def update_video_room_status(self):
    #        video_rooms = await VideoRoom.fetch_all()
    #        for video_room in video_rooms:
    #            channel = self.__bot.get_channel(video_room.channel_snowflake)
    #            if channel:
    #                try:
    #                    await channel.edit(status="Video-Only Room", reason="Reset video-only room status")
    #                except discord.Forbidden as e:
    #                    logger.info("Failed to enforce video room status")
    #            else:
    #                logger.info("Failed to enforce video room status")
    @tasks.loop(hours=8)
    async def check_expired_bugs(self):
        await self.__bug_service.clean_expired()
        # await self.__developer_service.ping_devs_about_bugs()
        self.__bot.logger.info("Sent developer log to developers.")

    @tasks.loop(minutes=1)
    async def check_expired_text_mutes(self):
        self.__bot.logger.info("Cleaned up expired text-mutes.")

    @tasks.loop(hours=8)
    async def check_guild_owners(self):
        await self.__guild_owner_service.update_guild_owners()
        self.__bot.logger.info("Updated active guild owners.")

    @tasks.loop(hours=1)
    async def temporarily_cleanup_overwrites(self):
        await self.__ban_service.clean_overwrites()
        await self.__text_mute_service.clean_overwrites()
        self.__bot.logger.info("Reset ban and text-mute overwrites.")

    @tasks.loop(hours=8)
    async def check_sysadmin(self):
        await self.__sysadmin_service.update_sysadmin()
        self.__bot.logger.info("Updated active sysadmin.")

    @tasks.loop(hours=24)
    async def backup_database(self):
        try:
            db = Database(config=self.__bot.config)
            db.create_backup_directory()
            db.execute_backup()
            logger.info("Backup completed successfully.")
        except Exception as e:
            logger.error(f"Error during database backup: {str(e).capitalize()}")

    @backup_database.before_loop
    async def before_backup(self):
        await self.__bot.wait_until_ready()

    @check_expired_bans.before_loop
    async def before_check_expired_bans(self):
        await self.__bot.wait_until_ready()

    @check_expired_voice_mutes.before_loop
    async def before_check_expired_voice_mutes(self):
        await self.__bot.wait_until_ready()

    @check_expired_text_mutes.before_loop
    async def before_check_expired_text_mutes(self):
        await self.__bot.wait_until_ready()

    # @check_expired_stages.before_loop
    # async def before_check_expired_stages(self):
    #     await self.__bot.wait_until_ready()
    #
    @check_guild_owners.before_loop
    async def before_check_guild_owners(self):
        await self.__bot.wait_until_ready()

    @check_sysadmin.before_loop
    async def before_check_sysadmin(self):
        await self.__bot.wait_until_ready()

    @temporarily_cleanup_overwrites.before_loop
    async def before_temporarily_cleanup_overwrites(self):
        await self.__bot.wait_until_ready()


#    @update_video_room_status.before_loop
#    async def before_update_video_room_status(self):
#        await self.__bot.wait_until_ready()


async def setup(bot: DiscordBot):
    await bot.add_cog(ScheduledTasks(bot))
