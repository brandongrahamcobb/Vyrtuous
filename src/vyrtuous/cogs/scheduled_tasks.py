"""scheduled_tasks.py A discord.py cog containing scheduled tasks for the Vyrtuous bot.

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
"""

from datetime import datetime, timedelta, timezone
from discord.ext import commands, tasks
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.actions.ban import Ban
from vyrtuous.database.database import Database
from vyrtuous.database.roles.developer import Developer
from vyrtuous.utils.properties.duration import DurationObject
from vyrtuous.database.logs.developer_log import DeveloperLog
from vyrtuous.utils.setup_logging import logger
from vyrtuous.database.actions.text_mute import TextMute
from vyrtuous.database.rooms.stage import Stage
from vyrtuous.database.rooms.video_room import VideoRoom
from vyrtuous.database.actions.voice_mute import VoiceMute

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
        if not self.check_expired_developer_logs.is_running():
            self.check_expired_developer_logs.start()

    #        if not self.update_video_room_status.is_running():
    #            self.update_video_room_status.start()

    @tasks.loop(minutes=5)
    async def check_expired_bans(self):
        expired_bans = await Ban.select(expired=True)
        if expired_bans:
            for expired_ban in expired_bans:
                member_snowflake = expired_ban.member_snowflake
                guild_snowflake = expired_ban.guild_snowflake
                channel_snowflake = expired_ban.channel_snowflake
                guild = self.bot.get_guild(guild_snowflake)
                if guild is None:
                    logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired ban."
                    )
                    await Ban.delete(guild_snowflake=guild_snowflake)
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}, cleaning up expired ban."
                    )
                    await Ban.delete(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        member_snowflake=member_snowflake,
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), cleaning up expired ban."
                    )
                    await Ban.delete(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        member_snowflake=member_snowflake,
                    )
                    continue
                await Ban.delete(
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                    member_snowflake=member_snowflake,
                )
                try:
                    await channel.set_permissions(member, overwrite=None)
                    logger.info(
                        f"Unbanned member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake})."
                    )
                except discord.Forbidden as e:
                    logger.warning(
                        f"Unable to unban member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}). {str(e).capitalize()}"
                    )
            logger.info("Cleaned up expired bans.")

    @tasks.loop(seconds=15)
    async def check_expired_voice_mutes(self):
        expired_voice_mutes = await VoiceMute.select(expired=True)
        if expired_voice_mutes:
            for expired_voice_mute in expired_voice_mutes:
                guild_snowflake = expired_voice_mute.guild_snowflake
                member_snowflake = expired_voice_mute.member_snowflake
                channel_snowflake = expired_voice_mute.channel_snowflake
                target = expired_voice_mute.target
                guild = self.bot.get_guild(guild_snowflake)
                if guild is None:
                    await VoiceMute.delete(
                        guild_snowflake=guild_snowflake, target=target
                    )
                    logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired voice-mute."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await VoiceMute.delete(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        target=target,
                    )
                    logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired voice-mute."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    await VoiceMute.delete(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        member_snowflake=member_snowflake,
                        target="user",
                    )
                    logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.name}), cleaning up expired voice-mute."
                    )
                    continue
                await VoiceMute.delete(
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                    member_snowflake=member.id,
                    target="user",
                )
                if (
                    member.voice
                    and member.voice.channel
                    and member.voice.channel.id == channel_snowflake
                ):
                    try:
                        await member.edit(mute=False)
                        logger.info(
                            f"Undone voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake})."
                        )
                    except discord.Forbidden as e:
                        logger.warning(
                            f"Unable to undo voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}). {str(e).capitalize()}"
                        )
                else:
                    logger.info(
                        f"Member {member.display_name} ({member.id}) is not in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), skipping undo voice-mute."
                    )
            logger.info("Cleaned up expired voice-mutes.")

    @tasks.loop(minutes=1)
    async def check_expired_stages(self):
        expired_stages = await Stage.select(expired=True)
        if expired_stages:
            for expired_stage in expired_stages:
                channel_snowflake = expired_stage.channel_snowflake
                guild_snowflake = expired_stage.guild_snowflake
                guild = self.bot.get_guild(guild_snowflake)
                if guild is None:
                    await Stage.delete(guild_snowflake=guild_snowflake)
                    logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired stage."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await Stage.delete(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                    )
                    logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired voice-mute."
                    )
                    continue
                voice_mutes = await VoiceMute.select(
                    channel_snowflake=channel.id,
                    guild_snowflake=guild.id,
                    target="room",
                )
                await Stage.delete(
                    channel_snowflake=channel.id, guild_snowflake=guild.id
                )
                for voice_mute in voice_mutes:
                    member_snowflake = voice_mute.member_snowflake
                    member = guild.get_member(member_snowflake)
                    if member is None:
                        await VoiceMute.delete(
                            channel_snowflake=channel.id,
                            member_snowflake=member_snowflake,
                            guild_snowflake=guild.id,
                            target="room",
                        )
                        logger.info(
                            f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}) from expired stage."
                        )
                        continue
                    await VoiceMute.delete(
                        channel_snowflake=channel.id,
                        member_snowflake=member.id,
                        guild_snowflake=guild.id,
                        target="room",
                    )
                    if (
                        member.voice
                        and member.voice.channel
                        and member.voice.mute
                        and member.voice_channel.id == channel.id
                    ):
                        try:
                            await member.edit(
                                mute=False, reason="Stage room closed automatically."
                            )
                            logger.info(
                                f"Undone voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in in guild {guild.name} ({guild_snowflake}) after stage expired."
                            )
                        except discord.Forbidden as e:
                            logger.warning(
                                f"Unable to undo voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}). {str(e).capitalize()}"
                            )
                    else:
                        logger.info(
                            f"Member {member.display_name} ({member.id}) is not in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), skipping undo voice-mute."
                        )
            logger.info("Cleaned up expired stages.")

    #    @tasks.loop(minutes=5)
    #    async def update_video_room_status(self):
    #        video_rooms = await VideoRoom.fetch_all()
    #        for video_room in video_rooms:
    #            channel = self.bot.get_channel(video_room.channel_snowflake)
    #            if channel:
    #                try:
    #                    await channel.edit(status="Video-Only Room", reason="Reset video-only room status")
    #                except discord.Forbidden as e:
    #                    logger.info("Failed to enforce video room status")
    #            else:
    #                logger.info("Failed to enforce video room status")
    @tasks.loop(hours=8)
    async def check_expired_developer_logs(self):
        now = datetime.now(timezone.utc)
        developer_logs = await DeveloperLog.select(resolved=True)
        if developer_logs:
            for developer_log in developer_logs:
                channel_snowflake = developer_log.channel_snowflake
                guild_snowflake = developer_log.guild_snowflake
                member_snowflake = developer_log.member_snowflake
                if developer_log.created_at < now - timedelta(weeks=1):
                    guild = self.bot.get_guild(developer_log.guild_snowflake)
                    if guild is None:
                        logger.info(
                            f"Unable to locate guild {guild_snowflake}, not sending developer log."
                        )
                        continue
                    embed = discord.Embed(
                        title=f"\U000026a0\U0000fe0f An issue is unresolved in {guild.name}",
                        color=discord.Color.red(),
                    )
                    channel = guild.get_channel(channel_snowflake)
                    if channel is None:
                        logger.info(
                            f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}, not sending developer log."
                        )
                        continue
                    member = channel.get_member(member_snowflake)
                    if member is None:
                        logger.info(
                            f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), not sending developer log."
                        )
                        continue
                    embed.set_thumbnail(url=member.display_avatar.url)
                    try:
                        msg = await channel.fetch_message(developer_log.message)
                    except Exception as e:
                        logger.warning(
                            f"Unable to locate a message {msg} in {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), not sending developer log. {str(e).capitalize()}"
                        )
                    time_since_updated = await DurationObject.from_expires_at(
                        developer_log.updated_at
                    )
                    assigned_developer_mentions = []
                    for developer_snowflake in developer_log.developer_snowflakes:
                        assigned_developer = self.bot.get_user(developer_snowflake)
                        assigned_developer_mentions.append(assigned_developer.mention)
                    embed.add_field(
                        name=f"Updated: {time_since_updated}",
                        value=f"**Link:** {msg.jump_url}\n**Developers:** {', '.join(assigned_developer_mentions)}\n**Notes:** {developer_log.notes}",
                        inline=False,
                    )
                    developers = await Developer.select()
                    for developer in developers:
                        member = guild.get_member(developer.member_snowflake)
                        if member is None:
                            logger.info(
                                f"Unable to locate member {member.id} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), not sending developer log."
                            )
                            continue
                        try:
                            await member.send(embed=embed)
                            logger.info(
                                f"Sent the issue to member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake})."
                            )
                        except Exception as e:
                            logger.warning(
                                f"Unable to send the issue to member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}). {str(e).capitalize()}"
                            )
            logger.info("Sent developer log to developers.")

    @tasks.loop(minutes=1)
    async def check_expired_text_mutes(self):
        expired_text_mutes = await TextMute.select(expired=True)
        if expired_text_mutes:
            for expired_text_mute in expired_text_mutes:
                channel_snowflake = expired_text_mute.channel_snowflake
                guild_snowflake = expired_text_mute.guild_snowflake
                member_snowflake = expired_text_mute.member_snowflake
                guild = self.bot.get_guild(guild_snowflake)
                if guild is None:
                    await TextMute.delete(guild_snowflake=guild_snowflake)
                    logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired text-mute."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await TextMute.delete(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                    )
                    logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired text-mute."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), cleaning up expired text-mute"
                    )
                    await TextMute.delete(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        member_snowflake=member_snowflake,
                    )
                    continue
                await TextMute.delete(
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                    member_snowflake=member_snowflake,
                )
                try:
                    await channel.set_permissions(member, send_messages=None)
                    logger.info(
                        f"Undone text-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake})."
                    )
                except discord.Forbidden as e:
                    logger.warning(
                        f"Unable to undo text-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}). {str(e).capitalize()}"
                    )
            logger.info("Cleaned up expired text-mutes.")

    @tasks.loop(hours=24)
    async def backup_database(self):
        try:
            db = Database()
            db.create_backup_directory()
            db.execute_backup()
            logger.info("Backup completed successfully.")
        except Exception as e:
            logger.error(f"Error during database backup: {str(e).capitalize()}")

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


#    @update_video_room_status.before_loop
#    async def before_update_video_room_status(self):
#        await self.bot.wait_until_ready()


async def setup(bot: DiscordBot):
    await bot.add_cog(ScheduledTasks(bot))
