"""scheduled_tasks.py A discord.py cog containing scheduled tasks for the Vyrtuous bot.

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

from vyrtuous.bot.discord_bot import DiscordBot

# from vyrtuous.db.rooms.video.video_room import VideoRoom
from vyrtuous.commands.fields.duration import DurationObject
from vyrtuous.db.database import Database
from vyrtuous.db.infractions.ban.ban import Ban
from vyrtuous.db.infractions.tmute.text_mute import TextMute
from vyrtuous.db.infractions.vmute.voice_mute import VoiceMute
from vyrtuous.db.mgmt.bug.bug import Bug
from vyrtuous.db.roles.dev.developer import Developer
from vyrtuous.db.roles.owner.guild_owner import GuildOwner
from vyrtuous.db.roles.sysadmin.sysadmin import Sysadmin
from vyrtuous.db.rooms.stage.stage import Stage
from vyrtuous.utils.logger import logger


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
        expired_bans = await Ban.select(expired=True)
        if expired_bans:
            for expired_ban in expired_bans:
                channel_snowflake = expired_ban.channel_snowflake
                guild_snowflake = expired_ban.guild_snowflake
                member_snowflake = expired_ban.member_snowflake
                kwargs = {
                    "channel_snowflake": int(channel_snowflake),
                    "guild_snowflake": int(guild_snowflake),
                    "member_snowflake": member_snowflake,
                }
                guild = self.bot.get_guild(guild_snowflake)
                if guild is None:
                    await Ban.delete(**kwargs)
                    logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired ban."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await Ban.delete(**kwargs)
                    logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}, cleaning up expired ban."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    await Ban.delete(**kwargs)
                    logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), cleaning up expired ban."
                    )
                    continue
                await Ban.delete(**kwargs)
                try:
                    await channel.set_permissions(
                        member, view_channel=None, reason="Cleaning up expired ban."
                    )
                except discord.Forbidden as e:
                    logger.error(str(e).capitalize())
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
                kwargs = {
                    "channel_snowflake": int(channel_snowflake),
                    "guild_snowflake": int(guild_snowflake),
                    "member_snowflake": member_snowflake,
                    "target": target,
                }
                if guild is None:
                    await VoiceMute.delete(**kwargs)
                    logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired voice-mute."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await VoiceMute.delete(**kwargs)
                    logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired voice-mute."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    await VoiceMute.delete(**kwargs)
                    logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.name}), cleaning up expired voice-mute."
                    )
                    continue
                await VoiceMute.delete(**kwargs)
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
                    await Stage.delete(
                        channel_snowflake=int(channel_snowflake),
                        guild_snowflake=int(guild_snowflake),
                    )
                    logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired stage."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await Stage.delete(
                        channel_snowflake=int(channel_snowflake),
                        guild_snowflake=int(guild_snowflake),
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
                            member_snowflake=int(member_snowflake),
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
    async def check_expired_bugs(self):
        now = datetime.now(timezone.utc)
        bugs = await Bug.select(resolved=True)
        if bugs:
            for bug in bugs:
                channel_snowflake = bug.channel_snowflake
                guild_snowflake = bug.guild_snowflake
                member_snowflakes = bug.member_snowflakes
                message_snowflake = bug.message_snowflake
                reference = bug.id
                if bug.created_at < now - timedelta(weeks=1):
                    guild = self.bot.get_guild(bug.guild_snowflake)
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
                    for member_snowflake in member_snowflakes:
                        member = channel.get_member(member_snowflake)
                        if member is None:
                            logger.info(
                                f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), not sending developer log."
                            )
                            continue
                        embed.set_thumbnail(url=member.display_avatar.url)
                        try:
                            msg = await channel.fetch_message(message_snowflake)
                        except Exception as e:
                            logger.warning(
                                f"Unable to locate a message {msg} in {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), deleting developer log. {str(e).capitalize()}"
                            )
                            return await Bug.delete(id=reference)
                        time_since_updated = await DurationObject.from_expires_at(
                            bug.updated_at
                        )
                        assigned_developer_mentions = []
                        for developer_snowflake in bug.member_snowflakes:
                            assigned_developer = self.bot.get_user(developer_snowflake)
                            assigned_developer_mentions.append(
                                assigned_developer.mention
                            )
                        embed.add_field(
                            name=f"Updated: {time_since_updated}",
                            value=f"**Link:** {msg.jump_url}\n**Developers:** {', '.join(assigned_developer_mentions)}\n**Notes:** {bug.notes}",
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
                kwargs = {
                    "channel_snowflake": int(channel_snowflake),
                    "guild_snowflake": int(guild_snowflake),
                    "member_snowflake": member_snowflake,
                }
                guild = self.bot.get_guild(guild_snowflake)
                if guild is None:
                    await TextMute.delete(**kwargs)
                    logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired text-mute."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await TextMute.delete(**kwargs)
                    logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired text-mute."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    await TextMute.delete(**kwargs)
                    logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), cleaning up expired text-mute."
                    )
                    continue
                await TextMute.delete(**kwargs)
                try:
                    await channel.set_permissions(
                        target=member,
                        send_messages=None,
                        add_reactions=None,
                        reason="Cleaning up expired text-mute",
                    )
                except discord.Forbidden as e:
                    logger.error(str(e).capitalize())
            logger.info("Cleaned up expired text-mutes.")

    @tasks.loop(hours=8)
    async def check_guild_owners(self):
        for guild in self.bot.guilds:
            guild_owner = await GuildOwner.select(
                guild_snowflake=guild.id, singular=True
            )
            if guild_owner and guild_owner.member_snowflake == guild.owner_id:
                logger.info(
                    f"Guild owner ({guild_owner.member_snowflake}) already in the db."
                )
                continue
            else:
                guild_owner = GuildOwner(
                    guild_snowflake=guild.id, member_snowflake=guild.owner_id
                )
                await guild_owner.create()
                logger.info(
                    f"Guild owner ({guild_owner.member_snowflake}) added to the db."
                )

    @tasks.loop(hours=1)
    async def temporarily_cleanup_overwrites(self):
        now = datetime.now(timezone.utc)
        text_mutes = await TextMute.select()
        for text_mute in text_mutes:
            channel_snowflake = text_mute.channel_snowflake
            guild_snowflake = text_mute.guild_snowflake
            member_snowflake = text_mute.member_snowflake
            where_kwargs = {
                "channel_snowflake": int(channel_snowflake),
                "guild_snowflake": int(guild_snowflake),
                "member_snowflake": member_snowflake,
            }
            set_kwargs = {"reset": True}
            if not text_mute.reset and text_mute.last_muted < now - timedelta(weeks=1):
                guild = self.bot.get_guild(guild_snowflake)
                if guild is None:
                    logger.info(
                        f"Unable to locate guild {guild_snowflake} for removing overwrite."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}) for removing overwrite."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}) for removing overwrite."
                    )
                    continue
                try:
                    await channel.set_permissions(
                        target=member,
                        overwrite=None,
                        reason="Resetting text-mute overwrite",
                    )
                except discord.Forbidden as e:
                    logger.error(str(e).capitalize())
                await TextMute.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        bans = await Ban.select()
        for ban in bans:
            channel_snowflake = ban.channel_snowflake
            guild_snowflake = ban.guild_snowflake
            member_snowflake = ban.member_snowflake
            where_kwargs = {
                "channel_snowflake": int(channel_snowflake),
                "guild_snowflake": int(guild_snowflake),
                "member_snowflake": member_snowflake,
            }
            set_kwargs = {"reset": True}
            if not ban.reset and ban.last_kicked < now - timedelta(weeks=1):
                guild = self.bot.get_guild(guild_snowflake)
                if guild is None:
                    logger.info(
                        f"Unable to locate guild {guild_snowflake} for removing overwrite."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}) for removing overwrite."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}) for removing overwrite."
                    )
                    continue
                try:
                    await channel.set_permissions(
                        target=member, overwrite=None, reason="Resetting ban overwrite."
                    )
                except discord.Forbidden as e:
                    logger.error(str(e).capitalize())
                await Ban.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        logger.info("Reset ban and text-mute overwrites.")

    @tasks.loop(hours=8)
    async def check_sysadmin(self):
        member_snowflake = self.bot.config.get("discord_owner_id", None)
        sysadmin = await Sysadmin.select(
            member_snowflake=int(member_snowflake), singular=True
        )
        if not sysadmin:
            sysadmin = Sysadmin(member_snowflake=int(member_snowflake))
            await sysadmin.create()
            logger.info(f"Sysadmin ({member_snowflake}) added to the db.")
        else:
            logger.info(f"Sysadmin ({member_snowflake}) already in the db.")

    @tasks.loop(hours=24)
    async def backup_database(self):
        try:
            db = Database(config=self.bot.config)
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

    @check_guild_owners.before_loop
    async def before_check_guild_owners(self):
        await self.bot.wait_until_ready()

    @check_sysadmin.before_loop
    async def before_check_sysadmin(self):
        await self.bot.wait_until_ready()


#    @update_video_room_status.before_loop
#    async def before_update_video_room_status(self):
#        await self.bot.wait_until_ready()


async def setup(bot: DiscordBot):
    await bot.add_cog(ScheduledTasks(bot))
