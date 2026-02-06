"""!/bin/python3
channel_event_listeners.py A discord.py cog containing channel event listeners for the Vyrtuous bot.

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

import time
from collections import defaultdict
from datetime import datetime, timedelta, timezone

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.fields.duration import DurationObject
from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.alias.alias import Alias
from vyrtuous.db.infractions.ban.ban import Ban
from vyrtuous.db.infractions.ban.ban_service import BanService
from vyrtuous.db.infractions.flag.flag import Flag
from vyrtuous.db.infractions.smute.server_mute import ServerMute
from vyrtuous.db.infractions.tmute.text_mute import TextMute
from vyrtuous.db.infractions.vmute.voice_mute import VoiceMute
from vyrtuous.db.mgmt.cap.cap import Cap
from vyrtuous.db.mgmt.stream.stream_service import StreamService
from vyrtuous.db.roles.coord.coordinator import Coordinator
from vyrtuous.db.roles.mod.moderator import Moderator
from vyrtuous.db.rooms.stage.stage import Stage
from vyrtuous.db.rooms.temp.temporary_room import TemporaryRoom
from vyrtuous.db.rooms.video.video_room import VideoRoom
from vyrtuous.db.rooms.video.video_room_service import VideoRoomService
from vyrtuous.utils.logger import logger


class ChannelEventListeners(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.flags = []
        self.join_log = defaultdict(list)
        self._ready_done = False
        self.deleted_rooms = {}

    async def cog_load(self):
        VideoRoomService.video_rooms = await VideoRoom.select()
        self.flags = await Flag.select()

    @commands.Cog.listener()
    async def on_guild_channel_grant(self, channel: discord.abc.GuildChannel):
        guild = channel.guild
        name = channel.name
        for c in guild.channels:
            if c.id != channel.id and c.name == name:
                return
        room = self.deleted_rooms.pop(name, None)
        if not room:
            room = await TemporaryRoom.select(
                guild_snowflake=guild.id, room_name=name, singular=True
            )
        if room:
            old_id = room.channel_snowflake
        set_kwargs = {"channel_snowflake": channel.id}
        where_kwargs = {
            "channel_snowflake": old_id,
            "guild_snowflake": channel.guild.id,
        }
        await Alias.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        await Ban.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        await Cap.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        await Coordinator.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        await Moderator.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        await Stage.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        await TemporaryRoom.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        await TextMute.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        await VoiceMute.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)

    @commands.Cog.listener()
    async def on_guild_channel_delete(self, channel: discord.abc.GuildChannel):
        room = await TemporaryRoom.select(
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
            singular=True,
        )
        if room:
            self.deleted_rooms[channel.name] = room

    @commands.Cog.listener()
    async def on_guild_channel_update(self, before, after):
        if before.name == after.name:
            return
        set_kwargs = {"room_name": after.id}
        where_kwargs = {"room_name": before.id}
        await TemporaryRoom.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)

    @commands.Cog.listener()
    async def on_voice_state_update(
        self,
        member: discord.Member,
        before: discord.VoiceState,
        after: discord.VoiceState,
    ):
        if member.bot:
            return
        await BanService.ban_overwrite(channel=after.channel, member=member)
        if before.channel == after.channel:
            if before.mute == after.mute:
                if before.self_mute == after.self_mute:
                    return
        await VideoRoomService.reinforce_video_room(
            member=member, before=before, after=after
        )
        target = "user"
        server_mute = await ServerMute.select(member_snowflake=member.id, singular=True)
        if server_mute:
            if member.guild.id == server_mute.guild_snowflake:
                if not after.mute:
                    try:
                        await member.edit(mute=True, reason="Server mute is active.")
                    except discord.Forbidden as e:
                        logger.warning(
                            f"No permission to "
                            f"edit mute for {member.display_name}. {str(e).capitalize()}"
                        )
                    except discord.HTTPException as e:
                        logger.warning(
                            f"Failed to edit mute for "
                            f"{member.display_name}: "
                            f"{str(e).capitalize()}"
                        )
                return
        if after.channel:
            should_be_muted = False
            stage = await Stage.select(channel_snowlfake=after.channel)
            if stage:
                Stage.send_stage_ask_to_speak_message(
                    join_log=self.join_log, member=member
                )
                default_kwargs = {
                    "channel_snowflake": after.channel.id,
                    "guild_snowflake": after.channel.guild.id,
                    "member_snowflake": member.id,
                }
                updated_kwargs = default_kwargs.copy()
                highest_role = await PermissionService.resolve_highest_role(
                    updated_kwargs=updated_kwargs
                )
                if highest_role == "Everyone":
                    should_be_muted = True
                    target = "room"
                    expires_in = stage.expires_in
            voice_mute = await VoiceMute.select(
                channel_snowflake=after.channel.id,
                guild_snowflake=after.channel.guild.id,
                member_snowflake=member.id,
                target="user",
                singular=True,
            )
            if voice_mute:
                should_be_muted = True
            if not before.mute and before.channel == after.channel and after.mute:
                if member.id in PermissionService.get_invincible_members():
                    embed = discord.Embed(
                        title=f"\u1f4aB {member.display_name} is a hero!",
                        description=f"{member.display_name} cannot be muted.",
                        color=discord.Color.gold(),
                    )
                    embed.set_thumbnail(url=member.display_avatar.url)
                    await after.channel.send(embed=embed)
                elif not should_be_muted:
                    expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
                    voice_mute = VoiceMute(
                        channel_snowflake=after.channel.id,
                        expires_in=expires_in,
                        guild_snowflake=after.channel.guild.id,
                        member_snowflake=member.id,
                        reason="No reason provided.",
                        target=target,
                    )
                    await voice_mute.create()
                    should_be_muted = True
                    duration = DurationObject.from_expires_in(expires_in=expires_in)
                    await StreamService.send_entry(
                        channel_snowflake=after.channel.id,
                        duration=duration,
                        identifier="vmute",
                        is_channel_scope=True,
                        member=member,
                        reason="Right-click voice-mute.",
                    )
            elif before.mute and not after.mute and before.channel == after.channel:
                await VoiceMute.delete(
                    channel_snowflake=before.channel.id,
                    guild_snowflake=before.channel.guild.id,
                    member_snowflake=member.id,
                    target=target,
                )
                should_be_muted = False
                duration = DurationObject("0")
                await StreamService.send_entry(
                    channel_snowflake=after.channel.id,
                    duration=duration,
                    identifier="unvmute",
                    is_channel_scope=True,
                    member=member,
                    reason="Right-click unvoice-mute.",
                )
            if after.mute != should_be_muted:
                try:
                    await member.edit(
                        mute=should_be_muted,
                        reason=f"Setting mute to {should_be_muted} "
                        f"in {after.channel.name}",
                    )
                except discord.Forbidden as e:
                    logger.warning(
                        f"No permission to edit "
                        f"mute for {member.display_name}. {str(e).capitalize()}"
                    )
                except discord.HTTPException as e:
                    logger.warning(
                        f"Failed to edit mute for "
                        f"{member.display_name}: {str(e).capitalize()}"
                    )
            await self.print_flags(member, after.channel)

    #                    explicit_deny_roles = []
    #                    for role in member.roles:
    #                        ow = after_channel.overwrites_for(role)
    #                        if ow.speak is False:
    #                            explicit_deny_roles.append(role)
    #                    if explicit_deny_roles:
    #                        try:
    #                            await member.move_to(after_channel)
    #                            await logger.warning(
    #                                f"🔇 Auto-muted {member.mention} in **{after_channel.name}** "
    #                                f"due to explicit speak deny from roles: "
    #                                f"{', '.join(r.name for r in explicit_deny_roles)}"
    #                            )
    #                        except Exception as e:
    #                            logger.warning(
    #                                f"⚠ Failed to auto-mute {member.mention} in "
    #                                f"**{after_channel.name}** — `{str(e).capitalize()}`"
    #                            )
    #                        except discord.HTTPException as e:
    #                            logger.warning(f'Failed to mute {member.display_name}: {str(e).capitalize()}')

    async def print_flags(
        self, member: discord.Member, after_channel: discord.abc.GuildChannel
    ):
        if after_channel.id == 1222056499959042108:
            for flag in self.flags:
                if flag.channel_snowflake == after_channel.id:
                    if flag.member_snowflake == member.id:
                        embed = discord.Embed(
                            title=f"\u26a0\ufe0f {member.display_name} " f"is flagged",
                            color=discord.Color.red(),
                        )
                        embed.set_thumbnail(url=member.display_avatar.url)
                        embed.add_field(
                            name=f"Channel: {after_channel.mention}",
                            value=f"Reason: {flag.reason}",
                            inline=False,
                        )
                        now = time.time()
                        self.join_log[member.id] = [
                            t for t in self.join_log[member.id] if now - t < 300
                        ]
                        if len(self.join_log[member.id]) < 1:
                            self.join_log[member.id].append(now)
                            await after_channel.send(embed=embed)


async def setup(bot: DiscordBot):
    await bot.add_cog(ChannelEventListeners(bot))
