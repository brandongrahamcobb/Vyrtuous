"""event_listeners.py A discord.py cog containing event listeners for the Vyrtuous bot.

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

from collections import defaultdict
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace
import asyncio
import time

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.actions.ban import Ban
from vyrtuous.db.actions.flag import Flag
from vyrtuous.db.actions.server_mute import ServerMute
from vyrtuous.db.actions.text_mute import TextMute
from vyrtuous.db.actions.voice_mute import VoiceMute
from vyrtuous.db.roles.administrator import Administrator, AdministratorRole
from vyrtuous.db.roles.coordinator import Coordinator
from vyrtuous.db.roles.guild_owner import GuildOwner
from vyrtuous.db.roles.moderator import Moderator
from vyrtuous.db.rooms.stage import Stage
from vyrtuous.db.rooms.temporary_room import TemporaryRoom
from vyrtuous.db.rooms.video_room import VideoRoom
from vyrtuous.db.mgmt.cap import Cap
from vyrtuous.db.mgmt.stream import Streaming
from vyrtuous.fields.duration import DurationObject
from vyrtuous.service.discord_object_service import DiscordObjectNotFound

from vyrtuous.utils.check import has_equal_or_lower_role_wrapper
from vyrtuous.utils.logger import logger
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.invincibility import Invincibility


class EventListeners(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.flags = []
        self.message_service = MessageService(self.bot)
        self.join_log = defaultdict(list)
        self._ready_done = False
        self.deleted_rooms = {}

    async def cog_load(self):
        VideoRoom.video_rooms = await VideoRoom.select()
        for room in VideoRoom.video_rooms:
            channel = self.bot.get_channel(room.channel_snowflake)
            if channel:
                try:
                    await channel.edit(
                        status="Video-Only Room",
                        reason="Enforce default video-only status",
                    )
                except discord.Forbidden as e:
                    logger.warning(str(e).capitalize())
        self.flags = await Flag.select()

    @commands.Cog.listener()
    async def on_guild_update(self, before: discord.Guild, after: discord.Guild):
        if before.owner_id != after.owner_id:
            where_kwargs = {
                "guild_snowflake": before.guild.id,
                "member_snowflake": before.owner_id,
            }
            set_kwargs = {
                "guild_snowflake": after.guild.id,
                "member_snowflake": after.owner_id,
            }
            await GuildOwner.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)

    @commands.Cog.listener()
    async def on_guild_channel_grant(self, channel: discord.abc.GuildChannel):
        guild = channel.guild
        name = channel.name
        for c in guild.channels:
            if c.id != channel.id and c.name == name:
                return
        room = self.deleted_rooms.pop(name, None)
        if not room:
            room = await TemporaryRoom.select(guild_snowflake=guild.id, room_name=name)
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
            channel_snowflake=channel.id, guild_snowflake=channel.guild.id
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

    # Done
    @commands.Cog.listener()
    async def on_voice_state_update(
        self,
        member: discord.Member,
        before: discord.VoiceState,
        after: discord.VoiceState,
    ):
        if member.id == self.bot.user.id:
            return
        if not after.channel:
            VideoRoom.cancel_task((member.guild.id, member.id))
            return
        for video_room in VideoRoom.video_rooms:
            if after.channel.id != video_room.channel_snowflake:
                continue
            if not after.self_video:
                if after.channel != before.channel:
                    if after.channel.permissions_for(
                        after.channel.guild.me
                    ).send_messages:
                        await VideoRoom.enforce_video_message(
                            channel_snowflake=after.channel.id,
                            member_snowflake=member.id,
                            message=f"{get_random_emoji()} "
                            f"Hi {member.mention}, "
                            f"{after.channel.mention} is a video "
                            f"only room. You have 5 minutes to turn "
                            f"on your camera!",
                        )
            key = (member.guild.id, member.id)
            if before.channel != after.channel:
                VideoRoom.cancel_task(key)
                if not after.self_video:
                    task = asyncio.create_task(
                        VideoRoom.enforce_video(member, after.channel, 300)
                    )
                    VideoRoom.video_tasks[key] = task
                break
            if before.self_video and not after.self_video:
                VideoRoom.cancel_task(key)
                task = asyncio.create_task(
                    VideoRoom.enforce_video(member, after.channel, 60)
                )
                VideoRoom.video_tasks[key] = task
                break
            if not before.self_video and after.self_video:
                VideoRoom.cancel_task(key)
                break
        allowed = True
        if before.channel == after.channel:
            if before.mute == after.mute:
                if before.self_mute == after.self_mute:
                    allowed = False
        if member.bot:
            allowed = False
        if not allowed:
            return
        # member_role = await role_check_with_specifics(after.channel, member)

        target = "user"
        # if after.channel:
        # stage = await Stage.fetch_stage_by_channel(after.channel)
        # temporary_stage_coordinator_ids = await stage.fetch_coordinator_temporary_stage_coordinator_ids(member, after.channel)
        # if stage:
        #     target = 'room'
        #     stage.send_stage_ask_to_speak_message(join_log=self.join_log, member=member)
        # else:
        #     target = 'user'
        # if stage and (member.id not in temporary_stage_coordinator_ids) and (member_role in ('Moderator', 'Everyone')) and (before.channel != after.channel):
        #      expires_in = stage.expires_in
        #      await conn.execute('''
        #          INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_in, target, room_name)
        #          VALUES ($1, $2, $3, $4, 'room', $5)
        #          ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
        #          DO UPDATE SET expires_in = EXCLUDED.expires_in
        #      ''', member.guild.id, member.id, after.channel.id, expires_in, after.channel.name)
        server_mute = await ServerMute.select(member_snowflake=member.id)
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
            voice_mute = await VoiceMute.select(
                channel_snowflake=after.channel.id,
                guild_snowflake=after.channel.guild.id,
                member_snowflake=member.id,
                target="user",
            )
            if voice_mute:
                should_be_muted = True
            if not before.mute and before.channel == after.channel and after.mute:
                if member.id in Invincibility.get_invincible_members():
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
                    alias = SimpleNamespace(alias_type="voice_mute")
                    duration = DurationObject("1h")
                    await Streaming.send_entry(
                        alias=alias,
                        channel=after.channel,
                        duration=duration,
                        executor_role="Role-specfic",
                        is_channel_scope=True,
                        is_modification=False,
                        member=member,
                        message=None,
                        reason="Right-click voice-mute.",
                    )
            elif before.mute and not after.mute and before.channel == after.channel:
                ban = await Ban.select(
                    channel_snowflake=before.channel.id,
                    guild_snowflake=before.channel.guild.id,
                    member_snowflake=member.id,
                    singular=True,
                )
                if not ban:
                    await VoiceMute.delete(
                        channel_snowflake=before.channel.id,
                        guild_snowflake=before.channel.guild.id,
                        member_snowflake=member.id,
                        target=target,
                    )
                    should_be_muted = False
                    alias = SimpleNamespace(alias_type="unvmute")
                    duration = DurationObject("0")
                    await Streaming.send_entry(
                        alias=alias,
                        channel=after.channel,
                        duration=duration,
                        executor_role="Role-specfic",
                        is_channel_scope=True,
                        is_modification=False,
                        member=member,
                        message=None,
                        reason="Right-click voice-mute.",
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

    @commands.Cog.listener()
    async def on_member_join(self, member: discord.Member):
        user_id = member.id
        guild = member.guild
        bans = await Ban.select(guild_snowflake=guild.id, member_snowflake=user_id)
        text_mutes = await TextMute.select(
            guild_snowflake=guild.id, member_snowflake=user_id
        )
        if bans:
            for ban in bans:
                channel = guild.get_channel(ban.channel_snowflake)
                role = guild.get_role(ban.role_snowflake)
                try:
                    await member.add_roles(role, reason="Reinstating channel ban")
                except discord.Forbidden as e:
                    logger.warning(
                        f"Unable to ban member "
                        f"{member.display_name} ({member.id}) "
                        f"in channel {channel.name} ({channel.id}) "
                        f"in guild {guild.name} ({guild.id}). "
                        f"{str(e).capitalize()}"
                    )
        if text_mutes:
            for text_mute in text_mutes:
                channel = guild.get_channel(text_mute.channel_snowflake)
                role = guild.get_role(text_mute.role_snowflake)
                try:
                    await member.add_roles(role, reason="Reinstating text-mute")
                except discord.Forbidden as e:
                    logger.warning(
                        f"Unable to text-mute member "
                        f"{member.display_name} ({member.id}) "
                        f"in channel {channel.name} ({channel.id}) "
                        f"in guild {guild.name} ({guild.id}). "
                        f"{str(e).capitalize()}"
                    )

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if after.author.bot:
            return
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)
            else:
                await self.on_message(after)

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message):
        args = (
            message.content[len(self.config["discord_command_prefix"]) :]
            .strip()
            .split()
        )
        if (
            not message.guild
            or not args
            or not message.content.startswith(self.config["discord_command_prefix"])
            or (self.config["release_mode"] and message.author.id == self.bot.user.id)
        ):
            return
        alias = await Alias.select(
            alias_name=args[0],
            guild_snowflake=message.guild.id,
            singular=True,
        )
        if not alias:
            return
        state = StateService(source=message)
        channel_obj = message.guild.get_channel(alias.channel_snowflake)
        member_obj = message.guild.get_member(int(args[1]))
        executor_role = await has_equal_or_lower_role_wrapper(
            source=message,
            member_snowflake=member_obj.id,
            sender_snowflake=message.author.id,
        )
        action_duration = (
            DurationObject(args[2]) if len(args) > 2 else DurationObject("8h")
        )
        action_reason = " ".join(args[3:]) if len(args) > 3 else "No reason provided."

        alias_class = alias.alias_class
        kwargs = {}
        primary_keys = await alias_class.primary_keys()
        if "channel_snowflake" in primary_keys:
            kwargs.update({"channel_snowflake": channel_obj.id})
        if "guild_snowflake" in primary_keys:
            kwargs.update({"guild_snowflake": message.guild.id})
        if "member_snowflake" in primary_keys:
            kwargs.update({"member_snowflake": member_obj.id})
        action_existing = await alias_class.select(
            **kwargs,
            singular=True,
        )
        action_modification = False

        if action_existing:
            action_modification = True
            await alias_class.delete(
                channel_snowflake=channel_obj.id,
                guild_snowflake=message.guild.id,
                member_snowflake=member_obj.id,
            )
            alias.alias_type = alias_class.UNDO

        action_channel_cap = await Alias.generate_cap_duration(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            moderation_type=alias_class.ACT,
        )
        action_information = {
            "alias_class": alias_class,
            "action_channel_cap": action_channel_cap,
            "action_channel_snowflake": channel_obj.id,
            "action_duration": action_duration,
            "action_executor_role": executor_role,
            "action_existing": action_existing,
            "action_guild_snowflake": message.guild.id,
            "action_member_snowflake": member_obj.id,
            "action_modification": action_modification,
            "action_reason": action_reason,
            "action_role_snowflake": (
                alias.role_snowflake if alias.role_snowflake else None
            ),
        }
        expires_in_timedelta = action_duration.to_timedelta()

        action_expires_in = datetime.now(timezone.utc) + expires_in_timedelta
        action_information["action_expires_in"] = action_expires_in

        await alias.handlers[alias.alias_type](
            alias=alias,
            action_information=action_information,
            channel=channel_obj,
            member=member_obj,
            message=message,
            state=state,
        )

    @commands.Cog.listener()
    async def on_command_error(self, ctx, error):
        state = StateService(source=ctx)
        logger.error(str(error))
        if isinstance(error, commands.BadArgument):
            return await state.end(error=str(error))
        elif isinstance(error, commands.CheckFailure):
            return await state.end(error=str(error))
        elif isinstance(error, commands.MissingRequiredArgument):
            missing = error.param.name
            return await state.end(error=f"Missing required argument: `{missing}`")

    @commands.Cog.listener()
    async def on_app_command_error(self, interaction, error):
        state = StateService(source=interaction)
        logger.error(str(error))
        if isinstance(error, app_commands.BadArgument):
            return await state.end(error=str(error))
        elif isinstance(error, app_commands.CheckFailure):
            return await state.end(error=str(error))
        elif isinstance(error, DiscordObjectNotFound()):
            return await state.end(error=str(error))

    #    @commands.Cog.listener()
    #    async def on_command(self, ctx):
    #        await ctx.send("Bot is currently down. Changes will not be saved permanently.")

    @commands.Cog.listener()
    async def on_ready(self):
        if getattr(self, "_ready_done", False):
            return
        self._ready_done = True
        method_names = [cmd.callback.__name__ for cmd in self.bot.commands]
        logger.info(method_names)

    @commands.Cog.listener()
    async def on_member_update(self, before: discord.Member, after: discord.Member):
        if before.roles == after.roles:
            return
        before_role_snowflakes = {r.id for r in before.roles}
        after_role_snowflakes = {r.id for r in after.roles}
        added_roles = after_role_snowflakes - before_role_snowflakes
        removed_roles = before_role_snowflakes - after_role_snowflakes
        administrator_role_snowflakes = []
        administrator_roles = await AdministratorRole.select(
            guild_snowflake=after.guild.id
        )
        for administrator_role in administrator_roles:
            administrator_role_snowflakes.append(administrator_role.role_snowflake)
        relevant_added_roles = added_roles & set(administrator_role_snowflakes)
        relevant_removed_roles = removed_roles & set(administrator_role_snowflakes)
        if not relevant_added_roles and not relevant_removed_roles:
            return
        administrator = await Administrator.select(
            guild_snowflake=after.guild.id, member_snowflake=after.id
        )
        if not administrator and relevant_added_roles:
            administrator = Administrator(
                guild_snowflake=after.guild.id,
                member_snowflake=after.id,
                role_snowflakes=list(after_role_snowflakes),
            )
            await administrator.create()
            return
        if administrator:
            where_kwargs = {
                "guild_snowflake": after.guild.id,
                "member_snowflake": after.id,
            }
            for role_snowflake in relevant_added_roles:
                if role_snowflake not in administrator.role_snowflakes:
                    set_kwargs = {
                        "role_snowflakes": administrator.role_snowflakes.append(
                            role_snowflake
                        )
                    }
                    await Administrator.update(
                        set_kwargs=set_kwargs, where_kwargs=where_kwargs
                    )
            for role_snowflake in relevant_removed_roles:
                if role_snowflake in administrator.role_snowflakes:
                    set_kwargs = {
                        "role_snowflakes": administrator.role_snowflakes.remove(
                            role_snowflake
                        )
                    }
                    await Administrator.update(
                        set_kwargs=set_kwargs, where_kwargs=where_kwargs
                    )
            remaining_admin_roles = (
                set(administrator.role_snowflakes) & after_role_snowflakes
            )
            if not remaining_admin_roles:
                await Administrator.delete(
                    guild_snowflake=after.guild.id, member_snowflake=after.id
                )

    @commands.Cog.listener()
    async def on_guild_role_delete(self, role: discord.Role):
        administrators = await Administrator.select(role_snowflakes=[role.id])
        for administrator in administrators:
            role_snowflakes = set(administrator.role_snowflake)
            guild_snowflakes = set(administrator.guild_snowflake)
            role_snowflakes.discard(role.id)
            if not role_snowflakes:
                guild_snowflakes.discard(role.guild)
            Administrator.update(
                guild_snowflake=list(guild_snowflakes),
                member_snowflake=administrator.member_snowflake,
                role_snowflakes=list(role_snowflakes),
            )

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
    await bot.add_cog(EventListeners(bot))
