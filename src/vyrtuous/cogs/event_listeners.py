''' event_listeners.py A discord.py cog containing event listeners for the Vyrtuous bot.

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
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from discord.ext import commands
from types import SimpleNamespace
from vyrtuous.inc.helpers import *
from vyrtuous.utils.setup_logging import logger
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.message_service import MessageService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.administrator import Administrator, AdministratorRole
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.ban import Ban
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.duration import DurationObject
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.flag import Flag
from vyrtuous.utils.history import History
from vyrtuous.utils.invincibility import Invincibility
from vyrtuous.utils.moderator import Moderator
from vyrtuous.utils.stage import Stage
from vyrtuous.utils.state import State
from vyrtuous.utils.server_mute import ServerMute
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.video_room import VideoRoom
from vyrtuous.utils.voice_mute import VoiceMute
import asyncio
import discord
import time

class EventListeners(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.emoji = Emojis()
        self.flags = []
        self.message_service = MessageService(self.bot, self.db_pool)
        self.join_log = defaultdict(list)
        self._ready_done = False
        self.deleted_rooms = {}
        self.member_service = MemberService()

    async def cog_load(self):
        VideoRoom.video_rooms = await VideoRoom.fetch_all()
        for room in VideoRoom.video_rooms:
            channel = self.bot.get_channel(room.channel_snowflake)
            if channel:
                try:
                    await channel.edit(status='Video-Only Room', reason='Enforce default video-only status')
                except discord.Forbidden as e:
                    pass
        self.flags = await Flag.fetch_all()

    
    @commands.Cog.listener()
    async def on_guild_channel_grant(self, channel: discord.abc.GuildChannel):
        guild = channel.guild
        name = channel.name
        for c in guild.channels:
            if c.id != channel.id and c.name == name:
                return
        async with self.bot.db_pool.acquire() as conn:
            room = self.deleted_rooms.pop(name, None)
            if not room:
                room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=guild.id, room_name=name)
            if room:
                old_id = room.channel_snowflake
            await Alias.update_by_source_and_target(source_channel_snowflake=old_id, target_channel_snowflake=channel.id)
            await Ban.update_by_source_and_target(source_channel_snowflake=old_id, target_channel_snowflake=channel.id)
            await Cap.update_by_source_and_target(source_channel_snowflake=old_id, target_channel_snowflake=channel.id)
            await Coordinator.update_by_source_and_target(source_channel_snowflake=old_id, target_channel_snowflake=channel.id)
            await Moderator.update_by_source_and_target(source_channel_snowflake=old_id, target_channel_snowflake=channel.id)
            await Stage.update_by_source_and_target(source_channel_snowflake=old_id, target_channel_snowflake=channel.id)
            await TemporaryRoom.update_by_source_and_target(source_channel_snowflake=old_id, target_channel_snowflake=channel.id)
            await TextMute.update_by_source_and_target(source_channel_snowflake=old_id, target_channel_snowflake=channel.id)
            await VoiceMute.update_by_source_and_target(source_channel_snowflake=old_id, target_channel_snowflake=channel.id)
            bans = await Ban.fetch_by_guild(guild_snowflake=guild.id)
            for ban in bans:
                member = self.bot.get_user(ban.member_snowflake)
                await channel.set_permissions(member, view_channel=False)
            text_mutes = await TextMute.fetch_by_guild(guild_snowflake=guild.id)
            for text_mute in text_mutes:
                member = self.bot.get_user(text_mute.member_snowflake)
                await channel.set_permissions(member, send_messages=False)
                    
    @commands.Cog.listener()                
    async def on_guild_channel_delete(self, channel: discord.abc.GuildChannel):
        room = await TemporaryRoom.fetch_by_channel_and_guild(channel_snowflake=channel.id, guild_snowflake=channel.guild.id)
        if room:
            self.deleted_rooms[channel.name] = room
    
    @commands.Cog.listener()
    async def on_guild_channel_update(self, before, after):
        if before.name == after.name:
            return
        await Ban.update_by_source_and_target(source_channel_snowflake=after.id, target_channel_snowflake=before.id)
        await Cap.update_by_source_and_target(source_channel_snowflake=after.id, target_channel_snowflake=before.id)
        await Coordinator.update_by_source_and_target(source_channel_snowflake=after.id, target_channel_snowflake=before.id)
        await Moderator.update_by_source_and_target(source_channel_snowflake=after.id, target_channel_snowflake=before.id)
        await Stage.update_by_source_and_target(source_channel_snowflake=after.id, target_channel_snowflake=before.id)
        await TextMute.update_by_source_and_target(source_channel_snowflake=after.id, target_channel_snowflake=before.id)
        await VoiceMute.update_by_source_and_target(source_channel_snowflake=after.id, target_channel_snowflake=before.id)

    # Done
    @commands.Cog.listener()
    async def on_voice_state_update(self, member: discord.Member, before: discord.VoiceState, after: discord.VoiceState):
        if member.id == self.bot.user.id:
            return
        if not after.channel:
            VideoRoom.cancel_task((member.guild.id, member.id))
            return
        for video_room in VideoRoom.video_rooms:
            if after.channel.id != video_room.channel_snowflake:
                continue
            if not after.self_video and after.channel != before.channel and after.channel.permissions_for(after.channel.guild.me).send_messages:
                await VideoRoom.enforce_video_message(channel_snowflake=after.channel.id, member_snowflake=member.id, message=f'{self.emoji.get_random_emoji()} Hi {member.mention}, {after.channel.mention} is a video only room. You have 5 minutes to turn on your camera!')
            key = (member.guild.id, member.id)
            if before.channel != after.channel:
                VideoRoom.cancel_task(key)
                if not after.self_video:
                    task = asyncio.create_task(VideoRoom.enforce_video(member, after.channel, 300))
                    VideoRoom.video_tasks[key] = task
                break
            if before.self_video and not after.self_video:
                VideoRoom.cancel_task(key)
                task = asyncio.create_task(VideoRoom.enforce_video(member, after.channel, 60))
                VideoRoom.video_tasks[key] = task
                break
            if not before.self_video and after.self_video:
                VideoRoom.cancel_task(key)
                break
        allowed = True
        if before.channel == after.channel and before.mute == after.mute and before.self_mute == after.self_mute:
            allowed = False
        if member.bot:
            allowed = False
        if not allowed:
            return
        # member_permission_role = await is_system_owner_developer_guild_owner_administrator_coordinator_moderator_via_channel_member(after.channel, member)

        target = 'user'
        # if after.channel:
            # stage = await Stage.fetch_stage_by_channel(after.channel)
            # temporary_stage_coordinator_ids = await stage.fetch_coordinator_temporary_stage_coordinator_ids(member, after.channel)
            # if stage:
            #     target = 'room'
            #     stage.send_stage_ask_to_speak_message(join_log=self.join_log, member=member)
            # else:
            #     target = 'user'
            # if stage and (member.id not in temporary_stage_coordinator_ids) and (member_permission_role in ('Moderator', 'Everyone')) and (before.channel != after.channel):
            #      expires_in = stage.expires_in
            #      await conn.execute('''
            #          INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_in, target, room_name)
            #          VALUES ($1, $2, $3, $4, 'room', $5)
            #          ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
            #          DO UPDATE SET expires_in = EXCLUDED.expires_in
            #      ''', member.guild.id, member.id, after.channel.id, expires_in, after.channel.name)
        server_mute = await ServerMute.fetch_by_member(member.id)
        if server_mute:
            if member.guild.id == server_mute.guild_snowflake:
                if not after.mute:
                    try:
                        await member.edit(mute=True, reason='Server mute is active.')
                    except discord.Forbidden as e:
                        logger.debug(f'No permission to edit mute for {member.display_name}')
                    except discord.HTTPException as e:
                        logger.debug(f'Failed to edit mute for {member.display_name}: {str(e).capitalize()}')
                return
        if after.channel:                    
            should_be_muted = False
            voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=after.channel.id, guild_snowflake=after.channel.guild.id, member_snowflake=member.id, target='user')
            if voice_mute:
                should_be_muted = True
            if not before.mute and before.channel == after.channel and after.mute:
                if member.id in Invincibility.get_invincible_members():
                    embed = discord.Embed(
                        title=f'\u1F4AB {member.display_name} is a hero!',
                        description=f'{member.display_name} cannot be muted.',
                        color=discord.Color.gold()
                    )
                    embed.set_thumbnail(url=member.display_avatar.url)
                    await after.channel.send(embed=embed)
                elif should_be_muted == False:
                    expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
                    voice_mute = VoiceMute(channel_snowflake=after.channel.id, expires_in=expires_in, guild_snowflake=after.channel.guild.id, member_snowflake=member.id, reason='No reason provided.', target=target)
                    await voice_mute.create()       
                    should_be_muted = True           
                    alias = SimpleNamespace(alias_type='voice_mute')
                    duration = DurationObject('1h')
                    await History.send_entry(alias, after.channel, duration, 'Role-specific', True, False, member, None, 'Right-click voice-mute.')
            if before.mute and not after.mute and before.channel == after.channel:
                await VoiceMute.delete_by_channel_guild_member_and_target(channel_snowflake=before.channel.id, guild_snowflake=before.channel.guild.id, member_snowflake=member.id, target=target)
                should_be_muted = False
                alias = SimpleNamespace(alias_type='unvoice_mute')
                duration = DurationObject('0')
                await History.send_entry(alias, after.channel, duration, 'Role-specific', True, False, member, None, 'Right-click undo voice-mute.')
            if after.mute != should_be_muted:
                try:
                    await member.edit(mute=should_be_muted, reason=f'Setting mute to {should_be_muted} in {after.channel.name}')
                except discord.Forbidden as e:
                    logger.debug(f'No permission to edit mute for {member.display_name}')
                except discord.HTTPException as e:
                    logger.debug(f'Failed to edit mute for {member.display_name}: {str(e).capitalize()}')
            await self.print_flags(member, after.channel)
                
#                    explicit_deny_roles = []
#                    for role in member.roles:
#                        ow = after_channel.overwrites_for(role)
#                        if ow.speak is False:
#                            explicit_deny_roles.append(role)
#                    if explicit_deny_roles:
#                        try:
#                            await member.move_to(after_channel)
#                            await logger.debug(
#                                f"🔇 Auto-muted {member.mention} in **{after_channel.name}** "
#                                f"due to explicit speak deny from roles: "
#                                f"{', '.join(r.name for r in explicit_deny_roles)}"
#                            )
#                        except Exception as e:
#                            logger.debug(
#                                f"⚠ Failed to auto-mute {member.mention} in "
#                                f"**{after_channel.name}** — `{str(e).capitalize()}`"
#                            )
#                        except discord.HTTPException as e:
#                            logger.debug(f'Failed to mute {member.display_name}: {str(e).capitalize()}')
            
    @commands.Cog.listener()
    async def on_member_join(self, member: discord.Member):
        user_id = member.id
        guild = member.guild
        async with self.db_pool.acquire() as conn:
            bans = await Ban.fetch_by_guild_and_member(guild_snowflake=guild.id, member_snowflake=user_id)
            text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=guild.id, member_snowflake=user_id)
            if bans:
                for ban in bans:
                    channel = guild.get_channel(ban.channel_snowflake)
                    try:
                        overwrite = channel.overwrites_for(member)
                        overwrite.view_channel = False
                        await channel.set_permissions(member, overwrite=overwrite, reason='Reinstating active channel ban')
                    except discord.Forbidden as e:
                        logger.warning(f'Missing permissions to ban in channel {channel.id}')
                    except discord.HTTPException as e:
                        logger.warning(f'Failed to apply ban for {member} in {channel.id}: {str(e).capitalize()}')
            if text_mutes:
                for text_mute in text_mutes:
                    channel = guild.get_channel(text_mute.channel_snowflake)
                    try:
                        overwrite = channel.overwrites_for(member)
                        overwrite.send_messages = False
                        await channel.set_permissions(member, overwrite=overwrite, reason='Reinstating text mute')
                    except discord.Forbidden as e:
                        logger.warning(f'Missing permissions to text mute in channel {channel.id}')
                    except discord.HTTPException as e:
                        logger.warning(f'Failed to apply text mute for {member} in {channel.id}: {str(e).capitalize()}')

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
        if not (message.guild and message.author.id != self.bot.user.id) or self.config['release_mode'] == False:
            return
        prefix = self.config['discord_command_prefix']
        if not message.content.startswith(prefix):
            return
        content = message.content[len(prefix):].strip()
        if not content:
            return
        parts = content.split()
        alias_name = parts[0]
        args = parts[1:]
        alias = await Alias.fetch_by_guild_and_name(alias_name=alias_name, guild_snowflake=message.guild.id)
        if not alias:
            return
        state = State(message)
        try:
            channel_obj = await self.channel_service.resolve_channel(message, alias.channel_snowflake)
            member = args[0] if len(args) > 0 else None
            member_obj = await self.member_service.resolve_member(message, member)
            await has_equal_or_higher_role(message, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, sender_snowflake=message.author.id)
            check_not_self(message, member_snowflake=member_obj.id)
            if not alias.handler:
                raise Exception(f'\U000026A0\U0000FE0F No alias handler exists for {alias.alias_name}.')
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        existing_guestroom_alias_event = await Alias.get_existing_guestroom_alias_event(alias=alias, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)
        target = args[1] if len(args) > 1 else '24h'
        is_reason_modification = target in ['+', '-', '=']
        is_duration_modification = target.startswith(('+', '-', '=')) and not is_reason_modification
        executor_role = await is_system_owner_developer_guild_owner_administrator_coordinator_moderator_via_channel_member(channel_snowflake=alias.channel_snowflake, guild_snowflake=message.guild.id, member_snowflake=message.author.id)
        if executor_role == 'Everyone':
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F You are not permitted to {alias.alias_type} users.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        if is_reason_modification and executor_role in ('Moderator', 'Everyone'):
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F You are not permitted to modify {alias.alias_type}s.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        await alias.handler(alias, args, channel_obj, executor_role, existing_guestroom_alias_event, is_duration_modification, is_reason_modification, member_obj, message, state)
        
    @commands.Cog.listener()
    async def on_command_error(self, ctx, error):
        state = State(ctx)
        try:
            match type(error):
                case commands.BadArgument:
                    return await state.end(error=f'\u274C {error}')    
                case commands.CheckFailure:
                    return await state.end(error=f'\u274C {error}')
                case commands.MissingRequiredArgument:
                    missing = error.param.name
                    return await state.end(error=f'\u274C Missing required argument: `{missing}`')
                case ValueError:
                    return await state.end(error=f'\u274C {error}')
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    @commands.Cog.listener()
    async def on_app_command_error(self, interaction, error):
        state = State(interaction)
        try:
            match type(error):
                case app_commands.BadArgument:
                    return await state.end(error=f'\u274C {error}') 
                case app_commands.CheckFailure:
                    return await state.end(error=f'\u274C {error}')
                case ValueError:
                    return await state.end(error=f'\u274C {error}') 
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
#    @commands.Cog.listener()
#    async def on_command(self, ctx):
#        await ctx.send("Bot is currently down. Changes will not be saved permanently.")

    @commands.Cog.listener()
    async def on_ready(self):
        if getattr(self, '_ready_done', False):
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
        administrator_roles = await AdministratorRole.fetch_by_guild(guild_snowflake=after.guild.id)
        for administrator_role in administrator_roles:
            administrator_role_snowflakes.append(administrator_role.role_snowflake)
        relevant_added_roles = added_roles & set(administrator_role_snowflakes)
        relevant_removed_roles = removed_roles & set(administrator_role_snowflakes)
        if not relevant_added_roles and not relevant_removed_roles:
            return
        administrator = await Administrator.fetch_by_guild_and_member(member_snowflake=after.id)
        if not administrator and relevant_added_roles:
            administrator = Administrator(guild_snowflake=after.guild.id, member_snowflake=after.id, role_snowflakes=list(after_role_snowflakes))
            await administrator.grant()
            return
        if administrator:
            for role_snowflake in relevant_added_roles:
                if role_snowflake not in administrator.role_snowflakes:
                    administrator.update_by_new_role(role_snowflake) 
            for role_snowflake in relevant_removed_roles:
                if role_snowflake in administrator.role_snowflakes:
                    administrator.update_by_removed_role(role_snowflake)
            remaining_admin_roles = set(administrator.role_snowflakes) & after_role_snowflakes
            if not remaining_admin_roles:
                await administrator.revoke()

    @commands.Cog.listener()
    async def on_guild_role_delete(self, role: discord.Role):
        async with self.bot.db_pool.acquire() as conn:
            administrators = await Administrator.fetch_by_role(role_snowflake=role.id)
            for administrator in administrators:
                role_snowflakes = set(administrator.role_snowflake)
                guild_snowflakes = set(administrator.guild_snowflake)
                role_snowflakes.discard(role.id)
                if not role_snowflakes:
                    guild_snowflakes.discard(role.guild)
                Administrator.update_guilds_and_roles_by_member(guild_snowflake=list(guild_snowflakes), member_snowflake=administrator.member_snowflake, role_snowflakes=list(role_snowflakes))

        
    async def print_flags(self, member: discord.Member, after_channel: discord.abc.GuildChannel):
        if after_channel.id == 1222056499959042108:
            for flag in self.flags:
                if flag.channel_snowflake == after_channel.id and flag.member_snowflake == member.id:
                    embed = discord.Embed(
                        title=f'\u26A0\uFE0F {member.display_name} is flagged',
                        color=discord.Color.red()
                    )
                    embed.set_thumbnail(url=member.display_avatar.url)
                    embed.add_field(name=f'Channel: {after_channel.mention}', value=f'Reason: {flag.reason}', inline=False)
                    now = time.time()
                    self.join_log[member.id] = [t for t in self.join_log[member.id] if now - t < 300]
                    if len(self.join_log[member.id]) < 1:
                        self.join_log[member.id].append(now)
                        await after_channel.send(embed=embed)
            
async def setup(bot: DiscordBot):
    await bot.add_cog(EventListeners(bot))
