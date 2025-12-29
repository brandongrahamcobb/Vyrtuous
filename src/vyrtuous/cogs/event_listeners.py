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
from typing import Union
from vyrtuous.inc.helpers import *
from vyrtuous.utils.setup_logging import logger
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.message_service import MessageService
from vyrtuous.utils.paginator import Paginator
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.ban import Ban
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.duration import Duration
from vyrtuous.utils.moderator import Moderator
from vyrtuous.utils.stage import Stage
from vyrtuous.utils.state import State
from vyrtuous.utils.server_mute import ServerMute
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.time_to_complete import TimeToComplete
from vyrtuous.utils.invincibility import Invincibility
from vyrtuous.utils.voice_mute import VoiceMute

import discord
import time

class EventListeners(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.handler = MessageService(self.bot, self.db_pool)
        self.join_log = defaultdict(list)
        self._ready_done = False
        self.deleted_rooms = {}
        self.member_service = MemberService()

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
        allowed = True
        if before.channel == after.channel and before.mute == after.mute and before.self_mute == after.self_mute:
            allowed = False
        if member.bot:
            allowed = False
        if not allowed:
            return
        # member_permission_role = await is_owner_developer_administrator_coordinator_moderator_via_channel_member(after.channel, member)

        async with self.db_pool.acquire() as conn:
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
                #      expires_at = stage.expires_at
                #      await conn.execute('''
                #          INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, target, room_name)
                #          VALUES ($1, $2, $3, $4, 'room', $5)
                #          ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
                #          DO UPDATE SET expires_at = EXCLUDED.expires_at
                #      ''', member.guild.id, member.id, after.channel.id, expires_at, after.channel.name)
            server_mute = await ServerMute.fetch_by_member(member.id)
            if server_mute:
                if member.guild.id == server_mute.guild_snowflake:
                    return
            if after.channel:                    
                should_be_muted = False
                if not before.mute and after.mute:
                    if member.id in Invincibility.get_invincible_members():
                        embed = discord.Embed(
                            title=f'\u1F4AB {member.display_name} is a hero!',
                            description=f'{member.display_name} cannot be muted.',
                            color=discord.Color.gold()
                        )
                        embed.set_thumbnail(url=member.display_avatar.url)
                        await after.channel.send(embed=embed)
                    else:
                        expires_at = datetime.utcnow() + timedelta(hours=1)
                        voice_mute = await VoiceMute(channel_snowflake=after.channel.id, expires_at=expires_at, guild_snowflake=after.channel.guild.id, member_snowflake=member.id, target=target)
                        await voice_mute.create()       
                        should_be_muted = True 
                if not should_be_muted:               
                    if before.mute and not after.mute and before.channel:
                        await VoiceMute.delete_by_channel_guild_member_and_target(channel_snowflake=before.channel.id, guild_snowflake=before.channel.guild.id, member_snowflake=member.id, target=target)
                    voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=after.channel.id, guild_snowflake=after.channel.guild.id, member_snowflake=member.id, target="user")
                if voice_mute:
                    should_be_muted = True
                if after.mute != should_be_muted:
                    try:
                        await member.edit(mute=should_be_muted, reason=f'Setting mute to {should_be_muted} in {after.channel.name}')
                    except discord.Forbidden:
                        logger.debug(f'No permission to edit mute for {member.display_name}')
                    except discord.HTTPException as e:
                        logger.debug(f'Failed to edit mute for {member.display_name}: {e}.')
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
#                                f"**{after_channel.name}** — `{e}`"
#                            )
#                        except discord.HTTPException as e:
#                            logger.debug(f'Failed to mute {member.display_name}: {e}.')
            
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
                    except discord.Forbidden:
                        logger.warning(f'Missing permissions to ban in channel {channel.id}')
                    except discord.HTTPException as e:
                        logger.warning(f'Failed to apply ban for {member} in {channel.id}: {e}.')
            if text_mutes:
                for text_mute in text_mutes:
                    channel = guild.get_channel(text_mute.channel_snowflake)
                    try:
                        overwrite = channel.overwrites_for(member)
                        overwrite.send_messages = False
                        await channel.set_permissions(member, overwrite=overwrite, reason='Reinstating text mute')
                    except discord.Forbidden:
                        logger.warning(f'Missing permissions to text mute in channel {channel.id}')
                    except discord.HTTPException as e:
                        logger.warning(f'Failed to apply text mute for {member} in {channel.id}: {e}.')

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if after.author.bot:
            return
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)
    
    @commands.Cog.listener()
    async def on_message(self, message: discord.Message):
        # if not message.guild or message.author.id == self.bot.user.id:
        #     return
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
            if member_obj.id == message.guild.me.id:
                raise Exception(f"\U000026A0\U0000FE0F You cannot {alias.alias_type} {message.guild.me.mention}.")
            if not alias.handler:
                raise Exception(f"\U000026A0\U0000FE0F No alias handler exists for {alias.alias_name}.")
        except Exception as e:
            try:
                return await state.end(warning=f"\U000026A0\U0000FE0F {e}")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        existing_guestroom_alias_event = await Alias.get_existing_guestroom_alias_event(alias=alias, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)
        target = args[1] if len(args) > 1 else '24h'
        is_reason_modification = target in ['+', '-', '=']
        is_duration_modification = target.startswith(('+', '-', '=')) and not is_reason_modification
        executor_role = await is_owner_developer_administrator_coordinator_moderator(message)
        if executor_role == 'Everyone' or (is_reason_modification and executor_role in ('Moderator', 'Everyone')):
            try:
                return await state.end(warning=f"\U000026A0\U0000FE0F You are not permitted to modify {alias.alias_type}s or {alias.alias_type} users.")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        await alias.handler(alias, args, channel_obj, executor_role, existing_guestroom_alias_event, is_duration_modification, is_reason_modification, member_obj, message, state)
        
    @commands.Cog.listener()
    async def on_command_error(self, ctx, error):
        if isinstance(error, (commands.BadArgument, commands.CheckFailure)):
            return await ctx.reply(f'\U000026A0\U0000FE0F {error}')
        if isinstance(error, commands.MissingRequiredArgument):
            missing = error.param.name
            return await ctx.reply(f'\U000026A0\U0000FE0F Missing required argument: `{missing}`')
    
    @commands.Cog.listener()
    async def on_app_command_error(self, interaction, error):
        if isinstance(error, (app_commands.BadArgument, app_commands.CheckFailure)):
            return await self.handler.send_message(interaction, f'\U000026A0\U0000FE0F {error}')
            
#    @commands.Cog.listener()
#    async def on_command(self, ctx):
#        await ctx.send("Bot is currently down. Changes will not be saved permanently.")

    @commands.Cog.listener()
    async def on_ready(self):
        if getattr(self, "_ready_done", False):
            return
        self._ready_done = True

    @commands.Cog.listener()
    async def on_member_update(self, before: discord.Member, after: discord.Member):
        if before.roles == after.roles:
            return
        after_role_ids = {r.id for r in after.roles}
        before_role_ids = {r.id for r in before.roles}
        added_roles = after_role_ids - before_role_ids
        removed_roles = before_role_ids - after_role_ids
        administrator = await Administrator.fetch_member(member_snowflake=after.id)
        if not administrator:
            return
        guild_snowflakes = set(administrator.guild_snowflakes)
        role_snowflakes = set(administrator.role_snowflakes)
        if added_roles and (added_roles & role_snowflakes):
            guild_snowflakes.add(after.guild.id)
        if removed_roles and (removed_roles & role_snowflakes):
            remaining_admin_roles = role_snowflakes & after_role_ids
            if not remaining_admin_roles:
                guild_snowflakes.discard(after.guild.id)
        await Administrator.update_guilds_and_roles_for_member(guild_snowflake=list(guild_snowflakes), member_snowflake=after.id, role_snowflake=list(role_snowflakes))

    @commands.Cog.listener()
    async def on_guild_role_delete(self, role: discord.Role):
        async with self.bot.db_pool.acquire() as conn:
            administrators = await Administrator.fetch_members_by_role(role_snowflake=role.id)
            for administrator in administrators:
                role_snowflakes = set(administrator.role_snowflake)
                guild_snowflakes = set(administrator.guild_snowflake)
                role_snowflakes.discard(role.id)
                if not role_snowflakes:
                    guild_snowflakes.discard(role.guild)
                Administrator.update_guilds_and_roles_by_member(guild_snowflake=list(guild_snowflakes), member_snowflake=administrator.member_snowflake, role_snowflakes=list(role_snowflakes))

        
    async def print_flags(self, member: discord.Member, after_channel: discord.abc.GuildChannel):
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_id, discord_snowflake, reason
                FROM active_flags
                WHERE guild_id = $1 AND discord_snowflake = $2
            ''', member.guild.id, member.id)
            if rows:
                grouped = {}
                for row in rows:
                    grouped.setdefault(row['channel_id'], []).append(row)
                context_records = grouped.get(after_channel.id)
                if after_channel.id == 1222056499959042108 and context_records:
                    if context_records and after_channel.id:
                        embeds = []
                        embed = discord.Embed(
                            title=f'\u26A0\uFE0F {member.display_name} is flagged',
                            color=discord.Color.red()
                        )
                        embed.set_thumbnail(url=member.display_avatar.url)
                        for record in context_records:
                            reason = record['reason'] or 'No reason provided'
                            embed.add_field(name=f'Channel: {after_channel.mention}', value=f'Reason: {reason}', inline=False)
                        other_channels = [ch_id for ch_id in grouped.keys() if ch_id != after_channel.id]
                        if other_channels:
                            ch_mentions = []
                            for ch_id in other_channels:
                                ch = member.guild.get_channel(ch_id)
                                if not ch:
                                    ch = await member.guild.fetch_channel(ch_id)
                                ch_mentions.append(ch.mention if ch else f'Channel ID `{ch_id}`')
                            embed.add_field(name='Other flagged channels', value='\n'.join(ch_mentions), inline=False)
                        embeds.append(embed)
                        for ch_id in other_channels:
                            records = grouped[ch_id]
                            ch = member.guild.get_channel(ch_id)
                            ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                            embed = discord.Embed(
                                title=f'\u26A0\uFE0F {member.display_name} is flagged in {ch_name}',
                                color=discord.Color.red()
                            )
                            embed.set_thumbnail(url=member.display_avatar.url)
                            for record in records:
                                reason = record['reason'] or 'No reason provided'
                                embed.add_field(name='Channel', value=f'{ch_name}\nReason: {reason}', inline=False)
                            embeds.append(embed)
                        now = time.time()
                        self.join_log[member.id] = [t for t in self.join_log[member.id] if now - t < 300]
                        if len(self.join_log[member.id]) < 1:
                            self.join_log[member.id].append(now)
                            if len(embeds) == 1:
                                await after_channel.send(embed=embeds[0])
                            else:
                                paginator = Paginator(self.bot, after_channel, embeds)
                                await paginator.start()
            
async def setup(bot: DiscordBot):
    await bot.add_cog(EventListeners(bot))
