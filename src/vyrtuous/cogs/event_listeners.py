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
from datetime import datetime, timezone
from discord.ext import commands
from vyrtuous.inc.helpers import *
from vyrtuous.utils.setup_logging import logger
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.discord_message_service import DiscordMessageService, Paginator
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.stage import Stage
from vyrtuous.utils.statistics import Statistics
from vyrtuous.utils.vegans import Vegans
from vyrtuous.utils.temporary_room import TemporaryRoom

import discord
import inspect
import time

class EventListeners(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.db_pool)
        self.join_log = defaultdict(list)
        self._ready_done = False
        self.deleted_rooms = {}
        
    async def fetch_active_bans(self, conn, guild_id: int, user_id: int):
        return await conn.fetch('''
            SELECT channel_id, expires_at
            FROM active_bans
            WHERE guild_id = $1
              AND discord_snowflake = $2
              AND (expires_at IS NULL OR expires_at > NOW())
        ''', guild_id, user_id)
        
    async def fetch_active_text_mutes(self, conn, guild_id: int, user_id: int):
        return await conn.fetch('''
            SELECT channel_id
            FROM active_text_mutes
            WHERE guild_id = $1
              AND discord_snowflake = $2
              AND (expires_at IS NULL OR expires_at > NOW())
        ''', guild_id, user_id)

    @commands.Cog.listener()
    async def on_guild_channel_create(self, channel: discord.abc.GuildChannel):
        guild = channel.guild
        name = channel.name
        for c in guild.channels:
            if c.id != channel.id and c.name == name:
                return
        async with self.bot.db_pool.acquire() as conn:
            room = self.deleted_rooms.pop(name, None)
            if not room:
                room = await TemporaryRoom.fetch_temporary_room_by_guild_and_room_name(guild=guild, room_name=name)
            if room:
                old_name = room.room_name
                old_id = room.channel_id
                await room.update_temporary_room_name_and_room_snowflake(channel=channel, room_name=channel.name)
                aliases = await Alias.fetch_command_aliases_by_channel_id(guild_id=channel.guild.id, channel_id=old_id)
                if aliases:
                    for alias in aliases:
                        await alias.update_command_aliases_with_channel(channel=channel)
                try:
                    old_stage = await Stage.fetch_stage_by_guild_id_and_channel_name(guild_id=guild.id, channel_name=old_name)
                except Exception:
                    old_stage = None
                if old_stage:
                    old_stage_temporary_coordinator_ids = await Stage.fetch_stage_temporary_coordinator_ids_by_guild_id_and_channel_name(guild_id=guild.id, channel_name=old_name)
                    await old_stage.update_stage_by_channel_id_name(channel_id=channel.id, channel_name=channel.name)
                await conn.execute('''
                    UPDATE users
                    SET coordinator_channel_ids = array_replace(coordinator_channel_ids, $1, $2),
                        updated_at = NOW()
                    WHERE $1 = ANY(coordinator_channel_ids)
                ''', old_id, channel.id)
                await conn.execute('''
                    UPDATE users
                    SET moderator_channel_ids = array_replace(moderator_channel_ids, $1, $2),
                        updated_at = NOW()
                    WHERE $1 = ANY(moderator_channel_ids)
                ''', old_id, channel.id)
            await conn.execute('UPDATE active_bans SET channel_id=$3 WHERE guild_id=$1 AND room_name=$2', guild.id, name, channel.id)
            await conn.execute('UPDATE active_text_mutes SET channel_id=$3 WHERE guild_id=$1 AND room_name=$2', guild.id, name, channel.id)
            await conn.execute('UPDATE active_voice_mutes SET channel_id=$3 WHERE guild_id=$1 AND room_name=$2', guild.id, name, channel.id)
            await conn.execute('UPDATE active_caps SET channel_id=$3 WHERE guild_id=$1 AND room_name=$2', guild.id, name, channel.id)
            banned_users = []
            rows = await conn.fetch('SELECT discord_snowflake FROM active_bans WHERE guild_id=$1 AND room_name=$2', guild.id, name)
            banned_users = [guild.get_member(r['discord_snowflake']) for r in rows if guild.get_member(r['discord_snowflake'])]
            for u in banned_users:
                if u:
                    await channel.set_permissions(u, view_channel=False)
            rows_mutes = await conn.fetch('SELECT discord_snowflake FROM active_text_mutes WHERE guild_id=$1 AND room_name=$2', guild.id, name)
            text_muted_users = [guild.get_member(r['discord_snowflake']) for r in rows_mutes if guild.get_member(r['discord_snowflake'])]
            for u in text_muted_users:
                if u:
                    await channel.set_permissions(u, send_messages=False)
                    
    @commands.Cog.listener()                
    async def on_guild_channel_delete(self, channel: discord.abc.GuildChannel):
        room = await TemporaryRoom.fetch_temporary_room_by_channel(channel=channel)
        if room:
            self.deleted_rooms[channel.name] = room
    
    @commands.Cog.listener()
    async def on_guild_channel_update(self, before, after):
        if before.name == after.name:
            return
        guild_id = after.guild.id
        channel_id = after.id
        new_name = after.name
        old_name = before.name
        async with self.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE active_stages
                SET room_name = $4
                WHERE guild_id = $1 AND channel_id = $2 AND room_name = $3
            ''', guild_id, channel_id, old_name, new_name)
            await conn.execute('''
                UPDATE stage_coordinators
                SET room_name = $4
                WHERE guild_id = $1 AND channel_id = $2 AND room_name = $3
            ''', guild_id, channel_id, old_name, new_name)
            await conn.execute('''
                UPDATE active_voice_mutes
                SET room_name = $4
                WHERE guild_id = $1 AND channel_id = $2 AND room_name = $3
            ''', guild_id, channel_id, old_name, new_name)
            await conn.execute('''
                UPDATE active_bans
                SET room_name = $4
                WHERE guild_id = $1 AND channel_id = $2 AND room_name = $3
            ''', guild_id, channel_id, old_name, new_name)
            await conn.execute('''
                UPDATE active_text_mutes
                SET room_name = $4
                WHERE guild_id = $1 AND channel_id = $2 AND room_name = $3
            ''', guild_id, channel_id, old_name, new_name)
            await conn.execute('''
                UPDATE active_caps
                SET room_name = $4
                WHERE guild_id = $1 AND channel_id = $2 AND room_name = $3
            ''', guild_id, channel_id, old_name, new_name)

    # Done
    @commands.Cog.listener()
    async def on_voice_state_update(self, member: discord.Member, before: discord.VoiceState, after: discord.VoiceState) -> None:
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

            user_data = await conn.fetchrow('SELECT server_mute_guild_ids FROM users WHERE discord_snowflake = $1', member.id)
            server_mute_guild_ids = user_data['server_mute_guild_ids'] or [] if user_data else []
            if member.guild.id in server_mute_guild_ids:
                return    

            if after.channel:                    
                should_be_muted = False
                if not before.mute and after.mute:
                    if member.id in Vegans.get_vegans():
                        embed = discord.Embed(
                            title=f'\u1F4AB {member.display_name} is a hero!',
                            description=f'{member.display_name} cannot be muted.',
                            color=discord.Color.gold()
                        )
                        embed.set_thumbnail(url=member.display_avatar.url)
                        await after.channel.send(embed=embed)
                    else:
                        await conn.execute('''
                            INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, room_name, target)
                            VALUES ($1, $2, $3, $4, $5)
                            ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
                            DO UPDATE SET expires_at = NOW() + interval '1 hour'
                        ''', member.guild.id, member.id, after.channel.id, after.channel.name, target)           
                        should_be_muted = True                
                if before.mute and not after.mute and before.channel:
                    result = await conn.execute('''
                        DELETE FROM active_voice_mutes
                        WHERE guild_id = $1
                            AND discord_snowflake = $2
                            AND channel_id = $3
                            AND room_name = $4
                            AND target = $5
                            AND expires_at IS NOT NULL
                    ''', member.guild.id, member.id, before.channel.id, after.channel.name, target)
                existing_mute_row = await conn.fetchrow('''
                    SELECT expires_at
                    FROM active_voice_mutes
                    WHERE guild_id = $1
                        AND discord_snowflake = $2
                        AND channel_id = $3
                        AND room_name = $4
                        AND target = $5
                ''', member.guild.id, member.id, after.channel.id, after.channel.name, target)
                if existing_mute_row:
                    should_be_muted = True
                if after.mute != should_be_muted:
                    try:
                        await member.edit(mute=should_be_muted, reason=f'Setting mute to {should_be_muted} in {after.channel.name}')
                    except discord.Forbidden:
                        logger.debug(f'No permission to edit mute for {member.display_name}')
                    except discord.HTTPException as e:
                        logger.debug(f'Failed to edit mute for {member.display_name}: {e}')
                
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
#                            logger.debug(f'Failed to mute {member.display_name}: {e}')
            
    @commands.Cog.listener()
    async def on_member_join(self, member: discord.Member) -> None:
        user_id = member.id
        guild = member.guild
        async with self.db_pool.acquire() as conn:
            bans = await self.fetch_active_bans(conn, guild.id, user_id)
            text_mutes = await self.fetch_active_text_mutes(conn, guild.id, user_id)
            for row in bans:
                channel = guild.get_channel(row['channel_id'])
                if not channel or not isinstance(channel, (discord.TextChannel, discord.VoiceChannel)):
                    continue
                if row['expires_at'] and row['expires_at'] < datetime.now(timezone.utc):
                    continue
                try:
                    overwrite = channel.overwrites_for(member)
                    overwrite.view_channel = False
                    await channel.set_permissions(member, overwrite=overwrite, reason='Reinstating active channel ban')
                except discord.Forbidden:
                    print(f'Missing permissions to ban in channel {channel.id}')
                except discord.HTTPException as e:
                    print(f'Failed to apply ban for {member} in {channel.id}: {e}')
            for row in text_mutes:
                channel = guild.get_channel(row['channel_id'])
                if not channel or not isinstance(channel, discord.TextChannel):
                    continue
                try:
                    overwrite = channel.overwrites_for(member)
                    overwrite.send_messages = False
                    await channel.set_permissions(member, overwrite=overwrite, reason='Reinstating text mute')
                except discord.Forbidden:
                    print(f'Missing permissions to text mute in channel {channel.id}')
                except discord.HTTPException as e:
                    print(f'Failed to apply text mute for {member} in {channel.id}: {e}')

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
        if not message.guild:
            return
        if message.author.id == self.bot.user.id:
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
        alias = await Alias.fetch_command_alias_by_guild_and_alias_name(message.guild, alias_name)
        if not alias:
            return
        await self.dispatch_alias(message, alias, args)
    
    async def dispatch_alias(self, message: discord.Message, alias: Alias, args):
        if not alias.handler:
            return
        await alias.handler(message, alias, args)
        
    @commands.Cog.listener()
    async def on_command_error(self, ctx, error):
        if isinstance(error, commands.MissingRequiredArgument):
            missing = error.param.name
            await ctx.reply(f'\U0001F6AB Missing required argument: `{missing}`')
            return
        if isinstance(error, commands.CheckFailure):
            return await send_check_failure_embed(ctx, error)
            
#    @commands.Cog.listener()
#    async def on_command(self, ctx):
#        await ctx.send("Bot is currently down. Changes will not be saved permanently.")

    @commands.Cog.listener()
    async def on_ready(self):
        await Statistics.load_channels()
        if getattr(self, "_ready_done", False):
            return
        self._ready_done = True

    @commands.Cog.listener()
    async def on_member_update(self, before: discord.Member, after: discord.Member):
        if before.roles == after.roles:
            return
        guild_id = after.guild.id
        added_roles = {r.id for r in after.roles} - {r.id for r in before.roles}
        removed_roles = {r.id for r in before.roles} - {r.id for r in after.roles}
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT administrator_role_ids, administrator_guild_ids FROM users WHERE discord_snowflake=$1',
                after.id
            )
            if not row:
                return
            admin_role_ids = set(row['administrator_role_ids'] or [])
            admin_guild_ids = set(row['administrator_guild_ids'] or [])
            if added_roles & admin_role_ids:
                admin_guild_ids.add(guild_id)
            if removed_roles & admin_role_ids:
                remaining_admin_roles = admin_role_ids & {r.id for r in after.roles}
                if not remaining_admin_roles:
                    admin_guild_ids.discard(guild_id)
            await conn.execute(
                'UPDATE users SET administrator_guild_ids=$2 WHERE discord_snowflake=$1',
                after.id,
                list(admin_guild_ids)
            )

    @commands.Cog.listener()
    async def on_guild_role_delete(self, role: discord.Role):
        guild = role.guild
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT discord_snowflake, administrator_role_ids, administrator_guild_ids FROM users WHERE $1 = ANY(administrator_role_ids)',
                role.id
            )
            for row in rows:
                role_ids = set(row['administrator_role_ids'] or [])
                guild_ids = set(row['administrator_guild_ids'] or [])
                role_ids.discard(role.id)
                if not role_ids:
                    guild_ids.discard(guild.id)
                await conn.execute(
                    'UPDATE users SET administrator_role_ids=$2, administrator_guild_ids=$3 WHERE discord_snowflake=$1',
                    row['discord_snowflake'],
                    list(role_ids),
                    list(guild_ids)
                )
        
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
