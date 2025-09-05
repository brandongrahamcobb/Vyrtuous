''' commands.py

    Copyright (C) 2024  github.com/brandongrahamcobb

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
import os
import random
import subprocess
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from typing import Any, Optional, Union

from discord.ext.commands import Command

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.check_service import *
from vyrtuous.service.discord_message_service import DiscordMessageService, Paginator
from vyrtuous.utils.setup_logging import logger

VEGAN_EMOJIS = [
    "\U0001F436",
    "\U0001F431",
    "\U0001F42D",
    "\U0001F439",
    "\U0001F430",
    "\U0001F98A",
    "\U0001F43B",
    "\U0001F43C",
    "\U0001F428",
    "\U0001F42F",
    "\U0001F981",
    "\U0001F42E",
    "\U0001F437",
    "\U0001F43D",
    "\U0001F438",
    "\U0001F435",
    "\U0001F412",
    "\U0001F98D",
    "\U0001F9A7",
    "\U0001F414",
    "\U0001F427",
    "\U0001F426",
    "\U0001F424",
    "\U0001F423",
    "\U0001F425",
    "\U0001F986",
    "\U0001F9A2",
    "\U0001F989",
    "\U0001F99A",
    "\U0001F99C",
    "\U0001F43A",
    "\U0001F99D",
    "\U0001F9A8",
    "\U0001F9A1",
    "\U0001F417",
    "\U0001F434",
    "\U0001F984",
    "\U0001F41D",
    "\U0001F41B",
    "\U0001F98B",
    "\U0001F40C",
    "\U0001F41E",
    "\U0001F40C",
    "\U0001FAB2",
    "\U0001F9F3",
    "\U0001F997",
    "\U0001F577",
    "\U0001F982",
    "\U0001F422",
    "\U0001F40D",
    "\U0001F98E",
    "\U0001F996",
    "\U0001F995",
    "\U0001F419",
    "\U0001F991",
    "\U0001F990",
    "\U0001F99E",
    "\U0001F980",
    "\U0001F421",
    "\U0001F420",
    "\U0001F41F",
    "\U0001F42C",
    "\U0001F988",
    "\U0001F433",
    "\U0001F40B",
    "\U0001F9AD",
    "\U0001F40A",
    "\U0001F406",
    "\U0001F405",
    "\U0001F403",
    "\U0001F402",
    "\U0001F42B",
    "\U0001F42A",
    "\U0001F999",
    "\U0001F992",
    "\U0001F98F",
    "\U0001F99B",
    "\U0001F418",
    "\U0001F998",
    "\U0001F9A5",
    "\U0001F9A6",
    "\U0001F9A8",
    "\U0001F9A9",
    "\U0001F54A"
]
PERMISSION_ORDER = ['Owner', 'Developer', 'Coordinator', 'Moderator', 'Everyone']

class Hybrid(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.bot.loop.create_task(self.load_server_muters())
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
        self.server_muters: dict[int, set[int]] = defaultdict(set)
#        self.backdoor = False
        
    def get_random_emoji(self):
        return random.choice(VEGAN_EMOJIS)
        
    
    @commands.command(name='admin', help='Grants server mute privileges to a member for the entire guild.')
    @is_owner()
    async def grant_server_muter(self, ctx, member: str):
        member, _ = await self.get_channel_and_member(ctx, member)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                               INSERT INTO users (user_id, server_muter_guild_ids)
                               VALUES ($1, ARRAY[$2]::BIGINT[]) ON CONFLICT (user_id) DO
                               UPDATE
                                   SET server_muter_guild_ids = (
                                   SELECT ARRAY(
                                   SELECT DISTINCT unnest(COALESCE (u.server_muter_guild_ids, '{}') || ARRAY[$2])
                                   )
                                   FROM users u WHERE u.user_id = EXCLUDED.user_id
                                   ),
                                   updated_at = NOW()
                               ''', member.id, ctx.guild.id)
        self.server_muters.setdefault(ctx.guild.id, set()).add(member.id)
        return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been granted server mute permissions.', allowed_mentions=discord.AllowedMentions.none())
        
#    @commands.command(name="toggle", hidden=True)
#    @is_owner()
#    async def toggle_feature(self, ctx: commands.Context):
#        self.backdoor = not self.backdoor
#        state = f"ON {self.get_random_emoji()}" if self.backdoor else f"OFF 🔥"
#        await ctx.send(f"{self.get_random_emoji()} Feature switched {state}.")
    
#    async def backdoor_start(self, ctx: commands.Context):
#        channel = ctx.channel
#        user_id = ctx.author.id
#        channel_id = channel.id
#        guild_id = ctx.guild.id
#        async with self.bot.db_pool.acquire() as conn:
#            row = await conn.fetchrow("""
#                SELECT coordinator_ids, coordinator_channel_ids
#                FROM users
#                WHERE user_id = $1
#            """, user_id)
#        if not row:
#            await ctx.send("❌ You are not registered in the database.")
#            return
#        coordinator_ids = row["coordinator_ids"] or []
#        coordinator_channel_ids = row["coordinator_channel_ids"] or []
#        is_coordinator = (user_id in coordinator_ids) and (channel_id in coordinator_channel_ids)
#        is_vegan_channel = "vegan" in channel.name.lower()
#        if is_coordinator and is_vegan_channel:
#            return True
#        else:
#            return False

    @commands.command(name="backup", description="Creates a backup of the database and uploads it")
    @is_owner_developer()
    async def backup(self, ctx: commands.Context):
        try:
            backup_dir = self.setup_backup_directory('/app/backups')
            backup_file = self.perform_backup(
                db_user=os.getenv("POSTGRES_USER"),
                db_name=os.getenv("POSTGRES_DATABASE"),
                db_host=os.getenv("POSTGRES_HOST"),
                db_password=os.getenv("POSTGRES_PASSWORD"),
                backup_dir=backup_dir
            )
            logger.info(f'Backup completed successfully: {backup_file}')
            await ctx.send(file=discord.File(backup_file))
        except Exception as e:
            logger.error(f'Error during database backup: {e}')
            await ctx.send(f"❌ Failed to create backup: {e}")


    @commands.command(name='alias', help='Set an alias for a cow, uncow, mute, unmute, ban, unban or flag action.')
    @is_owner_developer_coordinator()
    async def create_alias(
            self,
            ctx,
            alias_type: str = commands.parameter(description='One of: `cow`, `uncow`, `mute`, `unmute`, `ban`, `unban`, `flag`, `unflag`, `tmute`, `untmute`'),
            alias_name: str = commands.parameter(description='Alias/Pseudonym'),
            target: str = commands.parameter(description='Voice channel')
    ) -> None:
        cmd = None
        alias_type = alias_type.lower()
        valid_types = {'cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute'}
        if alias_type not in valid_types:
            return await self.handler.send_message(ctx, content=f'\U0001F525 Invalid alias type. Must be one of: {", ".join(valid_types)}')
        if not alias_name.strip():
            return await self.handler.send_message(ctx, content='\U0001F525 Alias name cannot be empty.')
        _, channel = await self.get_channel_and_member(ctx, target)
        is_owner_or_dev, _ = await check_owner_dev_coord_mod(ctx, channel)
        if not is_owner_or_dev:
            ctx._target_channel_id = channel.id
            try:
                await is_coordinator_in_channel.predicate(ctx)
            except commands.CheckFailure:
                return await self.handler.send_message(
                    ctx,
                    content=f'\U0001F525 You do not have permission for {channel.mention}.'
                )
        async with self.bot.db_pool.acquire() as conn:
            existing_alias = await conn.fetchrow(
                '''
                SELECT guild_id, channel_id
                FROM command_aliases
                WHERE alias_type = $1
                  AND alias_name = $2
                ''',
                alias_type, alias_name
            )
            if existing_alias:
                existing_channel = ctx.guild.get_channel(existing_alias['channel_id'])
                channel_mention = existing_channel.mention if existing_channel else f"<#{existing_alias['channel_id']}>"
                return await ctx.send(f'\U0001F525 Alias `{alias_name}` ({alias_type}) already exists and is set to {channel_mention}.', allowed_mentions=discord.AllowedMentions.none())
            if self.bot.get_command(alias_name):
                return await self.handler.send_message(ctx, content=f'\U0001F525 A command named `{alias_name}` already exists.')
            await conn.execute(
                '''
                INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id)
                VALUES ($1, $2, $3, $4)
                ''',
                ctx.guild.id, alias_type, alias_name, channel.id
            )
            await conn.execute('''
                               INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id, guild_id,
                                                            channel_id, reason)
                               VALUES ($1, $2, $3, $4, $5, $6)
                               ''', 'create_alias', None, ctx.author.id, ctx.guild.id, channel.id,
                               f'Created an alias: {alias_name}')
        self.bot.command_aliases.setdefault(ctx.guild.id, {}).setdefault(alias_type, {})[alias_name] = channel.id
        if alias_type == 'ban':
            cmd = self.create_ban_alias(alias_name)
        elif alias_type == 'flag':
            cmd = self.create_flag_alias(alias_name)
        elif alias_type == 'mute':
            cmd = self.create_voice_mute_alias(alias_name)
        elif alias_type == 'tmute':
            cmd = self.create_text_mute_alias(alias_name)
        elif alias_type == 'unban':
            cmd = self.create_unban_alias(alias_name)
        elif alias_type == 'cow':
            cmd = self.create_cow_alias(alias_name)
        elif alias_type == 'uncow':
            cmd = self.create_uncow_alias(alias_name)
        elif alias_type == 'unflag':
            cmd = self.create_unflag_alias(alias_name)
        elif alias_type == 'unmute':
            cmd = self.create_unmute_alias(alias_name)
        elif alias_type == 'untmute':
            cmd = self.create_untextmute_alias(alias_name)
        self.bot.add_command(cmd)
        return await ctx.send(f'{self.get_random_emoji()} Alias `{alias_name}` ({alias_type}) set to {channel.mention}.')

    def create_ban_alias(self, command_name: str) -> Command:
        @commands.command(
            name=command_name,
            help='Ban a user from a voice channel.'
        )
        @is_owner_developer_coordinator_moderator("ban")
        async def ban_alias(
            ctx,
            member: str = commands.parameter(description='Mention or user ID of the member to ban.'),
            duration_hours: Optional[str] = commands.parameter(default='24', description='Duration of ban in hours. Example 0 (permanent), 30m, 2h, 5d.'),
            *,
            reason: str = commands.parameter(default='', description='Reason for ban (required for permanent).')
        ) -> None:
            try:
                await is_owner_block(ctx, member)
            except commands.CheckFailure as e:
                logger.warning(e)
                return await self.handler.send_message(ctx, content='\U0001F525 You are not allowed to ban the owner.')
#            if self.backdoor:
#                if self.backdoor_start():
#                    return await ctx.send("🔥 You aren't vegan. Go vegan.")
            cmd = ctx.invoked_with
            member, _ = await self.get_channel_and_member(ctx, member)
            expires_at, duration_display = self.parse_duration(duration_hours)
            if (expires_at == '0' or expires_at is None) and (not is_owner_developer_coordinator("ban") or not reason.strip()):
                return await self.handler.send_message(ctx, content='\U0001F525 Reason required and coordinator-only for permanent bans.')
            static_channel_id = self.bot.command_aliases.get(ctx.guild.id, {}).get('ban', {}).get(cmd)
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F525 No channel alias mapping found for `{cmd}`.')
            async with self.bot.db_pool.acquire() as conn:
                existing_ban = await conn.fetchrow(
                    '''
                    SELECT expires_at
                    FROM active_bans
                    WHERE user_id = $1 AND channel_id = $2
                    ''',
                    member.id, static_channel_id
                )
            if existing_ban and not is_owner_developer_coordinator("ban"):
                if existing_ban['expires_at'] is None:
                    return await self.handler.send_message(
                        ctx,
                        content=f'\U0001F525 {member.mention} is already permanently banned from <#{static_channel_id}>.'
                    )
                else:
                    remaining = existing_ban['expires_at'] - discord.utils.utcnow()
                    if remaining.total_seconds() > 0:
                        hours_left = round(remaining.total_seconds() / 3600, 1)
                        return await self.handler.send_message(
                            ctx,
                            content=f'\U0001F525 {member.mention} is already banned from <#{static_channel_id}> for another {hours_left}h.'
                        )
            channel = ctx.guild.get_channel(static_channel_id)
            if not channel or not isinstance(channel, discord.VoiceChannel):
                return await self.handler.send_message(ctx, content=f'\U0001F525 Could not resolve a valid voice channel for ID `{static_channel_id}`.')
            try:
                await channel.set_permissions(
                    member,
                    view_channel=False,
                    reason=f"{self.get_random_emoji()} Banned from <#{channel.id}>: {reason or 'No reason provided'}"
                )
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content='\U0001F525 Missing permissions to deny channel access.')
            if member.voice and member.voice.channel and member.voice.channel.id == channel.id:
                try:
                    await member.move_to(None, reason="Banned from this channel")
                except discord.Forbidden:
                    await ctx.send(f"🔥️ Could not disconnect <@{member.id}> from <#{channel.id}>.", allowed_mentions=discord.AllowedMentions.none())
                except Exception as e:
                    logger.exception(f"Unexpected error while disconnecting user: {e}")
                    raise
            try:
                async with self.bot.db_pool.acquire() as conn:
                    await conn.execute(
                        '''
                        INSERT INTO active_bans (user_id, channel_id, expires_at)
                        VALUES ($1, $2, $3) ON CONFLICT (user_id, channel_id)
                        DO
                        UPDATE SET expires_at = EXCLUDED.expires_at
                        ''',
                        member.id,
                        channel.id,
                        expires_at
                    )
                    await conn.execute(
                        '''
                        UPDATE users
                        SET ban_channel_ids =
                                CASE
                                    WHEN NOT $2 = ANY (ban_channel_ids) THEN array_append(ban_channel_ids, $2)
                                    ELSE ban_channel_ids
                                    END,
                            updated_at = NOW()
                        WHERE user_id = $1
                        ''',
                        member.id, channel.id
                    )
                    await conn.execute(
                        '''
                        INSERT INTO ban_reasons (guild_id, user_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4) ON CONFLICT (guild_id, user_id, channel_id)
                        DO
                        UPDATE SET reason = EXCLUDED.reason
                        ''',
                        ctx.guild.id, member.id, channel.id, reason or 'No reason provided'
                    )
                    await conn.execute(
                        '''
                        INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                        ''',
                        'ban', member.id, ctx.author.id, ctx.guild.id, channel.id, reason or 'No reason provided'
                    )
            except Exception as e:
                logger.warning(f"Database error occurred: {e}")
                raise
            return await ctx.send(
                f'{self.get_random_emoji()} {member.mention} has been banned from <#{channel.id}> {duration_display} because: {reason or "No reason provided"}',
                allowed_mentions=discord.AllowedMentions.none()
            )
        return ban_alias
        
    @commands.command(name='coord', help='Grants coordinator access for a specific voice channel.')
    @is_owner_developer()
    async def create_coordinator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel: Optional[str] = commands.parameter(default=None, description='Mention a channel or provide its ID.'),
    ) -> None:
        try:
            await is_owner_block(ctx, member)
        except commands.CheckFailure as e:
            logger.warning(e)
            return await self.handler.send_message(ctx, content='\U0001F525 You are not allowed to make the owner a coordinator.')
        _, channel = await self.get_channel_and_member(ctx, channel)
        member, _ = await self.get_channel_and_member(ctx, member)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (user_id, coordinator_ids, coordinator_channel_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[], ARRAY[$3]::BIGINT[])
                ON CONFLICT (user_id) DO UPDATE
                SET 
                    coordinator_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(
                                COALESCE(users.coordinator_ids, ARRAY[]::BIGINT[]) || 
                                ARRAY[$2]::BIGINT[]
                            )
                        )
                    ),
                    coordinator_channel_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(
                                COALESCE(users.coordinator_channel_ids, ARRAY[]::BIGINT[]) || 
                                ARRAY[$3]::BIGINT[]
                            )
                        )
                    ),
                    updated_at = NOW()
            ''', member.id, ctx.guild.id, channel.id)
            await conn.execute('''
                               INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id, guild_id,
                                                            channel_id, reason)
                               VALUES ($1, $2, $3, $4, $5, $6)
                               ''', 'create_coordinator', member.id, ctx.author.id, ctx.guild.id, channel.id,
                               'Created a coordinator')
        return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been granted coordinator rights in {channel.mention}.', allowed_mentions=discord.AllowedMentions.none())
    
    def create_cow_alias(self, command_name: str) -> Command:
        @commands.command(
            name=command_name,
            help='Label a user as going vegan for tracking purposes.'
        )
        @is_owner_developer_coordinator_moderator("cow")
        async def going_vegan_alias(
                ctx,
                user: str = commands.parameter(description='Tag a user or include their snowflake ID.')
        ) -> None:
            cow_aliases = self.bot.command_aliases.get(ctx.guild.id, {}).get('cow', {})
            channel_id = cow_aliases.get("cow")
            member, _ = await self.get_channel_and_member(ctx, user)
            select_sql = '''
                         SELECT 1
                         FROM users
                         WHERE user_id = $1
                           AND $2 = ANY (going_vegan_channel_ids)
                         '''
            insert_sql = '''
                         INSERT INTO users (user_id, going_vegan_channel_ids)
                         VALUES ($1, ARRAY[$2]::BIGINT[]) ON CONFLICT (user_id)
                         DO
                         UPDATE SET going_vegan_channel_ids = (
                             SELECT ARRAY(
                             SELECT DISTINCT unnest(users.going_vegan_channel_ids || EXCLUDED.going_vegan_channel_ids)
                             )
                             )
                         '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    already_cowed = await conn.fetchval(select_sql, member.id, channel_id)
                    if already_cowed:
                        return await ctx.send(f'\U0001F525 <@{member.id}> is already going vegan.', allowed_mentions=discord.AllowedMentions.none())
                    await conn.execute(insert_sql, member.id, channel_id)
                    await conn.execute('''
                                       INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id,
                                                                    guild_id,
                                                                    channel_id, reason)
                                       VALUES ($1, $2, $3, $4, $5, $6)
                                       ''', 'cow', member.id, ctx.author.id, ctx.guild.id, channel_id,
                                       'Cowed a user')
                    await ctx.send(f'{self.get_random_emoji()} <@{member.id}> is going vegan!!! WOOO :)', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                return await self.handler.send_message(ctx, content=f'Database error: {e}')
                raise
        return going_vegan_alias
        
    @commands.command(name='dev', help='Elevates a user\'s permissions to a bot developer.')
    @is_owner()
    async def create_developer(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ) -> None:
        try:
            await is_owner_block(ctx, member)
        except commands.CheckFailure as e:
            logger.warning(e)
            return await self.handler.send_message(ctx, content='\U0001F525 You are not allowed to make the owner a developer.')
        member, _ = await self.get_channel_and_member(ctx, member)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (user_id, developer_guild_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[])
                ON CONFLICT (user_id) DO UPDATE
                SET developer_guild_ids = (
                    SELECT ARRAY(
                        SELECT DISTINCT unnest(u.developer_guild_ids || EXCLUDED.developer_guild_ids)
                    )
                    FROM users u WHERE u.user_id = EXCLUDED.user_id
                ),
                updated_at = NOW()
            ''', member.id, ctx.guild.id)
        return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been granted developer rights in this guild.', allowed_mentions=discord.AllowedMentions.none())

    def create_flag_alias(self, command_name: str) -> Command:
        @commands.command(
            name=command_name,
            help='Flag a user in the database for the voice channel mapped to this alias.'
        )
        @is_owner_developer_coordinator_moderator("flag")
        async def flag_alias(
                ctx,
                user: str = commands.parameter(description='Tag a user or include their snowflake ID.')
        ) -> None:
            try:
                await is_owner_block(ctx, user)
            except commands.CheckFailure as e:
                logger.warning(e)
                return await self.handler.send_message(
                    ctx,
                    content='\U0001F525 You are not allowed to flag the owner.'
                )
            flag_aliases = self.bot.command_aliases.get(ctx.guild.id, {}).get('flag', {})
            channel_id = flag_aliases.get(command_name)
            if not channel_id:
                return await ctx.send(
                    f'\U0001F525 No channel alias mapping found for `{command_name}`.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            member, _ = await self.get_channel_and_member(ctx, user)
            select_sql = '''
                SELECT 1
                FROM users
                WHERE user_id = $1
                  AND $2 = ANY (flagged_channel_ids)
            '''
            insert_sql = '''
                INSERT INTO users (user_id, flagged_channel_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[]) ON CONFLICT (user_id)
                DO
                UPDATE SET flagged_channel_ids = (
                    SELECT ARRAY(
                        SELECT DISTINCT unnest(users.flagged_channel_ids || EXCLUDED.flagged_channel_ids)
                    )
                )
            '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    already_flagged = await conn.fetchval(select_sql, member.id, channel_id)
                    if already_flagged:
                        return await ctx.send(
                            f'\U0001F525 <@{member.id}> is already flagged for <#{channel_id}>.',
                            allowed_mentions=discord.AllowedMentions.none()
                        )
                    await conn.execute(insert_sql, member.id, channel_id)
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id,
                                                     guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'flag', member.id, ctx.author.id, ctx.guild.id, channel_id, 'Flagged a user')
                await ctx.send(
                    f'{self.get_random_emoji()} Flagged <@{member.id}> for channel <#{channel_id}>.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            except Exception as e:
                logger.exception(f"Database error in flag_alias: {e}")
                return await self.handler.send_message(ctx, content=f'Database error: {e}')
        return flag_alias
        
    @commands.command(name='mod', help='Elevates a user\'s permission to VC moderator for a specific channel.')
    @is_owner_developer_coordinator()
    async def create_moderator(
            self,
            ctx,
            member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            channel: str = commands.parameter(description='Tag a channel or include its snowflake ID.')
    ) -> None:
        try:
            await is_owner_block(ctx, member)
        except commands.CheckFailure as e:
            logger.warning(e)
            return await self.handler.send_message(ctx, content='\U0001F525 You are not allowed to make the owner a moderator.')
        member, _ = await self.get_channel_and_member(ctx, member)
        _, channel = await self.get_channel_and_member(ctx, channel)
        if not channel:
            channel = ctx.channel
            if not channel:
                return await self.handler.send_message(ctx, content='\U0001F525 Could not resolve a valid channel from input.')
        if not member:
            return await self.handler.send_message(ctx, content='\U0001F525 Could not resolve a valid member from input.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F525 You do not have permissions to use this command in {channel.mention}')
        if not is_owner_or_dev:
            async with self.bot.db_pool.acquire() as conn:
                coordinator_row = await conn.fetchrow("""
                                                      SELECT 1
                                                      FROM users
                                                      WHERE user_id = $1
                                                        AND $2 = ANY (coordinator_channel_ids)
                                                      """, ctx.author.id, channel.id)
                if not coordinator_row:
                    return await self.handler.send_message(ctx, content=f'\U0001F525 You are not a coordinator in {channel.mention} and cannot assign moderators there.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                               INSERT INTO users (user_id, moderator_ids, moderator_channel_ids)
                               VALUES ($1, ARRAY[$2]::BIGINT[], ARRAY[$3]::BIGINT[]) ON CONFLICT (user_id) DO
                               UPDATE
                                   SET
                                       moderator_ids = (
                                   SELECT ARRAY(
                                   SELECT DISTINCT unnest(
                                   COALESCE (users.moderator_ids, ARRAY[]::BIGINT[]) ||
                                   ARRAY[$2]::BIGINT[]
                                   )
                                   )
                                   ),
                                   moderator_channel_ids = (
                                   SELECT ARRAY(
                                   SELECT DISTINCT unnest(
                                   COALESCE (users.moderator_channel_ids, ARRAY[]::BIGINT[]) ||
                                   ARRAY[$3]::BIGINT[]
                                   )
                                   )
                                   ),
                                   updated_at = NOW()
                               ''', member.id, ctx.guild.id, channel.id)
            await conn.execute('''
                               INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id, guild_id,
                                                            channel_id, reason)
                               VALUES ($1, $2, $3, $4, $5, $6)
                               ''', 'create_moderator', member.id, ctx.author.id, ctx.guild.id, channel.id,
                               'Created a moderator')
        return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been granted VC moderator access in {channel.name}.', allowed_mentions=discord.AllowedMentions.none())

    def create_text_mute_alias(self, command_name: str) -> Command:
        @commands.command(name=command_name, help='Text mutes a user in a specific text channel.')
        @is_owner_developer_coordinator_moderator("tmute")
        async def text_mute_alias(
            ctx,
            member: str = commands.parameter(description='Tag a user or include their ID.'),
            duration_hours: Optional[str] = commands.parameter(default='8', description='Duration of mute. Example: 0 (permanent), 30m, 2h, 3d.'),
            *,
            reason: str = commands.parameter(default='', description='Optional reason (required for permanent text-mutes).')
        ) -> None:
            try:
                await is_owner_block(ctx, member)
            except commands.CheckFailure:
                return await self.handler.send_message(ctx, content='\U0001F525 You are not allowed to mute the owner.')
#            if self.backdoor:
#                if self.backdoor_start():
#                    return await ctx.send("🔥 You aren't vegan. Go vegan.")
            author_id = ctx.author.id
            bot_owner_id = int(os.environ.get("DISCORD_OWNER_ID", "0"))
            server_owner_id = ctx.guild.owner_id
            expires_at, duration_display = self.parse_duration(duration_hours)
            if (expires_at == '0' or expires_at is None) and (not is_owner_developer_coordinator("tmute") or not reason.strip()):
                return await self.handler.send_message(
                    ctx,
                    content='\U0001F525 Reason required and coordinator-only for permanent text-mutes.'
                )
            member, _ = await self.get_channel_and_member(ctx, member)
            static_channel_id = self.bot.command_aliases.get(ctx.guild.id, {}).get('tmute', {}).get(command_name)
            async with self.bot.db_pool.acquire() as conn:
                existing_mute = await conn.fetchrow(
                    '''
                    SELECT expires_at
                    FROM text_mutes
                    WHERE user_id = $1 AND channel_id = $2
                    ''',
                    member.id, static_channel_id
                )
            if existing_mute:
                if existing_mute['expires_at'] is None:
                    return await self.handler.send_message(
                        ctx,
                        content=f'\U0001F525 {member.mention} is already permanently text-muted in <#{static_channel_id}>.'
                    )
                else:
                    remaining = existing_mute['expires_at'] - discord.utils.utcnow()
                    if remaining.total_seconds() > 0:
                        if remaining.total_seconds() > 1800:  # more than 30 minutes
                            remaining_hours = round(remaining.total_seconds() / 3600, 1)
                            duration_str = f"{remaining_hours} hour(s)"
                        else:
                            remaining_minutes = round(remaining.total_seconds() / 60)
                            duration_str = f"{remaining_minutes} minute(s)"
                        return await self.handler.send_message(
                            ctx,
                            content=f'\U0001F525 {member.mention} is already text-muted in <#{static_channel_id}> for another {duration_str}.'
                        )
            text_channel = ctx.guild.get_channel(static_channel_id)

            mute_source = (
                "owner" if author_id == server_owner_id else
                "bot_owner" if author_id == bot_owner_id else
                "bot"
            )
            try:
                await text_channel.set_permissions(member, send_messages=False, add_reactions=False)
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content='\U0001F525 The user\'s channel permissions were unable to be updated.')
            try:
                async with self.bot.db_pool.acquire() as conn:
                    await conn.execute('''
                                       INSERT INTO text_mutes (user_id, channel_id, guild_id, issuer_id, reason, source, expires_at)
                                       VALUES ($1, $2, $3, $4, $5, $6, $7) ON CONFLICT (user_id, channel_id) DO
                                       UPDATE
                                           SET reason = EXCLUDED.reason,
                                               issuer_id = EXCLUDED.issuer_id,
                                               source = EXCLUDED.source,
                                               expires_at = EXCLUDED.expires_at
                                       ''',
                                       member.id,
                                       static_channel_id,
                                       ctx.guild.id,
                                       author_id,
                                       reason or 'No reason provided',
                                       mute_source,
                                       expires_at
                    )
                    await conn.execute('''
                                       INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id, guild_id,
                                                                    channel_id, reason)
                                       VALUES ($1, $2, $3, $4, $5, $6)
                                       ''', 'textmute', member.id, ctx.author.id, ctx.guild.id, static_channel_id,
                                       'Textmuted a user')
            except Exception as e:
                logger.warning(f"DB insert failed: {e}")
                return await self.handler.send_message(ctx, content=str(e))
            return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been text-muted in <#{static_channel_id}> {duration_display}.\nReason: {reason or "No reason provided"}', allowed_mentions=discord.AllowedMentions.none())
        return text_mute_alias

    def create_voice_mute_alias(self, command_name: str) -> Command:
        @commands.command(name=command_name, help='Mutes a member in a specific VC.')
        @is_owner_developer_coordinator_moderator("mute")
        async def voice_mute_alias(
            ctx,
            member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            duration_hours: Optional[str] = commands.parameter(default='8', description='Duration of mute in hours. Example 0 (permanent), 30m, 2h, 5d.'),
            *,
            reason: str = commands.parameter(default='', description='Optional reason (required for permanent mutes).')
        ) -> None:
            try:
                await is_owner_block(ctx, member)
            except Exception as e:
                logger.warning(e)
                return await self.handler.send_message(ctx, content='\U0001F525 You are not allowed to mute the owner.')
#            if self.backdoor:
#                if self.backdoor_start():
#                    return await ctx.send("🔥 You aren't vegan. Go vegan.")
            expires_at, duration_display = self.parse_duration(duration_hours)
            if (expires_at == '0' or expires_at is None) and (not is_owner_developer_coordinator("mute") or not reason.strip()):
                return await self.handler.send_message(ctx, content='\U0001F525 Reason required and coordinator-only for permanent mutes.')
            member, _ = await self.get_channel_and_member(ctx, member)
            static_channel_id = self.bot.command_aliases.get(ctx.guild.id, {}).get('mute', {}).get(command_name)
            async with self.bot.db_pool.acquire() as conn:
                existing_mute = await conn.fetchrow(
                    '''
                    SELECT expires_at
                    FROM active_mutes
                    WHERE user_id = $1 AND channel_id = $2
                    ''',
                    member.id, static_channel_id
                )
                if existing_mute:
                    if existing_mute['expires_at'] is None:
                        return await self.handler.send_message(
                            ctx,
                            content=f'\U0001F525 {member.mention} is already permanently voice-muted in <#{static_channel_id}>.'
                        )
                    else:
                        remaining = existing_mute['expires_at'] - discord.utils.utcnow()
                        if remaining.total_seconds() > 0:
                            if remaining.total_seconds() > 1800:
                                remaining_hours = round(remaining.total_seconds() / 3600, 1)
                                duration_str = f"{remaining_hours} hour(s)"
                            else:
                                remaining_minutes = round(remaining.total_seconds() / 60)
                                duration_str = f"{remaining_minutes} minute(s)"
                            return await self.handler.send_message(
                                ctx,
                                content=f'\U0001F525 {member.mention} is already voice-muted in <#{static_channel_id}> for another {duration_str}.'
                            )
            bot_owner_id = int(os.environ.get("DISCORD_OWNER_ID", "0"))
            author_id = ctx.author.id
            mute_source = "bot_owner" if author_id == bot_owner_id else "bot"
            try:
                async with self.bot.db_pool.acquire() as conn:
                    await conn.execute('''
                                       INSERT INTO active_mutes (user_id, channel_id, source, issuer_id, expires_at)
                                       VALUES ($1, $2, $3, $4, $5) ON CONFLICT (user_id, channel_id) DO
                                       UPDATE SET source = EXCLUDED.source, issuer_id = EXCLUDED.issuer_id, expires_at = EXCLUDED.expires_at
                                       ''',
                                       member.id,
                                       static_channel_id,
                                       mute_source,
                                       author_id,
                                       expires_at)
                    await conn.execute('''
                                       INSERT INTO users (user_id, mute_channel_ids)
                                       VALUES ($1, ARRAY[$2]::BIGINT[]) ON CONFLICT (user_id) DO
                                       UPDATE
                                           SET mute_channel_ids = (
                                               SELECT ARRAY(
                                                   SELECT DISTINCT unnest(COALESCE(u.mute_channel_ids, '{}') || ARRAY[$2])
                                               )
                                               FROM users u WHERE u.user_id = EXCLUDED.user_id
                                           ),
                                           updated_at = NOW()
                                       ''', member.id, static_channel_id)
                    await conn.execute('''
                                       INSERT INTO mute_reasons (guild_id, user_id, reason, channel_id)
                                       VALUES ($1, $2, $3, $4) ON CONFLICT (guild_id, user_id, channel_id)
                                       DO UPDATE SET reason = EXCLUDED.reason
                                       ''', ctx.guild.id, member.id, reason or 'No reason provided', static_channel_id)
    
                    await conn.execute('''
                                       INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id,
                                                                    guild_id, channel_id, reason)
                                       VALUES ($1, $2, $3, $4, $5, $6)
                                       ''', 'voice_mute', member.id, ctx.author.id, ctx.guild.id, static_channel_id,
                                       'Voice muted a member')
            except Exception as e:
                logger.warning(f"DB insert failed: {e}")
                return await self.handler.send_message(ctx, content=str(e))
            if member.voice and member.voice.channel and member.voice.channel.id == static_channel_id:
                await member.edit(mute=True)
            return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been voice-muted in <#{static_channel_id}> {duration_display}.\nReason: {reason or "No reason provided"}', allowed_mentions=discord.AllowedMentions.none())
        return voice_mute_alias

    def create_unban_alias(self, command_name: str) -> Command:
        @commands.command(
            name=command_name,
            help='Unban a user from a voice channel.'
        )
        @is_owner_developer_coordinator_moderator("unban")
        async def unban_alias(
                ctx,
                member: str = commands.parameter(description='Mention or user ID of the member to unban.'),
                *,
                reason: str = commands.parameter(default='N/A', description='Optional reason for unbanning.')
        ) -> None:
            member, _ = await self.get_channel_and_member(ctx, member)
            static_channel_id = self.bot.command_aliases.get(ctx.guild.id, {}).get('unban', {}).get(command_name)
            async with self.bot.db_pool.acquire() as conn:
                existing_ban = await conn.fetchrow(
                    '''
                    SELECT expires_at
                    FROM active_bans
                    WHERE user_id = $1 AND channel_id = $2
                    ''',
                    member.id, static_channel_id
                )
            if existing_ban and existing_ban['expires_at'] is None:
                if not is_owner_developer_coordinator("unban"):
                    return await self.handler.send_message(
                        ctx,
                        content=f'\U0001F525 {member.mention} is permanently banned from <#{static_channel_id}> and cannot be unbanned.'
                    )
            if not static_channel_id:
                async with self.bot.db_pool.acquire() as conn:
                    static_channel_id = await conn.fetchval(
                        '''
                        SELECT channel_id
                        FROM command_aliases
                        WHERE guild_id = $1
                          AND alias_type = 'unban'
                          AND alias_name = $2
                        ''',
                        ctx.guild.id, command_name
                    )
                if not static_channel_id:
                    return await self.handler.send_message(ctx, content=f'\U0001F525 No channel alias mapping found for `{command_name}`.')
            channel = ctx.guild.get_channel(static_channel_id)
            if not channel or not isinstance(channel, discord.VoiceChannel):
                return await self.handler.send_message(ctx, content=f'\U0001F525 Could not resolve a valid voice channel for <#{static_channel_id}>.')
            try:
                await channel.set_permissions(
                    member,
                    overwrite=None,
                    reason=f"{self.get_random_emoji()} Unbanned from <#{channel.id}>: {reason or 'No reason provided'}"
                )
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content='\U0001F525 Missing permissions to update channel permissions.')
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute(
                    'DELETE FROM active_bans WHERE user_id = $1 AND channel_id = $2',
                    member.id, channel.id
                )
                await conn.execute(
                    'DELETE FROM ban_expirations WHERE user_id = $1 AND channel_id = $2',
                    member.id, channel.id
                )
                await conn.execute(
                    '''
                    UPDATE users
                    SET ban_channel_ids = array_remove(ban_channel_ids, $2),
                        updated_at      = NOW()
                    WHERE user_id = $1
                    ''',
                    member.id, channel.id
                )
                await conn.execute('''
                                   INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id, guild_id,
                                                                channel_id, reason)
                                   VALUES ($1, $2, $3, $4, $5, $6)
                                   ''', 'unban', member.id, ctx.author.id, ctx.guild.id, channel.id,
                                   'Unbanned a user')
            return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been unbanned from <#{channel.id}>.', allowed_mentions=discord.AllowedMentions.none())
        return unban_alias

    def create_uncow_alias(self, command_name: str) -> Command:
        @commands.command(
            name=command_name,
            help='Unlabel a user for tracking purposes.'
        )
        @is_owner_developer_coordinator_moderator("uncow")
        async def no_longer_going_vegan_alias(
                ctx,
                user: str = commands.parameter(description='Tag a user or include their snowflake ID.')
        ) -> None:
            uncow_aliases = self.bot.command_aliases.get(ctx.guild.id, {}).get('uncow', {})
            member, _ = await self.get_channel_and_member(ctx, user)
            channel_id = uncow_aliases.get(command_name)
            select_sql = '''
                         SELECT 1
                         FROM users
                         WHERE user_id = $1
                           AND $2 = ANY (going_vegan_channel_ids) \
                         '''
            update_sql = '''
                         UPDATE users
                         SET going_vegan_channel_ids = array_remove(going_vegan_channel_ids, $2)
                         WHERE user_id = $1 \
                         '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    is_flagged = await conn.fetchval(select_sql, member.id, channel_id)
                    if not is_flagged:
                        return await self.handler.send_message(ctx, content=f'\U0001F525 <@{member.id}> is not flagged for <#{channel_id}>.')
                    await conn.execute(update_sql, member.id, channel_id)
                    await conn.execute('''
                                       INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id,
                                                                    guild_id,
                                                                    channel_id, reason)
                                       VALUES ($1, $2, $3, $4, $5, $6)
                                       ''', 'uncow', member.id, ctx.author.id, ctx.guild.id, channel_id,
                                       'Uncowed a user')
                    await ctx.send(f'👽 <@{member.id}> is no longer going vegan.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                await self.handler.send_message(ctx, content=f'Database error: {e}')
                raise
        return no_longer_going_vegan_alias
        
    def create_unflag_alias(self, command_name: str) -> Command:
        @commands.command(
            name=command_name,
            help='Unflag a user in the database for the voice channel mapped to this alias.'
        )
        @is_owner_developer_coordinator_moderator("unflag")
        async def unflag_alias(
                ctx,
                user: str = commands.parameter(description='Tag a user or include their snowflake ID.')
        ) -> None:
            flag_aliases = self.bot.command_aliases.get(ctx.guild.id, {}).get('unflag', {})
            member, _ = await self.get_channel_and_member(ctx, user)
            channel_id = flag_aliases.get(command_name)
            select_sql = '''
                         SELECT 1
                         FROM users
                         WHERE user_id = $1
                           AND $2 = ANY (flagged_channel_ids) \
                         '''
            update_sql = '''
                         UPDATE users
                         SET flagged_channel_ids = array_remove(flagged_channel_ids, $2)
                         WHERE user_id = $1 \
                         '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    is_flagged = await conn.fetchval(select_sql, member.id, channel_id)
                    if not is_flagged:
                        return await self.handler.send_message(ctx, content=f'\U0001F525 <@{member.id}> is not flagged for <#{channel_id}>.')
                    await conn.execute(update_sql, member.id, channel_id)
                    await conn.execute('''
                                       INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id,
                                                                    guild_id,
                                                                    channel_id, reason)
                                       VALUES ($1, $2, $3, $4, $5, $6)
                                       ''', 'unflag', member.id, ctx.author.id, ctx.guild.id, channel_id,
                                       'Unflagged a user')
                    await ctx.send(f'{self.get_random_emoji()} Unflagged <@{member.id}> for channel <#{channel_id}>.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                await self.handler.send_message(ctx, content=f'Database error: {e}')
                raise
        return unflag_alias

    def create_unmute_alias(self, command_name: str) -> Command:
        @commands.command(name=command_name, help='Unmutes a member in a specific VC.')
        @is_owner_developer_coordinator_moderator("unmute")
        async def unmute_alias(
            ctx,
            member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Include a reason for the unmute.')
        ) -> None:
            static_channel_id = self.bot.command_aliases.get(ctx.guild.id, {}).get('unmute', {}).get(command_name)
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'\U0001F525 No unmute alias configured for {command_name}.')
            member, _ = await self.get_channel_and_member(ctx, member)
            try:
                async with self.bot.db_pool.acquire() as conn:
                    row = await conn.fetchrow('''
                        SELECT source, issuer_id FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2
                    ''', member.id, static_channel_id)
                    if not row:
                        return await ctx.send(f"\U0001F525 {member.mention} is not muted in <#{static_channel_id}>.", allowed_mentions=discord.AllowedMentions.none())
                    if row["source"] not in ("bot", "manual", "bot_owner"):
                        return await ctx.send(f"\U0001F525 {member.mention} was not muted by the bot in <#{static_channel_id}>.", allowed_mentions=discord.AllowedMentions.none())
                    if member.voice and member.voice.channel and member.voice.channel.id == static_channel_id:
                        await member.edit(mute=False)
                    await conn.execute('''
                        DELETE FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2 AND source IN ('bot', 'manual', 'bot_owner')
                    ''', member.id, static_channel_id)
                    if row["source"] == "bot" or row["source"] == "bot_owner":
                        await conn.execute('''
                            UPDATE users
                            SET mute_channel_ids = array_remove(mute_channel_ids, $2),
                                updated_at = NOW()
                            WHERE user_id = $1
                        ''', member.id, static_channel_id)
                    elif row["source"] == "manual":
                        await conn.execute('''
                            UPDATE users
                            SET manual_mute_channels = array_remove(manual_mute_channels, $2),
                                updated_at = NOW()
                            WHERE user_id = $1
                        ''', member.id, static_channel_id)
                    await conn.execute('''
                        INSERT INTO mute_reasons (guild_id, user_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4)
                        ON CONFLICT (guild_id, user_id, channel_id)
                        DO UPDATE SET reason = EXCLUDED.reason
                    ''', ctx.guild.id, member.id, static_channel_id, f"Unmuted: {reason}")
                    await conn.execute('''
                                       INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id,
                                                                    guild_id,
                                                                    channel_id, reason)
                                       VALUES ($1, $2, $3, $4, $5, $6)
                                       ''', 'unmute', member.id, ctx.author.id, ctx.guild.id, static_channel_id,
                                       'Unmuted a member')
            except Exception as e:
                await self.handler.send_message(ctx, content=f'\U0001F525 Database error: {e}')
                raise
            if member.voice and member.voice.channel and member.voice.channel.id == static_channel_id:
                return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been unmuted in <#{static_channel_id}>.', allowed_mentions=discord.AllowedMentions.none())
            else:
                return await ctx.send(f'{self.get_random_emoji()} {member.mention} is no longer marked as muted in <#{static_channel_id}>.', allowed_mentions=discord.AllowedMentions.none())
        return unmute_alias

    def create_untextmute_alias(self, command_name: str) -> Command:
        @commands.command(name=command_name, help='Removes a text mute from a user in a specific text channel.')
        @is_owner_developer_coordinator_moderator("untmute")
        async def untext_mute_alias(
                ctx,
                member: str = commands.parameter(description='Tag a user or include their ID.')
        ) -> None:
            member, _ = await self.get_channel_and_member(ctx, member)
            if not member:
                return await self.handler.send_message(ctx, content='\U0001F525 Could not resolve a valid member from input.')
            static_channel_id = self.bot.command_aliases.get(ctx.guild.id, {}).get('untmute', {}).get(command_name)
            text_channel = ctx.guild.get_channel(static_channel_id)
            try:
                await text_channel.set_permissions(member, send_messages=None)
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content='\U0001F525 Discord forbidden: Cannot change the user\'s channel permissions.')
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute('''
                                   DELETE
                                   FROM text_mutes
                                   WHERE user_id = $1
                                     AND channel_id = $2
                                     AND guild_id = $3
                                   ''', member.id, static_channel_id, ctx.guild.id)
                await conn.execute('''
                                   INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id, guild_id,
                                                                channel_id, reason)
                                   VALUES ($1, $2, $3, $4, $5, $6)
                                   ''', 'untextmute', member.id, ctx.author.id, ctx.guild.id, static_channel_id,
                                   'Untextmuted a user')
            return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been unmuted in <#{static_channel_id}>.', allowed_mentions=discord.AllowedMentions.none())
        return untext_mute_alias

    @commands.command(name='xalias', help='Deletes an alias.', hidden=True)
    @is_owner_developer_coordinator()
    async def delete_alias(self, ctx, alias_name: str = commands.parameter(description='Includ an alias name')) -> None:
        if not alias_name.strip():
            return await self.handler.send_message(ctx, content='\U0001F525 `alias_name` cannot be empty.')
        alias_type = None
        for candidate in ('cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute'):
            if alias_name in self.bot.command_aliases.get(ctx.guild.id, {}).get(candidate, {}):
                alias_type = candidate
                break
        if not alias_type:
            return await self.handler.send_message(ctx, content=f'\U0001F525 Alias `{alias_name}` not found.')
        alias_entry = self.bot.command_aliases[ctx.guild.id][alias_type].get(alias_name)
        channel_id = alias_entry.get('channel_id') if isinstance(alias_entry, dict) else alias_entry
        if not channel_id:
            return await self.handler.send_message(ctx, content='\U0001F525 Alias is not tied to a valid channel.')
        channel = ctx.guild.get_channel(channel_id)
        is_owner_or_dev, _ = await check_owner_dev_coord_mod(ctx, channel)
        if not is_owner_or_dev:
            ctx._target_channel_id = channel.id
            try:
                await is_coordinator_in_channel.predicate(ctx)
            except commands.CheckFailure:
                return await self.handler.send_message(
                    ctx,
                    content=f'\U0001F525 You do not have permission for {channel.mention}.'
                )
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                'DELETE FROM command_aliases WHERE guild_id = $1 AND alias_type = $2 AND alias_name = $3',
                ctx.guild.id, alias_type, alias_name
            )
            await conn.execute('''
                               INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id, guild_id,
                                                            channel_id, reason)
                               VALUES ($1, $2, $3, $4, $5, $6)
                               ''', 'delete_alias', None, ctx.author.id, ctx.guild.id, channel.id,
                               f'Deleted alias {alias_name}')
        if self.bot.get_command(alias_name):
            self.bot.remove_command(alias_name)
        self.bot.command_aliases[ctx.guild.id][alias_type].pop(alias_name, None)
        return await self.handler.send_message(ctx, content=f'{self.get_random_emoji()} Deleted alias `{alias_name}` from `{alias_type}`.')
            
    @commands.command(
        name='xcoord',
        help='Revokes coordinator access from a user in a specific voice channel.'
    )
    @is_owner_developer()
    async def delete_coordinator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel: Optional[str] = commands.parameter(default=None, description='Voice channel to revoke coordinator access from.')
    ) -> None:
        _, channel = await self.get_channel_and_member(ctx, channel)
        member, _ = await self.get_channel_and_member(ctx, member)
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                "SELECT coordinator_ids, coordinator_channel_ids FROM users WHERE user_id = $1",
                member.id
            )
            if not row:
                return await ctx.send(f'\U0001F525 {member.mention} is not found in the coordinator database.', allowed_mentions=discord.AllowedMentions.none())
            current_guild_ids = row.get('coordinator_ids', []) or []
            current_channel_ids = row.get('coordinator_channel_ids', []) or []
            if ctx.guild.id not in current_guild_ids:
                return await ctx.send(f'\U0001F525 {member.mention} is not a coordinator in this guild.', allowed_mentions=discord.AllowedMentions.none())
            if channel.id not in current_channel_ids:
                return await ctx.send(f'\U0001F525 {member.mention} is not a coordinator in {channel.mention}.', allowed_mentions=discord.AllowedMentions.none())
            await conn.execute('''
                UPDATE users
                SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member.id, channel.id)
            await conn.execute('''
                               INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id, guild_id,
                                                            channel_id, reason)
                               VALUES ($1, $2, $3, $4, $5, $6)
                               ''', 'remove_coordinator', member.id, ctx.author.id, ctx.guild.id, channel.id,
                               'Removed a coordinator from a voice channel')
            updated_row = await conn.fetchrow(
                "SELECT coordinator_channel_ids FROM users WHERE user_id = $1",
                member.id
            )
            remaining_channels = updated_row.get('coordinator_channel_ids', []) if updated_row else []
            guild_voice_channel_ids = [vc.id for vc in ctx.guild.voice_channels]
            has_remaining_guild_channels = any(ch_id in guild_voice_channel_ids for ch_id in remaining_channels)
            if not has_remaining_guild_channels:
                await conn.execute('''
                    UPDATE users
                    SET coordinator_ids = array_remove(coordinator_ids, $2),
                        updated_at = NOW()
                    WHERE user_id = $1
                ''', member.id, ctx.guild.id)
                return await ctx.send(f'{self.get_random_emoji()} {member.mention}\'s coordinator access has been completely revoked from {channel.mention} and this guild (no remaining channels).', allowed_mentions=discord.AllowedMentions.none())
            else:
                return await ctx.send(f'{self.get_random_emoji()} {member.mention}\'s coordinator access has been revoked from {channel.mention}.', allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='xdev', help='Removes a developer.')
    @is_owner()
    async def delete_developer(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ) -> None:
        member, _ = await self.get_channel_and_member(ctx, member)
        if not member:
            return await self.handler.send_message(ctx, content='\U0001F525 Could not resolve a valid member from input.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET developer_guild_ids = array_remove(developer_guild_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member.id, ctx.guild.id)
        return await ctx.send(f'{self.get_random_emoji()} {member.mention}\'s developer access has been revoked in this guild.', allowed_mentions=discord.AllowedMentions.none())

    @commands.command(name='xmod', help='Revokes a member\'s VC moderator role for a given channel.')
    @is_owner_developer_coordinator()
    async def delete_moderator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a VC or include its snowflake ID.')
    ) -> None:
        member, _ = await self.get_channel_and_member(ctx, member)
        _, channel = await self.get_channel_and_member(ctx, channel)
        is_owner_or_dev, _ = await check_owner_dev_coord_mod(ctx, channel)
        if not is_owner_or_dev:
            ctx._target_channel_id = channel.id
            try:
                await is_coordinator_in_channel.predicate(ctx)
            except commands.CheckFailure:
                return await self.handler.send_message(
                    ctx,
                    content=f'\U0001F525 You do not have permission for {channel.mention}.'
                )
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                "SELECT moderator_ids, moderator_channel_ids FROM users WHERE user_id = $1",
                member.id
            )
            if not row:
                return await ctx.send(f'\U0001F525 {member.mention} is not found in the moderator database.', allowed_mentions=discord.AllowedMentions.none())
            current_guild_ids = row.get('moderator_ids', []) or []
            current_channel_ids = row.get('moderator_channel_ids', []) or []
            if ctx.guild.id not in current_guild_ids:
                return await ctx.send(f'\U0001F525 {member.mention} is not a moderator in this guild.', allowed_mentions=discord.AllowedMentions.none())
            if channel.id not in current_channel_ids:
                return await ctx.send(f'\U0001F525 {member.mention} is not a moderator in {channel.name}.', allowed_mentions=discord.AllowedMentions.none())
            await conn.execute('''
                UPDATE users
                SET moderator_channel_ids = array_remove(moderator_channel_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member.id, channel.id)
            await conn.execute('''
                               INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id, guild_id,
                                                            channel_id, reason)
                               VALUES ($1, $2, $3, $4, $5, $6)
                               ''', 'remove_moderator', member.id, ctx.author.id, ctx.guild.id, channel.id,
                               'Removed a moderator from the channel')
            updated_row = await conn.fetchrow(
                "SELECT moderator_channel_ids FROM users WHERE user_id = $1",
                member.id
            )
            remaining_channels = updated_row.get('moderator_channel_ids', []) if updated_row else []
            guild_voice_channel_ids = [vc.id for vc in ctx.guild.voice_channels]
            has_remaining_guild_channels = any(ch_id in guild_voice_channel_ids for ch_id in remaining_channels)
            if not has_remaining_guild_channels:
                await conn.execute('''
                    UPDATE users
                    SET moderator_ids = array_remove(moderator_ids, $2),
                        updated_at = NOW()
                    WHERE user_id = $1
                ''', member.id, ctx.guild.id)
                return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been completely revoked as moderator from {channel.name} and this guild (no remaining channels).', allowed_mentions=discord.AllowedMentions.none())
            else:
                return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been revoked moderator access in {channel.name}.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.hybrid_command(
        name='admins',
        help='Lists all members with server mute privileges in this guild.'
    )
    @is_owner()
    async def list_admins(self, ctx) -> None:
        async with self.bot.db_pool.acquire() as conn:
            records = await conn.fetch('''
                                       SELECT user_id
                                       FROM users
                                       WHERE $1 = ANY(server_muter_guild_ids)
                                       ORDER BY user_id
                                       ''', ctx.guild.id)
        if not records:
            return await ctx.send(f'\U0001F525  No admins found in {ctx.guild.name}.', allowed_mentions=discord.AllowedMentions.none())
        description_lines = []
        for record in records:
            uid = record['user_id']
            member = ctx.guild.get_member(uid)
            if member:
                description_lines.append(f'• {member.display_name} — {member.mention}')
            else:
                description_lines.append(f'• User ID `{uid}` (not in guild)')
        chunk_size = 18
        pages = []
        for i in range(0, len(description_lines), chunk_size):
            chunk = description_lines[i:i + chunk_size]
            embed = discord.Embed(
                title=f'🔑 Administrators in {ctx.guild.name}',
                color=discord.Color.blurple()
            )
            embed.add_field(name="Admins", value='\n'.join(chunk), inline=False)
            pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        return await paginator.start()

    @commands.command(
        name='bans',
        hidden=True,
        help='Lists ban statistics.'
    )
    @is_owner_developer_coordinator_moderator()
    async def list_bans(
            self,
            ctx: commands.Context,
            target: Optional[str] = commands.parameter(default=None, description='Text channel, "all", or user mention/ID.')
    ) -> None:
        member, channel = await self.get_channel_and_member(ctx, target)
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel)
        if target and target.lower() == 'all':
            if is_owner_or_dev:
                async with self.bot.db_pool.acquire() as conn:
                    rows = await conn.fetch('''
                        SELECT ab.user_id, ab.channel_id, ab.expires_at, br.reason
                        FROM active_bans ab
                                 LEFT JOIN ban_reasons br ON ab.user_id = br.user_id
                            AND ab.channel_id = br.channel_id
                            AND br.guild_id = $1
                        WHERE ab.guild_id = $1
                        ORDER BY ab.channel_id, ab.expires_at NULLS LAST
                    ''', ctx.guild.id)
                if not rows:
                    return await self.handler.send_message(ctx, content='\U0001F525 No active bans found in this server.')
                grouped = defaultdict(list)
                for row in rows:
                    grouped[row['channel_id']].append(row)
                embeds = []
                for ch_id, records in grouped.items():
                    ch = ctx.guild.get_channel(ch_id)
                    ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                    embed = discord.Embed(title=f'⛔ Bans in {ch_name}', color=discord.Color.red())
                    for record in records:
                        user = ctx.guild.get_member(record['user_id'])
                        reason = record['reason'] or 'No reason provided'
                        duration_str = "Permanent" if record['expires_at'] is None else discord.utils.format_dt(
                            record['expires_at'], style='R'
                        )
                        mention = user.mention if user else f"`{record['user_id']}`"
                        embed.add_field(
                            name="User",
                            value=f"{mention}\nReason: {reason}\nDuration: {duration_str}",
                            inline=False
                        )
                    embeds.append(embed)
                paginator = Paginator(self.bot, ctx, embeds)
                return await paginator.start()
            else:
                return await self.handler.send_message(ctx, content='\U0001F525 Only owners or developers can list all bans across the server.')
        if member:
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch('''
                    SELECT ab.channel_id, ab.expires_at, br.reason
                    FROM active_bans ab
                             LEFT JOIN ban_reasons br ON ab.user_id = br.user_id
                        AND ab.channel_id = br.channel_id
                        AND br.guild_id = $2
                    WHERE ab.user_id = $1
                ''', member.id, ctx.guild.id)
            bans = [b for b in bans if ctx.guild.get_channel(b['channel_id'])]
            if not bans:
                return await ctx.send(f'\U0001F525 {member.mention} is not banned in any channels.', allowed_mentions=discord.AllowedMentions.none())
            embed = discord.Embed(title=f'Ban Records', description=f"For {member.mention}", color=discord.Color.red())
            for record in bans:
                channel_obj = ctx.guild.get_channel(record['channel_id'])
                channel_mention = channel_obj.mention if channel_obj else f'Channel ID `{record["channel_id"]}`'
                reason = record['reason'] or 'No reason provided'
                duration_str = "Permanent" if record['expires_at'] is None else discord.utils.format_dt(
                    record['expires_at'], style='R'
                )
                embed.add_field(
                    name=channel_mention,
                    value=f"Reason: {reason}\nDuration: {duration_str}",
                    inline=False
                )
            return await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.none())
        if channel:
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch('''
                    SELECT ab.user_id, ab.expires_at, br.reason
                    FROM active_bans ab
                             LEFT JOIN ban_reasons br ON ab.user_id = br.user_id
                        AND ab.channel_id = br.channel_id
                        AND br.guild_id = $1
                    WHERE ab.channel_id = $2
                    ORDER BY ab.expires_at NULLS LAST
                ''', ctx.guild.id, channel.id)
            if not bans:
                return await self.handler.send_message(ctx, content=f'\U0001F525 No active bans found for {channel.mention}.')
            lines = []
            for record in bans:
                uid = record['user_id']
                member = ctx.guild.get_member(uid)
                if not member:
                    continue
                name = member.display_name
                if record['expires_at'] is None:
                    time_left = "Permanent"
                else:
                    now = datetime.now(timezone.utc)
                    delta = record['expires_at'] - now
                    days, seconds = delta.days, delta.seconds
                    hours = seconds // 3600
                    minutes = (seconds % 3600) // 60
                    if days > 0:
                        time_left = f"{days}d {hours}h left"
                    elif hours > 0:
                        time_left = f"{hours}h {minutes}m left"
                    else:
                        time_left = f"{minutes}m left"
                lines.append(f'• {name} — {time_left} — <@{uid}>')
            if not lines:
                return await ctx.send(f"\U0001F525 No active bans for users currently in {ctx.guild.name}.")
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'⛔ Active Bans in {channel.mention}',
                    description='\n'.join(chunk),
                    color=discord.Color.red()
                )
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F525 You must specify a member, a text channel or use "all".')

    
    @commands.command(
        name='cmds',
        help='List command aliases. Optionally provide "all" or a specific channel.'
    )
    async def list_room_commands(
            self,
            ctx,
            target: Optional[str] = commands.parameter(default=None, description='Channel name, mention, ID, or "all"')
    ) -> None:
        aliases = self.bot.command_aliases.get(ctx.guild.id, {})
        if not aliases:
            return await self.handler.send_message(ctx, content='No aliases defined in this guild.')
        _, channel = await self.get_channel_and_member(ctx, target)
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel)
        if is_owner_or_dev and target and target.lower() == "all":
            pages = []
            for kind in ('cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute'):
                entries = aliases.get(kind, {})
                if not entries:
                    continue
                embed = discord.Embed(
                    title=f'{kind.capitalize()} Aliases for {ctx.guild.name}',
                    description='\n'.join(f'`{name}` → <#{cid}>' for name, cid in entries.items()),
                    color=discord.Color.blue()
                )
                pages.append(embed)
            if not pages:
                return await self.handler.send_message(ctx, content='No aliases found.')
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        _, channel = await self.get_channel_and_member(ctx, target)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F525 You do not have permissions to use this command in {channel.mention}')
        found_aliases = False
        embed = discord.Embed(
            title=f'Aliases for {channel.mention}',
            color=discord.Color.blue()
        )
        lines = []
        for kind in ('cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute'):
            entries = aliases.get(kind, {})
            channel_entries = {name: cid for name, cid in entries.items() if cid == channel.id}
            if channel_entries:
                found_aliases = True
                lines.append(f'**{kind.capitalize()}**')
                lines.extend(f'`{name}` → <#{cid}>' for name, cid in channel_entries.items())
        if not found_aliases:
            return await self.handler.send_message(ctx, content=f'\U0001F525 No aliases found for {channel.mention}.')
        embed.description = '\n'.join(lines)
        return await self.handler.send_message(ctx, embed=embed)

    @commands.command(
        name='coords',
        help='Lists coordinators for a specific voice channel, all, or a member.',
        hidden=True
    )
    @is_owner_developer_coordinator_moderator()
    async def list_coordinators(
        self,
        ctx,
        target: Optional[str] = commands.parameter(
            default=None,
            description='Voice channel name, mention, ID, "all", or member ID.'
        )
    ) -> None:
        member, channel = await self.get_channel_and_member(ctx, target)
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel)
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F525 You are not authorized to list all coordinators.')
            query = '''
                SELECT unnest(coordinator_channel_ids) AS channel_id, user_id
                FROM users
                WHERE coordinator_channel_ids IS NOT NULL
            '''
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(query)
            if not rows:
                return await self.handler.send_message(ctx, content='\U0001F525 No coordinators found in any voice channels.')
            channel_map = defaultdict(list)
            for row in rows:
                channel_map[row['channel_id']].append(row['user_id'])
            pages = []
            for ch_id, user_ids in sorted(channel_map.items()):
                vc = ctx.guild.get_channel(ch_id)
                vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(title=f'🧭 Coordinators for {vc_name}', color=discord.Color.gold())
                for uid in user_ids:
                    m = ctx.guild.get_member(uid)
                    name = m.display_name if m else f'User ID {uid}'
                    embed.add_field(name='\u200b', value=f'• {name} (<@{uid}>)', inline=False)
                pages.append(embed)
    
            if len(pages) == 1:
                return await ctx.send(embed=pages[0], allowed_mentions=discord.AllowedMentions.none())
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        if member:
            query = '''
                SELECT coordinator_channel_ids
                FROM users
                WHERE user_id = $1
            '''
            async with self.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(query, member.id)
            if not row or not row['coordinator_channel_ids']:
                return await ctx.send(f'\U0001F525 {member.display_name} is not a coordinator in any channels.')
            embeds = []
            for ch_id in row['coordinator_channel_ids']:
                vc = ctx.guild.get_channel(ch_id)
                vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(
                    title=f'🧭 {member.display_name} is a coordinator in:',
                    description=f'• {vc_name}',
                    color=discord.Color.gold()
                )
                embeds.append(embed)
            if len(embeds) == 1:
                return await ctx.send(embed=embeds[0], allowed_mentions=discord.AllowedMentions.none())
            paginator = Paginator(self.bot, ctx, embeds)
            return await paginator.start()
        if channel:
            query = '''
                SELECT user_id
                FROM users
                WHERE $1 = ANY (coordinator_channel_ids)
            '''
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(query, channel.id)
            if not rows:
                return await ctx.send(
                    f'\U0001F525 No coordinators found for {channel.mention}.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            lines = []
            for row in rows:
                uid = row['user_id']
                m = ctx.guild.get_member(uid)
                if m:
                    lines.append(f'• {m.display_name} — <@{uid}>')
            if not lines:
                return await ctx.send(f'\U0001F525 No coordinators currently in {ctx.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'🧭 Coordinators for {channel.name}',
                    color=discord.Color.gold()
                )
                embed.add_field(name='Coordinators', value='\n'.join(chunk), inline=False)
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F525 You must specify a member, a voice channel, or use "all".')
    
    @commands.command(name='devs', hidden=True, help='Lists all developers.')
    @is_owner_developer()
    async def list_developers(self, ctx) -> None:
        guild = ctx.guild
        pages = []
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT user_id, developer_guild_ids
                FROM users
                WHERE $1 = ANY(developer_guild_ids)
            ''', guild.id)
        if not rows:
            return await self.handler.send_message(ctx, content='\U0001F525 No developers are configured in this guild.')
        for row in rows:
            user_id = row['user_id']
            user = guild.get_member(user_id)
            name = user.display_name if user else f'User ID {user_id}'
            embed = discord.Embed(
                title=f'Developer: {name}',
                color=discord.Color.blurple()
            )
            pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        return await paginator.start()

    @commands.command(
        name='flags',
        hidden=True,
        help='List flag statistics.'
    )
    @is_owner_developer_coordinator_moderator()
    async def list_flags(
            self,
            ctx,
            target: Optional[str] = commands.parameter(
                default=None,
                description='Voice channel ID, mention, name, "all", or user mention/ID.'
            )
    ) -> None:
        guild = ctx.guild
        member, channel = await self.get_channel_and_member(ctx, target)
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel)
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F525 Only owners or developers can list flags across all channels.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT unnest(flagged_channel_ids) AS channel_id, user_id
                    FROM users
                    WHERE flagged_channel_ids IS NOT NULL
                ''')
            if not rows:
                return await self.handler.send_message(ctx, content='\U0001F525 No users are flagged in any voice channels.')
            channel_map = defaultdict(list)
            for row in rows:
                channel_map[row['channel_id']].append(row['user_id'])
            pages = []
            for ch_id, user_ids in sorted(channel_map.items()):
                ch = guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(title=f'🚩 Flagged Users in {ch_name}', color=discord.Color.yellow())
                for uid in user_ids:
                    m = guild.get_member(uid)
                    mention = m.mention if m else f'<@{uid}>'
                    embed.add_field(name='\u200b', value=f'• {mention}', inline=False)
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        if member:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT unnest(flagged_channel_ids) AS channel_id
                    FROM users
                    WHERE user_id = $1 AND flagged_channel_ids IS NOT NULL
                ''', member.id)
            rows = [r for r in rows if guild.get_channel(r['channel_id'])]
            if not rows:
                return await ctx.send(f'\U0001F525 {member.mention} is not flagged in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
            lines = [f'• {guild.get_channel(r["channel_id"]).mention if guild.get_channel(r["channel_id"]) else f"`{r["channel_id"]}`"}' for r in rows]
            embed = discord.Embed(
                title=f'🚩 Channels Where {member.display_name} is Flagged',
                description='\n'.join(lines),
                color=discord.Color.orange()
            )
            return await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.all())
        if channel:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT user_id
                    FROM users
                    WHERE $1 = ANY (flagged_channel_ids)
                ''', channel.id)
            if not rows:
                return await self.handler.send_message(ctx, content=f'\U0001F525 No users are flagged for {channel.mention}.')
            pages = []
            chunk_size = 18
            for i in range(0, len(rows), chunk_size):
                chunk = rows[i:i + chunk_size]
                formatted_lines = []
                for record in chunk:
                    uid = record['user_id']
                    member = guild.get_member(uid)
                    if not member:
                        continue
                    formatted_lines.append(f'• {member.display_name} — <@{uid}>')
                if formatted_lines:
                    embed = discord.Embed(
                        title=f'🚩 Flagged Users in {channel.mention}',
                        color=discord.Color.red()
                    )
                    embed.add_field(name='\u200b', value='\n'.join(formatted_lines), inline=False)
                    pages.append(embed)
            if not pages:
                return await ctx.send(f"\U0001F525 No flagged users currently in {guild.name}.")
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F525 You must specify a member, a voice channel, or use "all".')

    @commands.command(name='mods', hidden=True, help='Lists moderator statistics.')
    @is_owner_developer_coordinator_moderator()
    async def list_moderators(self, ctx, target: Optional[str] = commands.parameter(default=None, description='Voice channel name/mention/ID, "all", or member mention/ID.')) -> None:
        member, channel = await self.get_channel_and_member(ctx, target)
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel)
        if target and target.lower() == 'all':
            if not is_owner_or_dev:
                return await self.handler.send_message(ctx, content='\U0001F525 You are not authorized to list all moderators.')
            query = '''
                SELECT unnest(moderator_channel_ids) AS channel_id, user_id
                FROM users
                WHERE moderator_channel_ids IS NOT NULL
            '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    rows = await conn.fetch(query)
                if not rows:
                    return await self.handler.send_message(ctx, content='\U0001F525 No moderators found in any voice channels.')
                channel_map = defaultdict(list)
                for row in rows:
                    channel_map[row['channel_id']].append(row['user_id'])
                pages = []
                for ch_id, user_ids in sorted(channel_map.items()):
                    vc = ctx.guild.get_channel(ch_id)
                    vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
                    embed = discord.Embed(title=f'🛡️ Moderators for {vc_name}', color=discord.Color.magenta())
                    for uid in user_ids:
                        m = ctx.guild.get_member(uid)
                        name = m.display_name if m else f'User ID {uid}'
                        embed.add_field(name='\u200b', value=f'• {name} (<@{uid}>)', inline=False)
                    pages.append(embed)
                if len(pages) == 1:
                    return await ctx.send(embed=pages[0], allowed_mentions=discord.AllowedMentions.none())
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            except Exception as e:
                await self.handler.send_message(ctx, content=f'Database error: {e}')
                raise
        if member:
            query = '''
                SELECT moderator_channel_ids
                FROM users
                WHERE user_id = $1
            '''
            try:
                async with self.bot.db_pool.acquire() as conn:
                    row = await conn.fetchrow(query, member.id)
                if not row or not row['moderator_channel_ids']:
                    return await ctx.send(f'\U0001F525 {member.display_name} is not a moderator in any channels.')
                embeds = []
                for ch_id in row['moderator_channel_ids']:
                    vc = ctx.guild.get_channel(ch_id)
                    vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
                    embed = discord.Embed(
                        title=f'🛡️ {member.display_name} moderates:',
                        description=f'• {vc_name}',
                        color=discord.Color.magenta()
                    )
                    embeds.append(embed)
                if len(embeds) == 1:
                    return await ctx.send(embed=embeds[0], allowed_mentions=discord.AllowedMentions.none())
                paginator = Paginator(self.bot, ctx, embeds)
                return await paginator.start()
            except Exception as e:
                await self.handler.send_message(ctx, content=f'\U0001F525 Database error: {e}')
                raise
        query = '''
            SELECT user_id
            FROM users
            WHERE $1 = ANY (moderator_channel_ids)
        '''
        try:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(query, channel.id)
            if not rows:
                return await self.handler.send_message(ctx, content=f'\U0001F525 No moderators found for {channel.mention}.')
            lines = []
            for row in rows:
                uid = row['user_id']
                m = ctx.guild.get_member(uid)
                if not m:
                    continue
                lines.append(f'• {m.display_name} — <@{uid}>')
            if not lines:
                return await ctx.send(f"\U0001F525 No moderators currently in {ctx.guild.name}.")
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'🛡️ Moderators for {channel.name}',
                    description='\n'.join(chunk),
                    color=discord.Color.magenta()
                )
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        except Exception as e:
            await self.handler.send_message(ctx, content=f'\U0001F525 Database error: {e}')
            raise
        return await self.handler.send_message(ctx, content='\U0001F525 You must specify a member, a voice channel, or use "all".')

    @commands.command(
        name='mutes',
        hidden=True,
        help='Lists mute statistics.'
    )
    @is_owner_developer_coordinator()
    async def list_mutes(
            self,
            ctx,
            target: Optional[str] = commands.parameter(default=None, description='Voice channel, "all", or user mention/ID.')
    ) -> None:
        member, channel = await self.get_channel_and_member(ctx, target)
        if member:
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''
                                           SELECT am.channel_id,
                                                  am.expires_at,
                                                  COALESCE(mr.reason, 'No reason provided') as reason
                                           FROM active_mutes am
                                                    LEFT JOIN mute_reasons mr
                                                              ON am.user_id = mr.user_id
                                                                  AND am.channel_id = mr.channel_id
                                                                  AND mr.guild_id = $2
                                           WHERE am.user_id = $1
                                           ''', member.id, ctx.guild.id)
                guild_records = []
                for record in records:
                    channel_obj = ctx.guild.get_channel(record['channel_id'])
                    if channel_obj:
                        guild_records.append(record)
                records = guild_records
                if not records:
                    return await ctx.send(f'\U0001F525 {member.mention} is not muted in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
                description_lines = []
                for record in records:
                    channel_obj = ctx.guild.get_channel(record['channel_id'])
                    channel_mention = channel_obj.mention if channel_obj else f'`{record["channel_id"]}`'
                    reason = record['reason']
                    duration_str = "Permanent" if record['expires_at'] is None else discord.utils.format_dt(
                        record['expires_at'], style='R'
                    )
                    description_lines.append(f'• {channel_mention} — {reason} — {duration_str}')
                embed = discord.Embed(
                    title=f'Mute Records for {member.mention}',
                    description='\n'.join(description_lines),
                    color=discord.Color.orange()
                )
                return await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.none())
        if channel:
            is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel)
            if target:
                if target.lower() == 'all':
                    if is_owner_or_dev:
                        async with self.bot.db_pool.acquire() as conn:
                            records = await conn.fetch('''
                                                       SELECT am.user_id,
                                                              am.channel_id,
                                                              am.expires_at,
                                                              am.source,
                                                              COALESCE(mr.reason, 'No reason provided') as reason
                                                       FROM active_mutes am
                                                                LEFT JOIN mute_reasons mr
                                                                          ON am.user_id = mr.user_id
                                                                              AND am.channel_id = mr.channel_id
                                                                              AND mr.guild_id = $1
                                                       ORDER BY am.channel_id, am.user_id
                                                       ''', ctx.guild.id)
                        if not records:
                            return await self.handler.send_message(ctx, content='\U0001F525 No users are currently muted in the server.')
                        grouped = defaultdict(list)
                        for record in records:
                            grouped[record['channel_id']].append(record)
                        pages = []
                        for channel_id, user_entries in sorted(grouped.items()):
                            channel = ctx.guild.get_channel(channel_id)
                            channel_name = channel.mention if channel else f'Unknown Channel ({channel_id})'
                            embed = discord.Embed(
                                title=f'🔇 Active Mutes in {channel_name}',
                                color=discord.Color.orange()
                            )
                            for record in user_entries:
                                user_id = record['user_id']
                                member = ctx.guild.get_member(user_id)
                                name = member.display_name if member else f'User ID {user_id}'
                                mention = member.mention if member else f'`{user_id}`'
                                duration_str = "Permanent" if record['expires_at'] is None else discord.utils.format_dt(
                                    record['expires_at'], style='R'
                                )
                                embed.add_field(
                                    name=name,
                                    value=f'{mention}\nReason: {record["reason"]}\nDuration: {duration_str}',
                                    inline=False
                                )
                            pages.append(embed)
                        paginator = Paginator(self.bot, ctx, pages)
                        return await paginator.start()
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''
                                           SELECT am.user_id,
                                                  am.expires_at,
                                                  COALESCE(mr.reason, 'No reason provided') as reason,
                                                  am.source
                                           FROM active_mutes am
                                                    LEFT JOIN mute_reasons mr
                                                              ON am.user_id = mr.user_id
                                                                  AND am.channel_id = mr.channel_id
                                                                  AND mr.guild_id = $2
                                           WHERE am.channel_id = $1
                                           ''', channel.id, ctx.guild.id)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F525 No users are currently muted in {channel.mention}.')
                description_lines = []
                for record in records:
                    uid = record['user_id']
                    member = ctx.guild.get_member(uid)
                    if not member:
                        continue
                    name = member.display_name
                    duration_str = "Permanent" if record['expires_at'] is None else discord.utils.format_dt(
                        record['expires_at'], style='R'
                    )
                    description_lines.append(f'• {name} — <@{uid}> — {duration_str}')
                if not description_lines:
                    return await ctx.send(f"\U0001F525 No muted users currently in {ctx.guild.name}.")
                chunk_size = 18
                pages = []
                for i in range(0, len(description_lines), chunk_size):
                    chunk = description_lines[i:i + chunk_size]
                    embed = discord.Embed(
                        title=f'\U0001F507 Muted Users in {channel.mention}',
                        color=discord.Color.orange()
                    )
                    embed.add_field(name="Muted Users", value='\n'.join(chunk), inline=False)
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
        else:
            return await self.handler.send_message(ctx, content='\U0001F525 You must specify a member, a voice channel or be connected to a voice channel.')
            
    @commands.command(
        name='tmutes',
        hidden=True,
        help='Lists text-mute statistics.'
    )
    @is_owner_developer_coordinator()
    async def list_text_mutes(
            self,
            ctx,
            target: Optional[str] = commands.parameter(default=None, description='Optional: "all", channel name/ID/mention, or user mention/ID.')
    ) -> None:
        member, channel = await self.get_channel_and_member(ctx, target)
        if not member and not channel:
            return await self.handler.send_message(ctx, content='\U0001F525 Could not resolve a valid text channel or user.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F525 You do not have permissions to use this command in {channel.mention}')
        if target:
            if is_owner_or_dev:
                if target.lower() == 'all':
                    async with self.bot.db_pool.acquire() as conn:
                        records = await conn.fetch('''
                                                   SELECT user_id, channel_id, reason, source
                                                   FROM text_mutes
                                                   WHERE guild_id = $1
                                                   ''', ctx.guild.id)
                    if not records:
                        return await self.handler.send_message(ctx, content='\U0001F525 No users are currently text-muted in this server.')
                    grouped = defaultdict(list)
                    for record in records:
                        grouped[record['channel_id']].append(record)
                    pages = []
                    for channel_id, entries in sorted(grouped.items()):
                        channel = ctx.guild.get_channel(channel_id)
                        ch_name = channel.mention if channel else f'Unknown Channel ({channel_id})'
                        embed = discord.Embed(
                            title=f'🔇 Text Mutes in {ch_name}',
                            color=discord.Color.orange()
                        )
                        for entry in entries:
                            user = ctx.guild.get_member(entry['user_id'])
                            mention = user.mention if user else f'`{entry["user_id"]}`'
                            reason = entry['reason']
                            source = entry['source']
                            embed.add_field(
                                name=mention,
                                value=f'Reason: {reason}\nSource: `{source}`',
                                inline=False
                            )
                        pages.append(embed)
                    paginator = Paginator(self.bot, ctx, pages)
                    return await paginator.start()
            else:
                return await self.handler.send_message(ctx, content='\U0001F525 You are not authorized to specify a target. This command defaults to the current channel.')
        async with self.bot.db_pool.acquire() as conn:
            if channel and not member:
                records = await conn.fetch('''
                                           SELECT user_id, reason, source
                                           FROM text_mutes
                                           WHERE channel_id = $1
                                             AND guild_id = $2
                                           ''', channel.id, ctx.guild.id)
                if not records:
                    return await ctx.send(f'\U0001F525 No users are currently text-muted in {channel.mention}.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for record in records:
                    uid = record['user_id']
                    mention = f"<@{uid}>"
                    lines.append(f'• {mention} — {record["reason"]} (via `{record["source"]}`)')
                chunk_size = 18
                pages = []
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i + chunk_size]
                    embed = discord.Embed(
                        title=f'Text-Muted Users in {channel.mention}',
                        description='\n'.join(chunk),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            elif member:
                records = await conn.fetch('''
                                           SELECT channel_id, reason, source
                                           FROM text_mutes
                                           WHERE user_id = $1
                                             AND guild_id = $2
                                           ''', member.id, ctx.guild.id)
                if not records:
                    return await ctx.send(f'\U0001F525 {member.mention} is not text-muted in any channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for record in records:
                    channel_id = record['channel_id']
                    ch = ctx.guild.get_channel(channel_id)
                    if not ch:
                        continue
                    ch_mention = f'<#{channel_id}>'
                    lines.append(f'• {ch_mention} — {record["reason"]} (via `{record["source"]}`)')
                if not lines:
                    return await ctx.send(f"\U0001F525 No text mute records found for {member.mention}.")
                chunk_size = 18
                pages = []
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i + chunk_size]
                    embed = discord.Embed(
                        title=f'Text Mute Records for {member.mention}',
                        description='\n'.join(chunk),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()

    @commands.command(
        name='ls',
        help='List users cowed as going vegan in this guild.'
    )
    @is_owner_developer_coordinator_moderator()
    async def list_members(
        self,
        ctx,
        target: Optional[str] = commands.parameter(default=None, description='Channel ID.')
    ) -> None:
        member, channel = await self.get_channel_and_member(ctx, target)
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(ctx, channel)
        if target is not None and not (is_owner_or_dev or is_mod_or_coord):
            return await self.handler.send_message(ctx, content='\U0001F525 You are not authorized to specify a channel. Use this command while connected to a channel.')
        if not is_owner_or_dev and not is_mod_or_coord:
            return await self.handler.send_message(ctx, content=f'\U0001F525 You do not have permissions to use this command in {channel.mention}')
        guild = ctx.guild
        try:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                                        SELECT user_id
                                        FROM users
                                        WHERE going_vegan_channel_ids IS NOT NULL
                                          AND $1 = ANY (going_vegan_channel_ids)
                                        ''', channel.id)
            if not rows:
                return await self.handler.send_message(ctx, content='\U0001F525 No users are flagged as going vegan in this channel.')
            lines = []
            for row in rows:
                uid = row['user_id']
                member = guild.get_member(uid)
                if not member:
                    continue
                name = member.display_name
                lines.append(f'• {name} — <@{uid}>')
            if not lines:
                return await ctx.send(f"\U0001F525 No new vegans currently in {guild.name}.")
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'🐮 New Vegan in {guild.name}',
                    description='\n'.join(chunk),
                    color=discord.Color.green()
                )
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        except Exception as e:
            await self.handler.send_message(ctx, content=f'Database error: {e}')
            raise
            
    @commands.command(name='smute', help='Mutes a member throughout the entire guild.')
    @commands.check(lambda ctx: ctx.bot.get_cog("Hybrid").can_server_mute(ctx))
    async def server_mute(
            self,
            ctx,
            member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Optionally include a reason for the mute.')
    ) -> None:
        try:
            await is_owner_block(ctx, member)
        except commands.CheckFailure:
            return await self.handler.send_message(ctx, content='\U0001F525 You are not allowed to mute the owner.')
        member, _ = await self.get_channel_and_member(ctx, member)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                               INSERT INTO users (user_id, server_mute_guild_ids)
                               VALUES ($1, ARRAY[$2]::BIGINT[]) ON CONFLICT (user_id) DO
                               UPDATE
                                   SET server_mute_guild_ids = (
                                   SELECT ARRAY(
                                   SELECT DISTINCT unnest(COALESCE (u.server_mute_guild_ids, '{}') || ARRAY[$2])
                                   )
                                   FROM users u WHERE u.user_id = EXCLUDED.user_id
                                   ),
                                   updated_at = NOW()
                               ''', member.id, ctx.guild.id)
            await conn.execute('''
                               INSERT INTO server_mute_reasons (guild_id, user_id, reason)
                               VALUES ($1, $2, $3) ON CONFLICT (guild_id, user_id)
                               DO
                               UPDATE SET reason = EXCLUDED.reason
                               ''', ctx.guild.id, member.id, reason)
        if member.voice and member.voice.channel:
            await member.edit(mute=True)
            return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been server muted in <#{member.voice.channel.id}> for reason: {reason}', allowed_mentions=discord.AllowedMentions.none())
        else:
            return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been server muted. They are not currently in a voice channel.', allowed_mentions=discord.AllowedMentions.none())

    @commands.hybrid_command(name="rmute", description="Full room mute")
    @is_owner()
    async def room_mute(
        self,
        ctx,
        channel: discord.VoiceChannel,
        duration_hours: Optional[str] = commands.parameter(default="8", description="Duration of mute in hours. Example 0 (permanent), 30m, 2h, 5d."),
        *,
        reason: str = commands.parameter(default="", description="Optional reason (required for permanent mutes).")
    ) -> None:
        expires_at, duration_display = self.parse_duration(duration_hours)
        if (expires_at == "0" or expires_at is None) and (not is_owner_developer_coordinator("mute") or not reason.strip()):
            return await self.handler.send_message(ctx, content="\U0001F525 Reason required and coordinator-only for permanent mutes.")
        bot_owner_id = int(os.environ.get("DISCORD_OWNER_ID", "0"))
        author_id = ctx.author.id
        mute_source = "bot_owner" if author_id == bot_owner_id else "bot"
        muted_members = []
        failed_members = []
        async with self.bot.db_pool.acquire() as conn:
            for member in channel.members:
                if member.id == author_id:
                    continue
                try:
                    await conn.execute(
                        """
                        INSERT INTO active_mutes (user_id, channel_id, source, issuer_id, expires_at)
                        VALUES ($1, $2, $3, $4, $5)
                        ON CONFLICT (user_id, channel_id) DO UPDATE
                        SET source = EXCLUDED.source, issuer_id = EXCLUDED.issuer_id, expires_at = EXCLUDED.expires_at
                        """,
                        member.id,
                        channel.id,
                        mute_source,
                        author_id,
                        expires_at
                    )
                    await conn.execute(
                        """
                        INSERT INTO users (user_id, mute_channel_ids)
                        VALUES ($1, ARRAY[$2]::BIGINT[])
                        ON CONFLICT (user_id) DO UPDATE
                        SET mute_channel_ids = (
                            SELECT ARRAY(
                                SELECT DISTINCT unnest(COALESCE(u.mute_channel_ids, '{}') || ARRAY[$2])
                            )
                            FROM users u WHERE u.user_id = EXCLUDED.user_id
                        ),
                        updated_at = NOW()
                        """,
                        member.id,
                        channel.id
                    )
                    await conn.execute(
                        """
                        INSERT INTO mute_reasons (guild_id, user_id, reason, channel_id)
                        VALUES ($1, $2, $3, $4)
                        ON CONFLICT (guild_id, user_id, channel_id)
                        DO UPDATE SET reason = EXCLUDED.reason
                        """,
                        ctx.guild.id,
                        member.id,
                        reason or "No reason provided",
                        channel.id
                    )
                    await conn.execute(
                        """
                        INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id,
                                                     guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                        """,
                        "voice_mute",
                        member.id,
                        ctx.author.id,
                        ctx.guild.id,
                        channel.id,
                        reason or "No reason provided"
                    )
                    if member.voice and member.voice.channel and member.voice.channel.id == channel.id:
                        await member.edit(mute=True)
                    muted_members.append(member)
                except Exception as e:
                    logger.warning(f"Failed to mute {member}: {e}")
                    failed_members.append(member)
        summary = f"{self.get_random_emoji()} Muted {len(muted_members)} member(s) in {channel.mention} {duration_display}.\nReason: {reason or 'No reason provided'}"
        if failed_members:
            summary += f"\n⚠️ Failed to mute {len(failed_members)} member(s)."
        return await ctx.send(summary, allowed_mentions=discord.AllowedMentions.none())
    
    
    @commands.command(name="rmv")
    @is_owner()
    async def room_move_all(self, ctx: commands.Context, source_id: int, target_id: int):
        source_channel = ctx.guild.get_channel(source_id)
        target_channel = ctx.guild.get_channel(target_id)
        if not source_channel or not target_channel:
            await ctx.send("🔥 One or both channel IDs are invalid.")
            return
        await self.move_all_members(source_channel, target_channel)
        await ctx.send(f"{self.get_random_emoji()} Moved all members from `{source_channel.name}` to `{target_channel.name}`.")
        
    @commands.command(name="xrmute", help="Unmutes all members in a specific VC (except yourself).")
    @is_owner()
    async def room_unmute_all(
        self,
        ctx,
        channel: discord.VoiceChannel,
        *,
        reason: str = commands.parameter(default="N/A", description="Include a reason for the unmute.")
    ) -> None:
        unmuted_members = []
        skipped_members = []
        failed_members = []
        async with self.bot.db_pool.acquire() as conn:
            for member in channel.members:
                if member.id == ctx.author.id:
                    skipped_members.append(member)
                    continue
                try:
                    row = await conn.fetchrow(
                        """
                        SELECT source, issuer_id FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2
                        """,
                        member.id,
                        channel.id
                    )
                    if not row:
                        skipped_members.append(member)
                        continue
                    if row["source"] not in ("bot", "manual", "bot_owner"):
                        skipped_members.append(member)
                        continue
                    if member.voice and member.voice.channel and member.voice.channel.id == channel.id:
                        await member.edit(mute=False)
                    await conn.execute(
                        """
                        DELETE FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2 AND source IN ('bot', 'manual', 'bot_owner')
                        """,
                        member.id,
                        channel.id
                    )
                    if row["source"] in ("bot", "bot_owner"):
                        await conn.execute(
                            """
                            UPDATE users
                            SET mute_channel_ids = array_remove(mute_channel_ids, $2),
                                updated_at = NOW()
                            WHERE user_id = $1
                            """,
                            member.id,
                            channel.id
                        )
                    elif row["source"] == "manual":
                        await conn.execute(
                            """
                            UPDATE users
                            SET manual_mute_channels = array_remove(manual_mute_channels, $2),
                                updated_at = NOW()
                            WHERE user_id = $1
                            """,
                            member.id,
                            channel.id
                        )
                    await conn.execute(
                        """
                        INSERT INTO mute_reasons (guild_id, user_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4)
                        ON CONFLICT (guild_id, user_id, channel_id)
                        DO UPDATE SET reason = EXCLUDED.reason
                        """,
                        ctx.guild.id,
                        member.id,
                        channel.id,
                        f"Unmuted: {reason}"
                    )
                    await conn.execute(
                        """
                        INSERT INTO moderation_logs (action_type, target_user_id, executor_user_id,
                                                     guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                        """,
                        "unmute",
                        member.id,
                        ctx.author.id,
                        ctx.guild.id,
                        channel.id,
                        f"Unmuted via unmute_all: {reason}"
                    )
                    unmuted_members.append(member)
                except Exception as e:
                    logger.warning(f"Unmute failed for {member}: {e}")
                    failed_members.append(member)
        summary = f"{self.get_random_emoji()} Unmuted {len(unmuted_members)} member(s) in {channel.mention}."
        if skipped_members:
            summary += f"\n⏭️ Skipped {len(skipped_members)} (not muted or not bot-muted)."
        if failed_members:
            summary += f"\n⚠️ Failed to unmute {len(failed_members)}."
        return await ctx.send(summary, allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='xadmin', help='Revokes server mute privileges from a user.')
    @is_owner()
    async def revoke_server_muter(self, ctx, member: str):
        member, _ = await self.get_channel_and_member(ctx, member)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                               UPDATE users
                               SET server_muter_guild_ids = (SELECT ARRAY(
                                                                        SELECT unnest(server_muter_guild_ids)
                        EXCEPT SELECT $2
                                                                    ))
                               WHERE user_id = $1
                               ''', member.id, ctx.guild.id)
        self.server_muters.get(ctx.guild.id, set()).discard(member.id)
        return await ctx.send(f'{self.get_random_emoji()} {member.mention} no longer has server mute privileges.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.command(name='xsmute', help='Unmutes a member throughout the entire guild.')
    @commands.check(lambda ctx: ctx.bot.get_cog("Hybrid").can_server_mute(ctx))
    async def unsmute(
            self,
            ctx,
            member: str = commands.parameter(description='Tag a user or include their snowflake ID.')
    ) -> None:
        member, _ = await self.get_channel_and_member(ctx, member)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                               UPDATE users
                               SET server_mute_guild_ids = (SELECT ARRAY(
                                                                       SELECT unnest(server_mute_guild_ids)
                        EXCEPT SELECT $2
                                                                   ))
                               WHERE user_id = $1
                               ''', member.id, ctx.guild.id)
            await conn.execute('''
                               DELETE
                               FROM server_mute_reasons
                               WHERE user_id = $1
                                 AND guild_id = $2
                               ''', member.id, ctx.guild.id)
        if member.voice and member.voice.channel:
            await member.edit(mute=False)
        return await ctx.send(f'{self.get_random_emoji()} {member.mention} has been server unmuted.', allowed_mentions=discord.AllowedMentions.none())

    async def cog_load(self) -> None:
        if not hasattr(self, "_loaded_aliases"):
            self._loaded_aliases = set()
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT guild_id, alias_type, alias_name, channel_id FROM command_aliases'
            )
            for row in rows:
                guild_id = row['guild_id']
                alias_type = row['alias_type']
                alias_name = row['alias_name']
                channel_id = row['channel_id']
                if channel_id is None:
                    continue
                self.bot.command_aliases.setdefault(guild_id, {}).setdefault(alias_type, {})[alias_name] = channel_id
                if alias_name in self._loaded_aliases:
                    continue
                cmd = None
                if alias_type == 'mute':
                    cmd = self.create_voice_mute_alias(alias_name)
                elif alias_type == 'unmute':
                    cmd = self.create_unmute_alias(alias_name)
                elif alias_type == 'ban':
                    cmd = self.create_ban_alias(alias_name)
                elif alias_type == 'unban':
                    cmd = self.create_unban_alias(alias_name)
                elif alias_type == 'cow':
                    cmd = self.create_cow_alias(alias_name)
                elif alias_type == 'uncow':
                    cmd = self.create_uncow_alias(alias_name)
                elif alias_type == 'flag':
                    cmd = self.create_flag_alias(alias_name)
                elif alias_type == 'unflag':
                    cmd = self.create_unflag_alias(alias_name)
                elif alias_type == 'tmute':
                    cmd = self.create_text_mute_alias(alias_name)
                elif alias_type == 'untmute':
                    cmd = self.create_untextmute_alias(alias_name)
                if cmd and not self.bot.get_command(alias_name):
                    self.bot.add_command(cmd)
                    self._loaded_aliases.add(alias_name)


    async def get_channel_and_member(self, ctx: commands.Context, value: Optional[Union[str, discord.TextChannel, discord.VoiceChannel, discord.Member]]) -> tuple[Optional[discord.Member], Optional[Union[discord.TextChannel, discord.VoiceChannel]]]:
        member = None
        channel = None
        try:
            if isinstance(value, (discord.TextChannel, discord.VoiceChannel)):
                channel = value
            elif isinstance(value, discord.Member):
                member = value
            elif isinstance(value, str):
                if value.isdigit():
                    entity_id = int(value)
                    potential_channel = ctx.guild.get_channel(entity_id)
                    if isinstance(potential_channel, (discord.TextChannel, discord.VoiceChannel)):
                        channel = potential_channel
                    else:
                        member = ctx.guild.get_member(entity_id)
                        if member is None:
                            try:
                                member = await ctx.guild.fetch_member(entity_id)
                            except discord.NotFound:
                                member = None
                elif value.startswith('<#') and value.endswith('>'):
                    channel_id = int(value[2:-1])
                    channel = ctx.guild.get_channel(channel_id)
                elif value.startswith('<@') and value.endswith('>'):
                    member_id = int(value[2:-1].replace('!', ''))
                    member = ctx.guild.get_member(member_id)
                    if member is None:
                        try:
                            member = await ctx.guild.fetch_member(member_id)
                        except discord.NotFound:
                            member = None
                if channel is None:
                    channel = discord.utils.find(lambda c: c.name.lower() == value.lower(), ctx.guild.text_channels)
                if member is None:
                    member = discord.utils.find(
                        lambda m: m.name.lower() == value.lower() or (m.nick and m.nick.lower() == value.lower()),
                        ctx.guild.members)
            if channel is None:
                channel = ctx.channel
            if member is None and channel is None:
                await self.handler.send_message(ctx, content='\U0001F525 Could not resolve a valid guild channel or member from your input.')
        except (ValueError, AttributeError) as e:
            logger.warning(e)
            return None, None
        return member, channel
        
    async def load_server_muters(self):
        await self.bot.wait_until_ready()
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT user_id, server_muter_guild_ids FROM users WHERE server_muter_guild_ids IS NOT NULL')
            for row in rows:
                user_id = row['user_id']
                for guild_id in row['server_muter_guild_ids']:
                    self.server_muters[guild_id].add(user_id)
                    
    async def move_all_members(self, source_channel: discord.VoiceChannel, target_channel: discord.VoiceChannel):
        if not isinstance(source_channel, discord.VoiceChannel) or not isinstance(target_channel, discord.VoiceChannel):
            raise ValueError("🔥 Both source and target must be voice channels.")
        for member in source_channel.members:
            try:
                await member.move_to(target_channel)
            except discord.Forbidden:
                print(f"🔥 Missing permissions to move {member}.")
            except discord.HTTPException:
                print(f"⚠️ Failed to move {member} due to a network error.")
                
    def parse_duration(self, duration: Optional[str]) -> tuple[Optional[datetime], str]:
        if duration is None:
            delta = timedelta(hours=24)
            return datetime.now(timezone.utc) + delta, "for 24 hour(s)"
        duration = duration.strip().lower()
        if duration in ("0", "0h", "0d", "0m"):
            return None, "permanently"
        if duration.endswith("d"):
            days = int(duration[:-1])
            delta = timedelta(days=days)
            return datetime.now(timezone.utc) + delta, f"for {days} day(s)"
        if duration.endswith("h"):
            hours = int(duration[:-1])
            delta = timedelta(hours=hours)
            return datetime.now(timezone.utc) + delta, f"for {hours} hour(s)"
        if duration.endswith("m"):
            minutes = int(duration[:-1])
            delta = timedelta(minutes=minutes)
            return datetime.now(timezone.utc) + delta, f"for {minutes} minute(s)"
        hours = int(duration)
        delta = timedelta(hours=hours)
        return datetime.now(timezone.utc) + delta, f"for {hours} hour(s)"

    def perform_backup(self, db_user: str, db_name: str, db_host: str, db_password: str, backup_dir: str) -> str:
        timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        backup_file = os.path.join(backup_dir, f'backup_{timestamp}.sql')
        dump_command = [
            'pg_dump',
            '-U', db_user,
            '-h', db_host,
            '-d', db_name,
            '-F', 'p',
            '-f', backup_file,
        ]
        env = os.environ.copy()
        env["PGPASSWORD"] = db_password
        result = subprocess.run(
            dump_command,
            capture_output=True,
            text=True,
            env=env,
        )
        if result.returncode != 0:
            raise RuntimeError(f'Backup failed: {result.stderr}')
        return backup_file

    def setup_backup_directory(self, backup_dir: str) -> str:
        os.makedirs(backup_dir, exist_ok=True)
        return backup_dir

    async def vegan_db(self, ctx: commands.Context) -> bool:
        channel = ctx.channel
        user_id = ctx.author.id
        channel_id = channel.id
        guild_id = ctx.guild.id
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow("""
                SELECT coordinator_ids, coordinator_channel_ids
                FROM users
                WHERE user_id = $1
            """, user_id)
        if not row:
            await ctx.send("🔥 You are not registered in the database.")
            return
        coordinator_ids = row["coordinator_ids"] or []
        coordinator_channel_ids = row["coordinator_channel_ids"] or []
        is_coordinator = (user_id in coordinator_ids) and (channel_id in coordinator_channel_ids)
        is_vegan_channel = "vegan" in channel.name.lower()
        if is_coordinator and is_vegan_channel:
            return True
        else:
            return False

    @staticmethod
    def can_server_mute(ctx):
        cog = ctx.bot.get_cog('Hybrid')
        return cog and ctx.author.id in cog.server_muters.get(ctx.guild.id, set())
        
    def coord_mod_check_with_channel(channel_arg: str, alias_type: Optional[str] = None):
        def decorator(func):
            async def wrapper(ctx, *args, **kwargs):
                channel = kwargs.get(channel_arg) if channel_arg in kwargs else None
                if channel is None and args:
                    channel = args[0] if channel_arg == "channel" else None
                if channel:
                    ctx._target_channel_id = channel.id if hasattr(channel, "id") else int(channel)
                return await func(ctx, *args, **kwargs)
            return wrapper
        return lambda f: is_owner_developer_coordinator_moderator(alias_type)(wrapper)

async def setup(bot: DiscordBot):
    cog = Hybrid(bot)
    await bot.add_cog(cog)
