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
import datetime
import discord
import inspect
import os
import re
from typing import Any, Coroutine, Optional

from discord.ext.commands import Command

from vyrtuous.inc.helpers import *
from vyrtuous.service.check_service import *
from vyrtuous.service.discord_message_service import DiscordMessageService, Paginator
from vyrtuous.utils.setup_logging import logger

PERMISSION_ORDER = ['Owner', 'Developer', 'Coordinator', 'Moderator', 'Everyone']
class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
  
    @commands.hybrid_command(
        name='alias',
        help='Set an alias for a mute, unmute, ban, unban or flag action.'
    )
    @commands.check(is_owner_developer_coordinator)
    async def create_alias(
        self,
        ctx,
        alias_type: str = commands.parameter(description='One of: `mute`, `unmute`, `ban`, `unban`, `flag`, `unflag`'),
        alias_name: str = commands.parameter(description='Alias/Pseudonym'),
        channel_id: str = commands.parameter(description='Voice channel')
    ) -> None:
        alias_type = alias_type.lower()
        guild_id = ctx.guild.id
        valid_types = {'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag'}
        if alias_type not in valid_types:
            await self.handler.send_message(ctx, f'❌ Invalid alias type. Must be one of: {", ".join(valid_types)}', ephemeral=True)
            return
        if not alias_name.strip():
            await self.handler.send_message(ctx, '❌ Alias name cannot be empty.', ephemeral=True)
            return
        def resolve_channel(value: str):
            if value.isdigit():
                return ctx.guild.get_channel(int(value))
            if value.startswith('<#') and value.endswith('>'):
                return ctx.guild.get_channel(int(value.strip('<#>')))
            return discord.utils.get(ctx.guild.voice_channels, name=value)
        channel = resolve_channel(channel_id)
        if not channel or channel.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, '❌ Could not resolve a valid voice channel.', ephemeral=True)
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
                return await self.handler.send_message(
                    ctx,
                    f'❌ Alias `{alias_name}` ({alias_type}) already exists and is set to {channel_mention}.',
                    ephemeral=True
                )
            if self.bot.get_command(alias_name):
                return await self.handler.send_message(
                    ctx,
                    f'❌ A command named `{alias_name}` already exists.',
                    ephemeral=True
                )
            await conn.execute(
                '''
                INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id)
                VALUES ($1, $2, $3, $4)
                ''',
                guild_id, alias_type, alias_name, channel.id
            )
        self.bot.command_aliases.setdefault(guild_id, {}).setdefault(alias_type, {})[alias_name] = channel.id
        if alias_type == 'mute':
            cmd = self.create_voice_mute_alias(alias_name)
        elif alias_type == 'unmute':
            cmd = self.create_unmute_alias(alias_name)
        elif alias_type == 'ban':
            cmd = self.create_ban_alias(alias_name)
        elif alias_type == 'unban':
            cmd = self.create_unban_alias(alias_name)
        elif alias_type == 'flag':
            cmd = self.create_flag_alias(alias_name)
        elif alias_type == 'unflag':
            cmd = self.create_unflag_alias(alias_name)
        self.bot.add_command(cmd)
        await self.handler.send_message(
            ctx,
            content=f'✅ Alias `{alias_name}` ({alias_type}) set to {channel.mention}.'
        )

    def create_ban_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(
            name=command_name,
            help='Ban a user from a voice channel.'
        )
        @commands.check(is_owner_developer_coordinator_moderator)
        async def ban_alias(
                ctx,
                member: str = commands.parameter(description='Mention or user ID of the member to ban.'),
                duration_hours: Optional[int] = commands.parameter(default=24, description='Duration of ban in hours (0 = permanent).'),
                *,
                reason: str = commands.parameter(default='', description='Reason for ban (required for permanent).')
        ) -> Coroutine[Any, Any, None]:
            try:
                await is_owner_block(ctx, member)
            except commands.CheckFailure as e:
                return await self.handler.send_message(ctx, content='❌ You are not allowed to ban the owner.')
            command_name = ctx.invoked_with
            guild_id = ctx.guild.id
            member_id = None
            if member.isdigit():
                member_id = int(member)
            elif member.startswith('<@') and member.endswith('>'):
                try:
                    member_id = int(member.strip('<@!>'))
                except ValueError:
                    pass
            member_object = ctx.guild.get_member(member_id) if member_id else None
            if not member_object:
                return await self.handler.send_message(ctx, content='❌ Could not resolve a valid member.')
            if duration_hours == 0 and (not await is_owner_developer_coordinator(ctx) or not reason.strip()):
                return await self.handler.send_message(
                    ctx,
                    content='❌ Reason required and coordinator-only for permanent bans.'
                )
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('ban', {}).get(command_name)
            if not static_channel_id:
                return await self.handler.send_message(
                    ctx,
                    content=f'❌ No channel alias mapping found for `{command_name}`.'
                )
            channel = ctx.guild.get_channel(static_channel_id)
            if not channel or not isinstance(channel, discord.VoiceChannel):
                return await self.handler.send_message(
                    ctx,
                    content=f'❌ Could not resolve a valid voice channel for ID `{static_channel_id}`.'
                )
            try:
                await channel.set_permissions(
                    member_object,
                    view_channel=False,
                    reason=f"Banned from <#{channel.id}>: {reason or 'No reason provided'}"
                )
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content='❌ Missing permissions to deny channel access.')
            if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == channel.id:
                try:
                    await member_object.move_to(None, reason="Banned from this channel")
                except discord.Forbidden:
                    await ctx.send(f"⚠️ Could not disconnect <@{member_object.id}> from <#{channel.id}>.", allowed_mentions=discord.AllowedMentions.none())
                except Exception as e:
                    logger.exception(f"⚠️ Unexpected error while disconnecting user: {e}")
            try:
                async with self.bot.db_pool.acquire() as conn:
                    await conn.execute(
                        '''
                        INSERT INTO active_bans (user_id, channel_id, expires_at)
                        VALUES ($1, $2, $3) ON CONFLICT (user_id, channel_id)
                        DO
                        UPDATE SET expires_at = EXCLUDED.expires_at
                        ''',
                        member_object.id,
                        channel.id,
                        None if duration_hours == 0 else discord.utils.utcnow() + datetime.timedelta(
                            hours=duration_hours)
                    )
                    await conn.execute(
                        '''
                        UPDATE users
                        SET ban_channel_ids =
                                CASE
                                    WHEN NOT $2 = ANY (ban_channel_ids) THEN array_append(ban_channel_ids, $2)
                                    ELSE ban_channel_ids
                                    END,
                            updated_at      = NOW()
                        WHERE user_id = $1
                        ''',
                        member_object.id, channel.id
                    )
                    await conn.execute(
                        '''
                        INSERT INTO ban_reasons (guild_id, user_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4) ON CONFLICT (guild_id, user_id, channel_id)
                        DO
                        UPDATE SET reason = EXCLUDED.reason
                        ''',
                        guild_id, member_object.id, channel.id, reason or 'No reason provided'
                    )
            except Exception as e:
                logger.warning(f"🔥 Database error occurred: {e}")
            await ctx.send(f'🔇 {member_object.mention} has been banned from <#{channel.id}> {"permanently" if duration_hours == 0 else f"for {duration_hours} hour(s)"}', allowed_mentions=discord.AllowedMentions.none())
            if duration_hours > 0:
                expires_at = datetime.datetime.utcnow() + datetime.timedelta(hours=duration_hours)
                try:
                    async with self.bot.db_pool.acquire() as conn:
                        await conn.execute(
                            '''
                            INSERT INTO ban_expirations (user_id, channel_id, expires_at)
                            VALUES ($1, $2, $3) ON CONFLICT (user_id, channel_id)
                            DO
                            UPDATE SET expires_at = EXCLUDED.expires_at
                            ''',
                            member_object.id, channel.id, expires_at
                        )
                except Exception as e:
                    logger.warning(f"🔥 Database error occurred while setting ban expiration: {e}")
            return None
        return ban_alias
        
    @commands.hybrid_command(name='coord', help='Grants coordinator access for a specific voice channel.')
    @commands.check(is_owner_developer)
    async def create_coordinator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel: str = commands.parameter(description='Mention a channel or provide its ID.'),
    ) -> None:
        channel_id = None
        channel_obj = None
        if channel.isdigit():
            channel_id = int(channel)
            channel_obj = ctx.guild.get_channel(channel_id)
        elif channel.startswith('<#') and channel.endswith('>'):
            try:
                channel_id = int(channel[2:-1])
                channel_obj = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            channel_obj = discord.utils.get(ctx.guild.voice_channels, name=channel)
            if channel_obj:
                channel_id = channel_obj.id
        if not channel_obj or not isinstance(channel_obj, discord.VoiceChannel):
            return await self.handler.send_message(
                ctx,
                content='❌ Could not find a valid voice channel from your input.'
            )
        member_id = None
        member_obj = None
        if member.isdigit():
            member_id = int(member)
            member_obj = ctx.guild.get_member(member_id)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member[2:-1].lstrip('!'))
                member_obj = ctx.guild.get_member(member_id)
            except ValueError:
                pass
        if not member_obj:
            return await self.handler.send_message(
                ctx,
                content='❌ Could not find a valid guild member from your input.'
            )
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
            ''', member_obj.id, ctx.guild.id, channel_id)
        await ctx.send(f'✅ {member_obj.mention} has been granted coordinator rights in {channel_obj.mention}.',
            allowed_mentions=discord.AllowedMentions.none()
        )

    
    @commands.hybrid_command(name='dev', help='Elevates a user\'s permissions to a bot developer.')
    @commands.check(is_owner)
    async def create_developer(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ) -> None:
        guild_id = ctx.guild.id
        member_id = None
        member_object = None
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
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
            ''', member_object.id, guild_id)
        await ctx.send(f'{member_object.mention} has been granted developer rights in this guild.',
        allowed_mentions=discord.AllowedMentions.none())

    def create_flag_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(
            name=command_name,
            help='Flag a user in the database for the voice channel mapped to this alias.'
        )
        @commands.check(is_owner_developer_coordinator_moderator)
        async def flag_alias(
                ctx,
                user: str = commands.parameter(description='Tag a user or include their snowflake ID.')
        ) -> None:
            guild_id = ctx.guild.id
            flag_aliases = self.bot.command_aliases.get(guild_id, {}).get('flag', {})
            channel_id = flag_aliases.get(command_name)
            if not channel_id:
                await self.handler.send_message(ctx, content=f'❌ No flag alias configured for `{command_name}`.')
                return
            if not user:
                await self.handler.send_message(ctx, content='❌ You must provide a user ID or mention.')
                return
            if re.fullmatch(r'<@!?\d+>', user):
                user_id = int(re.sub(r'\D', '', user))
            elif user.isdigit():
                user_id = int(user)
            else:
                await self.handler.send_message(ctx, content='❌ Invalid user ID or mention.')
                return
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
                    already_flagged = await conn.fetchval(select_sql, user_id, channel_id)
                    if already_flagged:
                        await self.handler.send_message(ctx,content=f'ℹ️ <@{user_id}> is already flagged for <#{channel_id}>.')
                        return
                    await conn.execute(insert_sql, user_id, channel_id)
                    await ctx.send(f'✅ Flagged <@{user_id}> for channel <#{channel_id}>.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                await self.handler.send_message(ctx, content=f'❌ Database error: {e}')
                raise
        return flag_alias

    @commands.hybrid_command(name='mod', help='Elevates a user\'s permission to VC moderator for a specific channel.')
    @commands.check(is_owner_developer_coordinator)
    async def create_moderator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel: str = commands.parameter(description='Tag a channel or include its snowflake ID.')
    ) -> None:
        member_id = None
        member_object = None
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        resolved_channel = None
        if channel.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel))
        elif channel.startswith('<#') and channel.endswith('>'):
            try:
                channel_id = int(channel.strip('<#>'))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel.lower():
                    resolved_channel = vc
                    break
        if not isinstance(resolved_channel, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content='Could not resolve a valid **voice** channel from your input.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO users (user_id, moderator_ids, moderator_channel_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[], ARRAY[$3]::BIGINT[])
                ON CONFLICT (user_id) DO UPDATE
                SET 
                    moderator_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(
                                COALESCE(users.moderator_ids, ARRAY[]::BIGINT[]) || 
                                ARRAY[$2]::BIGINT[]
                            )
                        )
                    ),
                    moderator_channel_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(
                                COALESCE(users.moderator_channel_ids, ARRAY[]::BIGINT[]) || 
                                ARRAY[$3]::BIGINT[]
                            )
                        )
                    ),
                    updated_at = NOW()
            ''', member_object.id, ctx.guild.id, resolved_channel.id)
        await ctx.send(f'{member_object.mention} has been granted VC moderator access in {resolved_channel.name}.',
            allowed_mentions=discord.AllowedMentions.none()
        )

    def create_voice_mute_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(name=command_name, help='Mutes a member in a specific VC.')
        @commands.check(is_owner_developer_coordinator_moderator)
        async def voice_mute_alias(
            ctx,
            member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Optionally include a reason for the mute.')
        ) -> None:
            try:
                await is_owner_block(ctx, member)
            except commands.CheckFailure as e:
                return await self.handler.send_message(ctx, content='❌ You are not allowed to ban the owner.')
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member.isdigit():
                member_id = int(member)
            elif member.startswith('<@') and member.endswith('>'):
                try:
                    member_id = int(member.strip('<@!>'))
                except ValueError:
                    pass
            if member_id:
                member_object = ctx.guild.get_member(member_id)
            if not member_object:
                return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('mute', {}).get(command_name)
            bot_owner_id = int(os.environ.get("DISCORD_OWNER_ID", "0"))
            server_owner_id = ctx.guild.owner_id
            author_id = ctx.author.id
            if author_id == server_owner_id:
                mute_source = "owner"
            elif author_id == bot_owner_id:
                mute_source = "bot_owner"
            else:
                mute_source = "bot"
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute('''
                    INSERT INTO active_mutes (user_id, channel_id, source, issuer_id)
                    VALUES ($1, $2, $3, $4)
                    ON CONFLICT (user_id, channel_id) DO UPDATE 
                    SET source = EXCLUDED.source, issuer_id = EXCLUDED.issuer_id
                ''', member_object.id, static_channel_id, mute_source, author_id)
                if mute_source == "owner":
                    await conn.execute('''
                        INSERT INTO users (user_id, server_mute_channel_ids)
                        VALUES ($1, ARRAY[$2]::BIGINT[])
                        ON CONFLICT (user_id) DO UPDATE
                        SET server_mute_channel_ids = (
                            SELECT ARRAY(
                                SELECT DISTINCT unnest(COALESCE(u.server_mute_channel_ids, '{}') || ARRAY[$2])
                            )
                            FROM users u WHERE u.user_id = EXCLUDED.user_id
                        ),
                        updated_at = NOW()
                    ''', member_object.id, static_channel_id)
                else:
                    await conn.execute('''
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
                    ''', member_object.id, static_channel_id)
                await conn.execute('''
                    INSERT INTO mute_reasons (guild_id, user_id, reason, channel_id)
                    VALUES ($1, $2, $3, $4)
                    ON CONFLICT (guild_id, user_id, channel_id)
                    DO UPDATE SET reason = EXCLUDED.reason
                ''', guild_id, member_object.id, reason, static_channel_id)
            if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                await member_object.edit(mute=True)
                mute_type = "server muted" if mute_source == "owner" else "muted"
                await ctx.send(f'{member_object.mention} has been {mute_type} in <#{static_channel_id}> with reason {reason}.'),
                allowed_mentions=discord.AllowedMentions.none()
            else:
                mute_type = "server muted" if mute_source == "owner" else "muted"
                await ctx.send(f'{member_object.mention} has been {mute_type} in <#{static_channel_id}> with reason {reason}.'),
                allowed_mentions=discord.AllowedMentions.none()
        return voice_mute_alias

    def create_unban_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(
            name=command_name,
            help='Unban a user from a voice channel.'
        )
        @commands.check(is_owner_developer_coordinator_moderator)
        async def unban_alias(
                ctx,
                member: str = commands.parameter(description='Mention or user ID of the member to unban.'),
                *,
                reason: str = commands.parameter(default='N/A', description='Optional reason for unbanning.')
        ) -> None:
            guild_id = ctx.guild.id
            member_id = None
            if member.isdigit():
                member_id = int(member)
            elif member.startswith('<@') and member.endswith('>'):
                try:
                    member_id = int(member.strip('<@!>'))
                except ValueError:
                    pass
            member_object = ctx.guild.get_member(member_id) if member_id else None
            if not member_object:
                return await self.handler.send_message(ctx, content='❌ Could not resolve a valid member.')
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('unban', {}).get(command_name)
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
                        guild_id, command_name
                    )
                if not static_channel_id:
                    return await self.handler.send_message(ctx, content=f'❌ No channel alias mapping found for `{command_name}`.')
            channel = ctx.guild.get_channel(static_channel_id)
            if not channel or not isinstance(channel, discord.VoiceChannel):
                return await self.handler.send_message(ctx, content=f'❌ Could not resolve a valid voice channel for <#{static_channel_id}>.')
            try:
                await channel.set_permissions(
                    member_object,
                    overwrite=None,
                    reason=f"Unbanned from <#{channel.id}>: {reason or 'No reason provided'}"
                )
            except discord.Forbidden:
                return await self.handler.send_message(ctx, content='❌ Missing permissions to update channel permissions.')
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute(
                    'DELETE FROM active_bans WHERE user_id = $1 AND channel_id = $2',
                    member_object.id, channel.id
                )
                await conn.execute(
                    'DELETE FROM ban_expirations WHERE user_id = $1 AND channel_id = $2',
                    member_object.id, channel.id
                )
                await conn.execute(
                    '''
                    UPDATE users
                    SET ban_channel_ids = array_remove(ban_channel_ids, $2),
                        updated_at      = NOW()
                    WHERE user_id = $1
                    ''',
                    member_object.id, channel.id
                )
            await ctx.send(f'🔊 {member_object.mention} has been unbanned from <#{channel.id}>.',
                           allowed_mentions=discord.AllowedMentions.none()
                           )

        return unban_alias

    def create_unflag_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(
            name=command_name,
            help='Unflag a user in the database for the voice channel mapped to this alias.'
        )
        @commands.check(is_owner_developer_coordinator_moderator)
        async def unflag_alias(
                ctx,
                user: str = commands.parameter(description='Tag a user or include their snowflake ID.')
        ) -> None:
            guild_id = ctx.guild.id
            flag_aliases = self.bot.command_aliases.get(guild_id, {}).get('unflag', {})
            channel_id = flag_aliases.get(command_name)
            if not channel_id:
                await self.handler.send_message(ctx, content=f'❌ No unflag alias configured for `{command_name}`.')
                return
            if not user:
                await self.handler.send_message(ctx, content='❌ You must provide a user ID or mention.')
                return
            if re.fullmatch(r'<@!?\d+>', user):
                user_id = int(re.sub(r'\D', '', user))
            elif user.isdigit():
                user_id = int(user)
            else:
                await self.handler.send_message(ctx, content='❌ Invalid user ID or mention.')
                return
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
                    is_flagged = await conn.fetchval(select_sql, user_id, channel_id)
                    if not is_flagged:
                        await self.handler.send_message(ctx, content=f'ℹ️ <@{user_id}> is not flagged for <#{channel_id}>.')
                        return
                    await conn.execute(update_sql, user_id, channel_id)
                    await ctx.send(
                        f'✅ Unflagged <@{user_id}> for channel <#{channel_id}>.',
                        allowed_mentions=discord.AllowedMentions.none()
                    )
            except Exception as e:
                await self.handler.send_message(ctx, content=f'❌ Database error: {e}')
                raise
        return unflag_alias

    def create_unmute_alias(self, command_name: str) -> Command:
        @commands.hybrid_command(name=command_name, help='Unmutes a member in a specific VC.')
        @commands.check(is_owner_developer_coordinator_moderator)
        async def unmute_alias(
            ctx,
            member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Include a reason for the unmute.')
        ) -> None:
            guild_id = ctx.guild.id
            static_channel_id = self.bot.command_aliases.get(guild_id, {}).get('unmute', {}).get(command_name)
            if not static_channel_id:
                return await self.handler.send_message(ctx, content=f'❌ No unmute alias configured for {command_name}.')
            member_id = None
            member_object = None
            if member.isdigit():
                member_id = int(member)
            elif member.startswith('<@') and member.endswith('>'):
                try:
                    member_id = int(member.strip('<@!>'))
                except ValueError:
                    pass
            if member_id:
                member_object = ctx.guild.get_member(member_id)
            if not member_object:
                return await self.handler.send_message(ctx, content='❌ Could not resolve a valid guild member from your input.')
            try:
                async with self.bot.db_pool.acquire() as conn:
                    row = await conn.fetchrow('''
                        SELECT source, issuer_id FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2
                    ''', member_object.id, static_channel_id)
                    if not row:
                        return await ctx.send(f"❌ {member_object.mention} is not muted in <#{static_channel_id}>.", allowed_mentions=discord.AllowedMentions.none())
                    bot_owner_id = int(os.environ.get("DISCORD_OWNER_ID", "0"))
                    server_owner_id = ctx.guild.owner_id
                    command_author_id = ctx.author.id
                    if row["source"] == "owner":
                        if command_author_id == bot_owner_id:
                            await ctx.send(f"⚠️ Overriding server mute for {member_object.mention} as bot owner.", allowed_mentions=discord.AllowedMentions.none())
                        elif command_author_id == server_owner_id:
                            return await ctx.send(f"❌ You cannot unmute {member_object.mention} — server mutes can't be overridden by the server owner.", allowed_mentions=discord.AllowedMentions.none())
                        else:
                            return await ctx.send(f"❌ {member_object.mention} has a server mute in <#{static_channel_id}> that cannot be removed.")
                    elif row["source"] not in ("bot", "manual", "bot_owner"):
                        return await ctx.send(f"❌ {member_object.mention} was not muted by the bot in <#{static_channel_id}>.", allowed_mentions=discord.AllowedMentions.none())
                    if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                        await member_object.edit(mute=False)
                    await conn.execute('''
                        DELETE FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2 AND source IN ('bot', 'manual', 'bot_owner')
                    ''', member_object.id, static_channel_id)
                    if row["source"] == "bot" or row["source"] == "bot_owner":
                        await conn.execute('''
                            UPDATE users
                            SET mute_channel_ids = array_remove(mute_channel_ids, $2),
                                updated_at = NOW()
                            WHERE user_id = $1
                        ''', member_object.id, static_channel_id)
                    elif row["source"] == "manual":
                        await conn.execute('''
                            UPDATE users
                            SET manual_mute_channels = array_remove(manual_mute_channels, $2),
                                updated_at = NOW()
                            WHERE user_id = $1
                        ''', member_object.id, static_channel_id)
                    await conn.execute('''
                        INSERT INTO mute_reasons (guild_id, user_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4)
                        ON CONFLICT (guild_id, user_id, channel_id)
                        DO UPDATE SET reason = EXCLUDED.reason
                    ''', guild_id, member_object.id, static_channel_id, f"Unmuted: {reason}")
            except Exception as e:
                await self.handler.send_message(ctx, content=f'❌ Database error: {e}')
                raise
            if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                await ctx.send(f'✅ {member_object.mention} has been unmuted in <#{static_channel_id}>.', allowed_mentions=discord.AllowedMentions.none())
            else:
                await ctx.send(f'✅ {member_object.mention} is no longer marked as muted in <#{static_channel_id}>.', allowed_mentions=discord.AllowedMentions.none())
        return unmute_alias
        
    @commands.hybrid_command(name='xalias', help='Deletes an alias.')
    @commands.check(is_owner_developer_coordinator)
    async def delete_alias(
        self,
        ctx,
        alias_name: str = commands.parameter(description='Includ an alias name')
    ) -> None:
        if not alias_name.strip():
            await self.handler.send_message(ctx, '❌ `alias_name` cannot be empty.', ephemeral=True)
            return
        guild_id = ctx.guild.id
        alias_type = None
        for candidate in ('mute', 'unmute', 'ban', 'unban', 'flag'):
            if alias_name in self.bot.command_aliases.get(guild_id, {}).get(candidate, {}):
                alias_type = candidate
                break
        if alias_type.lower() not in {'mute', 'unmute', 'ban', 'unban', 'flag'}:
            await self.handler.send_message(ctx, '❌ `alias_type` must be either `mute`, `unmute`, `ban`, `unban`, or `flag`.', ephemeral=True)
            return
        if not alias_type:
            await self.handler.send_message(ctx, f'❌ Alias `{alias_name}` not found.', ephemeral=True)
            return
        alias_map = self.bot.command_aliases.get(guild_id, {}).get(alias_type.lower(), {})
        if alias_name not in alias_map:
            await self.handler.send_message(ctx, f'❌ Alias `{alias_name}` not found in `{alias_type}` for guild `{guild_id}`.', ephemeral=True)
            return
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                'DELETE FROM command_aliases WHERE guild_id = $1 AND alias_type = $2 AND alias_name = $3',
                guild_id, alias_type.lower(), alias_name
            )
        if self.bot.get_command(alias_name):
            self.bot.remove_command(alias_name)
        self.bot.command_aliases[guild_id][alias_type.lower()].pop(alias_name, None)
        await self.handler.send_message(ctx, content=f'✅ Deleted alias `{alias_name}` from `{alias_type}`.')
            
    @commands.hybrid_command(
        name='xcoord',
        help='Revokes coordinator access from a user in a specific voice channel.'
    )
    @commands.check(is_owner_developer)
    async def delete_coordinator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel: discord.VoiceChannel = commands.parameter(description='Voice channel to revoke coordinator access from.')
    ) -> None:
        member_id = None
        member_object = None
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='❌ Could not resolve a valid guild member from your input.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                "SELECT coordinator_ids, coordinator_channel_ids FROM users WHERE user_id = $1",
                member_object.id
            )
            if not row:
                return ctx.send(f'❌ {member_object.mention} is not found in the coordinator database.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            current_guild_ids = row.get('coordinator_ids', []) or []
            current_channel_ids = row.get('coordinator_channel_ids', []) or []
            if ctx.guild.id not in current_guild_ids:
                return ctx.send(f'❌ {member_object.mention} is not a coordinator in this guild.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            if channel.id not in current_channel_ids:
                return await ctx.send(f'❌ {member_object.mention} is not a coordinator in {channel.mention}.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            await conn.execute('''
                UPDATE users
                SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member_object.id, channel.id)
            updated_row = await conn.fetchrow(
                "SELECT coordinator_channel_ids FROM users WHERE user_id = $1",
                member_object.id
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
                ''', member_object.id, ctx.guild.id)
                await ctx.send(f'✅ {member_object.mention}\'s coordinator access has been completely revoked from {channel.mention} and this guild (no remaining channels).', allowed_mentions=discord.AllowedMentions.none())
            else:
                await ctx.send(f'✅ {member_object.mention}\'s coordinator access has been revoked from {channel.mention}.', allowed_mentions=discord.AllowedMentions.none())

    @commands.hybrid_command(name='xdev', help='Removes a developer.')
    @commands.check(is_owner)
    async def delete_developer(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ) -> None:
        member_id = None
        member_object = None
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        guild_id = ctx.guild.id
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET developer_guild_ids = array_remove(developer_guild_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member_object.id, guild_id)
        await ctx.send(f'{member_object.mention}\'s developer access has been revoked in this guild.', allowed_mentions=discord.AllowedMentions.none())
        
    @commands.hybrid_command(name='xmod', help='Revokes a member\'s VC moderator role for a given channel.')
    @commands.check(is_owner_developer_coordinator)
    async def delete_moderator(
        self,
        ctx,
        member: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel: str = commands.parameter(description='Tag a VC or include its snowflake ID.')
    ) -> None:
        member_id = None
        member_object = None
        if member.isdigit():
            member_id = int(member)
        elif member.startswith('<@') and member.endswith('>'):
            try:
                member_id = int(member.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        resolved_channel = None
        if channel.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel))
        elif channel.startswith('<#') and channel.endswith('>'):
            try:
                channel_id = int(channel.strip('<#>'))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel.lower():
                    resolved_channel = vc
                    break
        if not isinstance(resolved_channel, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content='Could not resolve a valid **voice** channel from your input.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                "SELECT moderator_ids, moderator_channel_ids FROM users WHERE user_id = $1",
                member_object.id
            )
            if not row:
                return await ctx.send(f'{member_object.mention} is not found in the moderator database.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            current_guild_ids = row.get('moderator_ids', []) or []
            current_channel_ids = row.get('moderator_channel_ids', []) or []
            if ctx.guild.id not in current_guild_ids:
                return await ctx.send(f'{member_object.mention} is not a moderator in this guild.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            if resolved_channel.id not in current_channel_ids:
                return await ctx.send(f'{member_object.mention} is not a moderator in {resolved_channel.name}.',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            await conn.execute('''
                UPDATE users
                SET moderator_channel_ids = array_remove(moderator_channel_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member_object.id, resolved_channel.id)
            updated_row = await conn.fetchrow(
                "SELECT moderator_channel_ids FROM users WHERE user_id = $1",
                member_object.id
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
                ''', member_object.id, ctx.guild.id)
                await ctx.send(f'{member_object.mention} has been completely revoked as moderator from {resolved_channel.name} and this guild (no remaining channels).',
                    allowed_mentions=discord.AllowedMentions.none()
                )
            else:
                await ctx.send(f'{member_object.mention} has been revoked moderator access in {resolved_channel.name}.',
                    allowed_mentions=discord.AllowedMentions.none()
                )


    @commands.hybrid_command(name='aliases', help='List all the aliases in the current guild, or filter by channel.')
    @commands.check(is_owner_developer_coordinator_moderator)
    async def list_aliases(self, ctx, channel: str = None) -> None:
        guild_id = ctx.guild.id
        aliases = self.bot.command_aliases.get(guild_id, {})
        if not aliases:
            await self.handler.send_message(ctx, content='No aliases defined in this guild.')
            return
        pages = []
        if channel:
            channel_obj = None
            channel_id = None
            mention_match = re.match(r'<#(\d+)>', channel)
            if mention_match:
                channel_id = int(mention_match.group(1))
            elif channel.isdigit():
                channel_id = int(channel)
            else:
                for ch in ctx.guild.text_channels:
                    if ch.name.lower() == channel.lower():
                        channel_id = ch.id
                        break
            if channel_id:
                channel_obj = ctx.guild.get_channel(channel_id)
            if not channel_obj:
                await self.handler.send_message(
                    ctx,
                    content=f'Channel "{channel}" not found. Please use a channel mention, ID, or name.'
                )
                return
            found_aliases = False
            for kind in ('mute', 'unmute', 'ban', 'unban', 'flag'):
                entries = aliases.get(kind, {})
                if not entries:
                    continue
                channel_entries = {name: cid for name, cid in entries.items() if cid == channel_id}
                if channel_entries:
                    found_aliases = True
                    embed = discord.Embed(
                        title=f'{kind.capitalize()} Aliases for {channel_obj.mention}',
                        description='\n'.join(f'`{name}` → <#{cid}>' for name, cid in channel_entries.items()),
                        color=discord.Color.blue()
                    )
                    pages.append(embed)
            if not found_aliases:
                await self.handler.send_message(
                    ctx,
                    content=f'No aliases found for {channel_obj.mention}.'
                )
                return
        else:
            for kind in ('mute', 'unmute', 'ban', 'unban', 'flag'):
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
            await self.handler.send_message(ctx, content='No aliases found.')
            return
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()
        
    @commands.hybrid_command(name='bans', help='Lists all users banned from a specific channel.')
    @commands.check(is_owner_developer_coordinator_moderator)
    async def list_bans(
        self,
        ctx: commands.Context,
        channel: str = commands.parameter(description='Mention or ID of the channel to check bans for.')
    ) -> None:
        guild = ctx.guild
        channel_id = None
        if channel.isdigit():
            channel_id = int(channel)
        elif channel.startswith('<#') and channel.endswith('>'):
            try:
                channel_id = int(channel.strip('<#>'))
            except ValueError:
                pass
        if not channel_id or not guild.get_channel(channel_id):
            return await self.handler.send_message(ctx, content='❌ Could not resolve a valid channel.')
        async with self.bot.db_pool.acquire() as conn:
            bans = await conn.fetch('''
                SELECT ab.user_id, ab.expires_at, br.reason
                FROM active_bans ab
                LEFT JOIN ban_reasons br ON ab.user_id = br.user_id AND ab.channel_id = br.channel_id AND br.guild_id = $1
                WHERE ab.channel_id = $2
                ORDER BY ab.expires_at NULLS LAST
            ''', guild.id, channel_id)
        if not bans:
            return await self.handler.send_message(ctx, content=f'✅ No active bans found for <#{channel_id}>.')
        embeds = []
        for record in bans:
            user = guild.get_member(record['user_id'])
            name = user.mention if user else f"User ID <!{record['user_id']}>"
            reason = record['reason'] or 'No reason provided'
            expires_at = record['expires_at']
            duration_str = (
                "Permanent" if expires_at is None
                else discord.utils.format_dt(expires_at, style='R')
            )
            embed = discord.Embed(
                title="⛔ Channel Ban",
                description=f"{name} is banned from <#{channel_id}>.",
                color=discord.Color.red()
            )
            embed.add_field(name="Reason", value=reason, inline=False)
            embed.add_field(name="Duration", value=duration_str, inline=True)
            embed.set_footer(text=f"Channel ID: {channel_id}")
            embeds.append(embed)
        paginator = Paginator(self.bot, ctx, embeds)
        await paginator.start()
    
    @commands.hybrid_command(name='coords', help='Lists coordinators for a specific voice channel.')
    @commands.check(is_owner_developer_coordinator)
    async def list_coordinators(
        self,
        ctx,
        channel: str = commands.parameter(description='Voice channel ID, mention, or name.')
    ) -> None:
        await self.list_role_members(
            ctx=ctx,
            channel_str=channel,
            column_name="coordinator_channel_ids",
            label="Coordinator",
            color=discord.Color.gold()
        )
    
    @commands.hybrid_command(name='devs', help='Lists all developers.')
    @commands.check(is_owner_developer)
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
            await self.handler.send_message(ctx, content='No developers are configured in this guild.')
            return
        for row in rows:
            user_id = row['user_id']
            user = guild.get_member(user_id)
            name = user.display_name if user else f'User ID {user_id}'
            embed = discord.Embed(
                title=f'Developer: {name}',
                color=discord.Color.blue()
            )
            pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()
        
    @commands.hybrid_command(
        name='flags',
        help='List users flagged in the database for a specific voice channel.'
    )
    @commands.check(is_owner_developer_coordinator_moderator)
    async def list_flags(
        self,
        ctx,
        channel: str = commands.parameter(description='Voice channel ID, mention, or name.')
    ) -> None:
        guild = ctx.guild
        def resolve_channel(value: str) -> Optional[discord.VoiceChannel]:
            if value.isdigit():
                channel = guild.get_channel(int(value))
            elif value.startswith('<#') and value.endswith('>'):
                channel = guild.get_channel(int(value.strip('<#>')))
            else:
                channel = discord.utils.get(guild.voice_channels, name=value)
            return channel if isinstance(channel, discord.VoiceChannel) else None
        channel = resolve_channel(channel)
        if not channel:
            return await self.handler.send_message(ctx, content='❌ Could not resolve a valid voice channel.')
        channel_id = channel.id
        sql = '''
            SELECT user_id FROM users
            WHERE $1 = ANY(flagged_channel_ids)
        '''
        try:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(sql, channel_id)
            if not rows:
                return await self.handler.send_message(ctx, content=f'✅ No users are flagged for <#{channel_id}>.')
            mentions = [f'<@{row["user_id"]}>' for row in rows]
            message = f'🚩 Users flagged for <#{channel_id}>:\n' + '\n'.join(mentions)
            await ctx.send(message, allowed_mentions=discord.AllowedMentions.none())
        except Exception as e:
            await self.handler.send_message(ctx, content=f'❌ Database error: {e}')
            raise
            
    @commands.hybrid_command(name='mods', help='Lists moderators for a specific voice channel.')
    @commands.check(is_owner_developer_coordinator_moderator)
    async def list_moderators(
        self,
        ctx,
        channel: str = commands.parameter(description='Voice channel ID, mention, or name.')
    ) -> None:
        await self.list_role_members(
            ctx=ctx,
            channel_str=channel,
            column_name="moderator_channel_ids",
            label="Moderator",
            color=discord.Color.green()
        )
        
    @commands.hybrid_command(
        name='mutes',
        help='Lists all users currently muted in a specific voice channel.'
    )
    @commands.check(is_owner_developer_coordinator)
    async def list_mutes(
        self,
        ctx,
        channel: str = commands.parameter(description='Voice channel mention or ID.')
    ) -> None:
        guild_id = ctx.guild.id
        try:
            if isinstance(channel, discord.VoiceChannel):
                channel_id = channel.id
            elif channel.isdigit():
                channel_id = int(channel)
            elif channel.startswith('<#') and channel.endswith('>'):
                channel_id = int(channel.strip('<#>'))
            else:
                raise ValueError
        except ValueError:
            return await self.handler.send_message(ctx, content='❌ Invalid channel mention or ID.')
        async with self.bot.db_pool.acquire() as conn:
            records = await conn.fetch('''
                SELECT am.user_id, mr.reason, am.source
                FROM active_mutes am
                JOIN mute_reasons mr ON am.user_id = mr.user_id AND am.channel_id = mr.channel_id
                WHERE am.channel_id = $1 AND mr.guild_id = $2
            ''', channel_id, guild_id)
        if not records:
            return await self.handler.send_message(ctx, content=f'ℹ️ No users are currently muted in <#{channel_id}>.')
        description_lines = []
        for record in records:
            user = ctx.guild.get_member(record['user_id'])
            mention = user.mention if user else f'`{record["user_id"]}`'
            reason = record['reason']
            description_lines.append(f'• {mention} — {reason}')
        embed = discord.Embed(
            title=f'Muted Users in <#{channel_id}>',
            description='\n'.join(description_lines),
            color=discord.Color.orange()
        )
        await ctx.send(embed=embed, allowed_mentions=discord.AllowedMentions.none())

    # @commands.hybrid_command(name='reason', help='Get the reason for a mute, unmute, ban, or unban.')
    # @commands.check(is_owner_developer_coordinator_moderator)
    # async def get_summary(
    #     ctx,
    #     member: str = commands.parameter(description='Tag a user or include their snowflake ID.')
    # ) -> None:
    #     guild_id = ctx.guild.id
    #     member_id = None
    #     member_object = None
    #     if member.isdigit():
    #         member_id = int(member)
    #     elif member.startswith('<@') and member.endswith('>'):
    #         try:
    #             member_id = int(member.strip('<@!>'))
    #         except ValueError:
    #             pass
    #     if member_id:
    #         member_object = ctx.guild.get_member(member_id)
    #     if not member_object:
    #         return await self.handler.send_message(ctx, content='❌ Could not resolve a valid guild member from your input.')
    #     async with self.bot.db_pool.acquire() as conn:
    #         mute_rows = await conn.fetch('''
    #             SELECT channel_id, reason, action
    #             FROM mute_reasons
    #             WHERE guild_id = $1 AND user_id = $2
    #         ''', guild_id, member_object.id)
    #         ban_rows = await conn.fetch('''
    #             SELECT channel_id, reason, action
    #             FROM ban_reasons
    #             WHERE guild_id = $1 AND user_id = $2
    #         ''', guild_id, member_object.id)
    #     if not mute_rows and not ban_rows:
    #         return await ctx.send(f'ℹ️ No mute, unmute, ban, or unban history found for {member_object.mention}.', allowed_mentions=discord.AllowedMentions.none())
    #     lines = []
    #     def format_row(row, kind):
    #         channel_id = row['channel_id']
    #         reason = row['reason'] or '*No reason provided*'
    #         action = row.get('action', kind)
    #         return f'• [{action}] <#{channel_id}>: `{reason}`'
    #     lines += [format_row(row, 'mute') for row in mute_rows]
    #     lines += [format_row(row, 'ban') for row in ban_rows]
    #     content = f'📄 Disciplinary reasons for {member_object.mention}:\n' + '\n'.join(lines)
    #     await ctx.send(content, allowed_mentions=discord.AllowedMentions.none())
    
    async def cog_load(self) -> None:
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT guild_id, alias_type, alias_name, channel_id FROM command_aliases')
            for row in rows:
                guild_id = row['guild_id']
                alias_type = row['alias_type']
                alias_name = row['alias_name']
                channel_id = row['channel_id']
                self.bot.command_aliases.setdefault(guild_id, {}).setdefault(alias_type, {})[alias_name] = channel_id
                cmd = None
                if alias_type == 'mute':
                    cmd = self.create_voice_mute_alias(alias_name)
                elif alias_type == 'unmute':
                    cmd = self.create_unmute_alias(alias_name)
                elif alias_type == 'ban':
                    cmd = self.create_ban_alias(alias_name)
                elif alias_type == 'unban':
                    cmd = self.create_unban_alias(alias_name)
                elif alias_type == 'flag':
                    cmd = self.create_flag_alias(alias_name)
                elif alias_type == 'unflag':
                    cmd = self.create_unflag_alias(alias_name)
                if cmd:
                    self.bot.add_command(cmd)
                    
    async def list_role_members(
        self,
        ctx: commands.Context,
        *,
        channel_str: str,
        column_name: str,
        label: str,
        color: discord.Color,
    ) -> None:
        guild = ctx.guild
        voice_channel = self.resolve_voice_channel(guild, channel_str)
        if not voice_channel:
            return await self.handler.send_message(ctx, content='❌ Could not resolve a valid voice channel.')
        async with self.bot.db_pool.acquire() as conn:
            query = f'''
                SELECT user_id FROM users
                WHERE $1 = ANY({column_name})
            '''
            rows = await conn.fetch(query, voice_channel.id)
        if not rows:
            return await self.handler.send_message(
                ctx, content=f'ℹ️ No {label.lower()}s found for <#{voice_channel.id}>.'
            )
        pages = []
        for row in rows:
            user_id = row['user_id']
            user = guild.get_member(user_id)
            display_name = user.display_name if user else f'User ID {user_id}'
            embed = discord.Embed(
                title=f'{label} for {voice_channel.name}',
                description=f'{display_name} (<@{user_id}>)',
                color=color
            )
            embed.set_footer(text=f'User ID: {user_id}')
            pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()
    
    def resolve_voice_channel(self, guild: discord.Guild, value: str) -> Optional[discord.VoiceChannel]:
        try:
            if value.isdigit():
                channel = guild.get_channel(int(value))
            elif value.startswith('<#') and value.endswith('>'):
                channel_id_str = value[2:-1]
                if channel_id_str.isdigit():
                    channel = guild.get_channel(int(channel_id_str))
                else:
                    return None
            else:
                channel = discord.utils.get(guild.voice_channels, name=value)
            return channel if isinstance(channel, discord.VoiceChannel) else None
        except (ValueError, AttributeError):
            return None
    
async def setup(bot: commands.Bot):
    cog = Hybrid(bot)
    await bot.add_cog(cog)
