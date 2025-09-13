''' event_listeners.py

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
import discord
from datetime import datetime, timezone
import os
from discord.ext import commands
from vyrtuous.inc.helpers import *
from vyrtuous.utils.setup_logging import logger
from vyrtuous.service.check_service import *
from vyrtuous.service.discord_message_service import DiscordMessageService
from vyrtuous.bot.discord_bot import DiscordBot
import asyncio

class EventListeners(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.db_pool)
    
    async def fetch_active_bans(conn, guild_id: int, user_id: int):
        return await conn.fetch('''
            SELECT channel_id, expires_at
            FROM active_bans
            WHERE guild_id = $1 AND discord_snowflake = $2
              AND (expires_at IS NULL OR expires_at > NOW())
        ''', guild_id, user_id)

    async def fetch_active_text_mutes(conn, guild_id: int, user_id: int):
        return await conn.fetch('''
            SELECT channel_id
            FROM active_text_mutes
            WHERE guild_id = $1 AND discord_snowflake = $2
              AND (expires_at IS NULL OR expires_at > NOW())
        ''', guild_id, user_id)
        
    @commands.Cog.listener()
    async def on_voice_state_update(self, member: discord.Member, before: discord.VoiceState, after: discord.VoiceState) -> None:
        if before.channel == after.channel and before.mute == after.mute and before.self_mute == after.self_mute:
            return
        if member.bot:
            return
        after_channel = after.channel
        before_channel = before.channel
        user_id = member.id
        should_be_muted = False
        just_manual_unmute = False
        async with self.db_pool.acquire() as conn:
            user_data = await conn.fetchrow('SELECT server_mute_guild_ids FROM users WHERE discord_snowflake = $1', user_id)
            server_mute_guild_ids = user_data['server_mute_guild_ids'] or [] if user_data else []
            if member.guild.id in server_mute_guild_ids:
                return
            if before.mute != after.mute:
                if before.mute and not after.mute and before_channel:
                    await conn.execute('DELETE FROM active_voice_mutes WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3', member.guild.id, user_id, before_channel.id)
                    just_manual_unmute = True
            if after_channel:
                row = await conn.fetchrow('SELECT expires_at FROM active_voice_mutes WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3', member.guild.id, user_id, after_channel.id)
                if row:
                    if row['expires_at'] and row['expires_at'] <= datetime.now(timezone.utc):
                        await conn.execute('DELETE FROM active_voice_mutes WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3', member.guild.id, user_id, after_channel.id)
                        row = None
                        should_be_muted = False
                if not before.mute and after.mute and after_channel:
                    existing_row = await conn.fetchrow('SELECT guild_id FROM active_voice_mutes WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3', member.guild.id, user_id, after_channel.id)
                    if not existing_row:
                        await conn.execute('INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id) VALUES ($1, $2, $3) ON CONFLICT (guild_id, discord_snowflake, channel_id) DO NOTHING', member.guild.id, user_id, after_channel.id)
                mute_row = await conn.fetchrow('''
                    SELECT expires_at
                    FROM active_voice_mutes
                    WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                ''', member.guild.id, user_id, after_channel.id)
                should_be_muted = False
                if mute_row:
                    if not mute_row['expires_at'] or mute_row['expires_at'] > datetime.now(timezone.utc):
                        should_be_muted = True
                if just_manual_unmute:
                    should_be_muted = False
                if should_be_muted and not after.mute:
                    await conn.execute('INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id) VALUES ($1, $2, $3) ON CONFLICT (guild_id, discord_snowflake, channel_id) DO NOTHING', member.guild.id, user_id, after_channel.id)
                    try:
                        await member.edit(mute=True, reason=f'Enforcing mute in {after_channel.name} (found in arrays)')
                    except discord.Forbidden:
                        logger.debug(f'No permission to mute {member.display_name}')
                    except discord.HTTPException as e:
                        logger.debug(f'Failed to mute {member.display_name}: {e}')
                elif not should_be_muted and after.mute:
                    try:
                        await member.edit(mute=False, reason=f'Auto-unmuting in {after_channel.name} (no mute record)')
                    except discord.Forbidden:
                        logger.debug(f'No permission to unmute {member.display_name}')
                    except discord.HTTPException as e:
                        logger.debug(f'Failed to unmute {member.display_name}: {e}')
                if before.channel is None and after_channel is not None:
                    try:
                        is_flagged = await conn.fetchval('''
                            SELECT 1
                            FROM active_flags
                            WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
                        ''', member.guild.id, user_id, after_channel.id)
                        if is_flagged and isinstance(after_channel, discord.VoiceChannel):
                            await after_channel.send(f'⚠️ <@{user_id}> has joined voice channel <#{after_channel.id}> and is flagged.', allowed_mentions=discord.AllowedMentions.none())
                    except Exception as e:
                        logger.exception('Error in on_voice_state_update', exc_info=e)
                        
    @commands.Cog.listener()
    async def on_member_join(self, member: discord.Member) -> None:
        user_id = member.id
        guild = member.guild
        async with self.db_pool.acquire() as conn:
            bans = await fetch_active_bans(conn, user_id)
            text_mutes = await fetch_active_text_mutes(conn, user_id)
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
    async def on_command_error(self, ctx, error):
        if isinstance(error, commands.MissingRequiredArgument):
            missing = error.param.name
            await ctx.reply(f'❌ Missing required argument: `{missing}`')
            return
        if isinstance(error, commands.CheckFailure):
            return await send_check_failure_embed(ctx, error)
    
async def setup(bot: DiscordBot):
    await bot.add_cog(EventListeners(bot))
