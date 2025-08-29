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
        source = 'bot'
        just_manual_unmute = False
        async with self.db_pool.acquire() as conn:
            user_data = await conn.fetchrow(
                "SELECT server_mute_guild_ids FROM users WHERE user_id = $1",
                user_id
            )
            server_mute_guild_ids = user_data['server_mute_guild_ids'] or [] if user_data else []
            if member.guild.id in server_mute_guild_ids:
                return
            if before.mute != after.mute:
                if before.mute and not after.mute and before_channel:
                    row = await conn.fetchrow("SELECT source FROM active_mutes WHERE user_id = $1 AND channel_id = $2", user_id, before_channel.id)
                    await conn.execute("DELETE FROM active_mutes WHERE user_id = $1 AND channel_id = $2", user_id, before_channel.id)
                    if row and row['source'] == 'bot_owner':
                        await conn.execute("UPDATE users SET server_mute_guild_ids = array_remove(server_mute_guild_ids, $2), updated_at = NOW() WHERE user_id = $1", user_id,     member.guild.id)
                    else:
                        await conn.execute("UPDATE users SET manual_mute_channels = array_remove(manual_mute_channels, $2), mute_channel_ids = array_remove(mute_channel_ids, $2),     updated_at = NOW() WHERE user_id = $1", user_id, before_channel.id)
                        just_manual_unmute = True
            if after_channel:
                row = await conn.fetchrow("SELECT source, expires_at FROM active_mutes WHERE user_id = $1 AND channel_id = $2", user_id, after_channel.id)
                if row:
                    if row['expires_at'] and row['expires_at'] <= datetime.now(timezone.utc):
                        logger.debug(f"Mute expired for {member.display_name} in {after_channel.name}")
                        await conn.execute("DELETE FROM active_mutes WHERE user_id = $1 AND channel_id = $2", user_id, after_channel.id)
                        if row['source'] == 'owner':
                            await conn.execute("UPDATE users SET server_mute_guild_ids = array_remove(server_mute_guild_ids, $2), updated_at = NOW() WHERE user_id = $1", user_id,     member.guild.id)
                        else:
                            await conn.execute("UPDATE users SET manual_mute_channels = array_remove(manual_mute_channels, $2), mute_channel_ids = array_remove(mute_channel_ids, $2),     updated_at = NOW() WHERE user_id = $1", user_id, after_channel.id)
                        row = None
                        should_be_muted = False
                if not before.mute and after.mute and after_channel:
                    existing_row = await conn.fetchrow("SELECT source FROM active_mutes WHERE user_id = $1 AND channel_id = $2", user_id, after_channel.id)
                    if not existing_row:
                        await conn.execute("INSERT INTO active_mutes (user_id, channel_id, source, issuer_id) VALUES ($1, $2, 'manual', $3) ON CONFLICT (user_id, channel_id) DO NOTHING", user_id, after_channel.id, member.guild.owner_id)
                        await conn.execute(
                            "INSERT INTO users (user_id, manual_mute_channels) VALUES ($1, ARRAY[$2]::BIGINT[]) ON CONFLICT (user_id) DO UPDATE SET manual_mute_channels = (SELECT     ARRAY(SELECT DISTINCT unnest(COALESCE (u.manual_mute_channels, '{}') || ARRAY[$2])) FROM users u WHERE u.user_id = EXCLUDED.user_id), updated_at = NOW()",
                            user_id, after_channel.id
                        )
                mute_row = await conn.fetchrow("SELECT source FROM active_mutes WHERE user_id = $1 AND channel_id = $2", user_id, after_channel.id)
                if mute_row:
                    should_be_muted = True
                    source = mute_row['source']
                else:
                    user_mute_data = await conn.fetchrow("SELECT mute_channel_ids, server_mute_guild_ids, manual_mute_channels FROM users WHERE user_id = $1", user_id)
                    if user_mute_data:
                        mute_channels = user_mute_data['mute_channel_ids'] or []
                        server_mute_guild_ids = user_mute_data['server_mute_guild_ids'] or []
                        manual_mute_channels = user_mute_data['manual_mute_channels'] or []
                        if str(member.id) == os.getenv("DISCORD_OWNER_ID"):
                            should_be_muted = True
                            source = 'bot_owner'
                        elif member.guild.id in server_mute_guild_ids:
                            should_be_muted = True
                            source = 'owner'
                        elif after_channel.id in manual_mute_channels:
                            should_be_muted = True
                            source = 'manual'
                        elif after_channel.id in mute_channels:
                            should_be_muted = True
                            source = 'bot'
                if just_manual_unmute:
                    should_be_muted = False
                if should_be_muted and not after.mute:
                    await conn.execute("INSERT INTO active_mutes (user_id, channel_id, source, issuer_id) VALUES ($1, $2, $3, $4) ON CONFLICT (user_id, channel_id) DO NOTHING", user_id, after_channel.id, source, member.guild.owner_id)
                    try:
                        await member.edit(mute=True, reason=f"Enforcing mute in {after_channel.name} (found in arrays)")
                    except discord.Forbidden:
                        logger.debug(f"No permission to mute {member.display_name}")
                    except discord.HTTPException as e:
                        logger.debug(f"Failed to mute {member.display_name}: {e}")
                elif not should_be_muted and after.mute:
                    try:
                        await member.edit(mute=False, reason=f"Auto-unmuting in {after_channel.name} (no mute record)")
                    except discord.Forbidden:
                        logger.debug(f"No permission to unmute {member.display_name}")
                    except discord.HTTPException as e:
                        logger.debug(f"Failed to unmute {member.display_name}: {e}")
                if before.channel is None and after.channel is not None:
                    try:
                        is_flagged = await conn.fetchval("SELECT 1 FROM users WHERE user_id = $1 AND $2 = ANY (flagged_channel_ids)", member.id, after_channel.id)
                        if not is_flagged:
                            return
                        if isinstance(after_channel, discord.VoiceChannel):
                            await after_channel.send(f'⚠️ <@{member.id}> has joined voice channel <#{after_channel.id}> and is flagged.', allowed_mentions=discord.AllowedMentions.none())
                    except Exception as e:
                        logger.exception("Error in on_voice_state_update", exc_info=e)
                        
    @commands.Cog.listener()
    async def on_member_join(self, member: discord.Member) -> None:
        user_id = member.id
        guild = member.guild
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch("""
                SELECT channel_id, expires_at
                FROM active_bans
                WHERE user_id = $1
            """, user_id)
        for row in rows:
            channel_id = row['channel_id']
            expires_at = row['expires_at']
            if expires_at and expires_at < datetime.now(timezone.utc):
                continue
            channel = guild.get_channel(channel_id)
            if not channel or not isinstance(channel, (discord.TextChannel, discord.VoiceChannel)):
                continue
            try:
                overwrite = channel.overwrites_for(member)
                overwrite.view_channel = False
                await channel.set_permissions(member, overwrite=overwrite, reason="Reinstating active channel ban")
            except discord.Forbidden:
                print(f"Missing permissions to ban in channel {channel_id}")
            except discord.HTTPException as e:
                print(f"Failed to apply ban for {member} in {channel_id}: {e}")

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
            await ctx.reply(f"❌ Missing required argument: `{missing}`")
            return
        if isinstance(error, commands.CheckFailure):
            return await send_check_failure_embed(ctx, error)
    
async def setup(bot: DiscordBot):
    await bot.add_cog(EventListeners(bot))


