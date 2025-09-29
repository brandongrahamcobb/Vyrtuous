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
from vyrtuous.service.discord_message_service import DiscordMessageService, ChannelPaginator
from vyrtuous.bot.discord_bot import DiscordBot
import asyncio

class EventListeners(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.db_pool)
    
    async def fetch_active_bans(self, conn, guild_id: int, user_id: int):
        return await conn.fetch('''
            SELECT channel_id, expires_at
            FROM active_bans
            WHERE guild_id = $1 AND discord_snowflake = $2
              AND (expires_at IS NULL OR expires_at > NOW())
        ''', guild_id, user_id)

    async def fetch_active_text_mutes(self, conn, guild_id: int, user_id: int):
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
                    expires_at = mute_row['expires_at'] if mute_row else None
                    if expires_at:
                        await conn.execute('''
                            INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at)
                            VALUES ($1, $2, $3, $4)
                            ON CONFLICT (guild_id, discord_snowflake, channel_id)
                            DO UPDATE SET expires_at = EXCLUDED.expires_at
                        ''', member.guild.id, user_id, after_channel.id, expires_at)
                    else:
                        await conn.execute('''
                            INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id)
                            VALUES ($1, $2, $3)
                            ON CONFLICT (guild_id, discord_snowflake, channel_id) DO NOTHING
                        ''', member.guild.id, user_id, after_channel.id)
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
                        rows = await conn.fetch('''
                            SELECT channel_id, discord_snowflake, reason
                            FROM active_flags
                            WHERE guild_id = $1 AND discord_snowflake = $2
                        ''', member.guild.id, user_id)
                        if not rows:
                            return
                        grouped = {}
                        for row in rows:
                            grouped.setdefault(row['channel_id'], []).append(row)
                        embeds = []
                        embed = discord.Embed(
                            title=f'\u26A0\uFE0F Flags for {member.display_name}',
                            color=discord.Color.red()
                        )
                        embed.set_thumbnail(url=member.display_avatar.url)
                        context_records = grouped.get(after_channel.id, [])
                        for record in context_records:
                            reason = record['reason'] or 'No reason provided'
                            embed.add_field(name=f'Channel: {after_channel.mention}', value=f'Reason: {reason}', inline=False)
                        other_channels = [ch_id for ch_id in grouped.keys() if ch_id != after_channel.id]
                        if other_channels:
                            ch_mentions = []
                            for ch_id in other_channels:
                                ch = member.guild.get_channel(ch_id)
                                ch_mentions.append(ch.mention if ch else f'Channel ID `{ch_id}`')
                            embed.add_field(name='Other channels', value='\n'.join(ch_mentions), inline=False)
                        embeds.append(embed)
                        for ch_id in other_channels:
                            records = grouped[ch_id]
                            ch = member.guild.get_channel(ch_id)
                            ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                            embed = discord.Embed(
                                title=f'\u26A0\uFE0F Flags for {member.display_name} in {ch_name}',
                                color=discord.Color.red()
                            )
                            embed.set_thumbnail(url=member.display_avatar.url)
                            for record in records:
                                reason = record['reason'] or 'No reason provided'
                                embed.add_field(name='Channel', value=f'{ch_name}\nReason: {reason}', inline=False)
                            embeds.append(embed)
                        paginator = ChannelPaginator(self.bot, after_channel, embeds)
                        await paginator.start()
                    except Exception as e:
                        logger.exception('Error in on_voice_state_update', exc_info=e)
                        
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
    async def on_command_error(self, ctx, error):
        if isinstance(error, commands.MissingRequiredArgument):
            missing = error.param.name
            await ctx.reply(f'\U0001F6AB Missing required argument: `{missing}`')
            return
        if isinstance(error, commands.CheckFailure):
            return await send_check_failure_embed(ctx, error)
    
async def setup(bot: DiscordBot):
    await bot.add_cog(EventListeners(bot))
