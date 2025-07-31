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
import os
from discord.ext import commands
from vyrtuous.inc.helpers import *

from vyrtuous.service.discord_message_service import DiscordMessageService


class EventListeners(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.db_pool)

    @commands.Cog.listener()
    async def on_voice_state_update(
        self,
        member: discord.Member,
        before: discord.VoiceState,
        after: discord.VoiceState
    ) -> None:
        if member.bot:
            return
    
        user_id = member.id
        before_channel = before.channel
        after_channel = after.channel
    
        if os.getenv("DEVELOPMENT") == "False":
            async with self.db_pool.acquire() as conn:
                # Handle when someone gets unmuted (was muted, now unmuted)
                if before.mute and not after.mute and before_channel:
                    row = await conn.fetchrow("""
                        SELECT source, expires_at FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2
                    """, user_id, before_channel.id)
    
                    if row and row['source'] == 'manual':
                        # Remove manual mutes when manually unmuted
                        await conn.execute("""
                            DELETE FROM active_mutes
                            WHERE user_id = $1 AND channel_id = $2
                        """, user_id, before_channel.id)
    
                        await conn.execute("""
                            UPDATE users
                            SET manual_mute_channels = array_remove(manual_mute_channels, $2),
                                updated_at = NOW()
                            WHERE user_id = $1
                        """, user_id, before_channel.id)
                    elif row and row['source'] in ('bot', 'owner', 'bot_owner'):
                        # Convert bot/owner mutes to 'unmuted' when manually unmuted
                        await conn.execute("""
                            UPDATE active_mutes
                            SET source = 'unmuted'
                            WHERE user_id = $1 AND channel_id = $2
                        """, user_id, before_channel.id)
    
                # Handle when someone gets muted (was unmuted, now muted)
                if not before.mute and after.mute and after_channel:
                    row = await conn.fetchrow("""
                        SELECT source FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2
                    """, user_id, after_channel.id)
    
                    if not row:
                        # No record exists - create manual mute
                        await conn.execute("""
                            INSERT INTO active_mutes (user_id, channel_id, source, issuer_id)
                            VALUES ($1, $2, 'manual', $3)
                            ON CONFLICT (user_id, channel_id) DO NOTHING
                        """, user_id, after_channel.id, member.guild.owner_id)
    
                        await conn.execute("""
                            INSERT INTO users (user_id, manual_mute_channels)
                            VALUES ($1, ARRAY[$2]::BIGINT[])
                            ON CONFLICT (user_id) DO UPDATE
                            SET manual_mute_channels = (
                                SELECT ARRAY(
                                    SELECT DISTINCT unnest(COALESCE(u.manual_mute_channels, '{}') || ARRAY[$2])
                                )
                                FROM users u WHERE u.user_id = EXCLUDED.user_id
                            ),
                            updated_at = NOW()
                        """, user_id, after_channel.id)
                    elif row['source'] == 'unmuted':
                        # Convert unmuted back to manual when manually muted
                        await conn.execute("""
                            UPDATE active_mutes
                            SET source = 'manual'
                            WHERE user_id = $1 AND channel_id = $2
                        """, user_id, after_channel.id)
    
                        await conn.execute("""
                            INSERT INTO users (user_id, manual_mute_channels)
                            VALUES ($1, ARRAY[$2]::BIGINT[])
                            ON CONFLICT (user_id) DO UPDATE
                            SET manual_mute_channels = (
                                SELECT ARRAY(
                                    SELECT DISTINCT unnest(COALESCE(u.manual_mute_channels, '{}') || ARRAY[$2])
                                )
                                FROM users u WHERE u.user_id = EXCLUDED.user_id
                            ),
                            updated_at = NOW()
                        """, user_id, after_channel.id)
    
                # Handle when someone joins/switches to a channel - enforce mute state
                if after_channel:
                    print(f"DEBUG: {member.display_name} joined {after_channel.name} (ID: {after_channel.id})")
    
                    row = await conn.fetchrow("""
                        SELECT source, expires_at FROM active_mutes
                        WHERE user_id = $1 AND channel_id = $2
                    """, user_id, after_channel.id)
    
                    print(f"DEBUG: Active mute record for {member.display_name} in {after_channel.name}: {row}")
    
                    # Check if mute has expired
                    if row and row['expires_at'] and row['expires_at'] <= datetime.utcnow():
                        print(f"DEBUG: Mute expired for {member.display_name} in {after_channel.name}")
                        # Mute has expired, clean it up
                        await conn.execute("""
                            DELETE FROM active_mutes
                            WHERE user_id = $1 AND channel_id = $2
                        """, user_id, after_channel.id)
    
                        # Remove from user arrays based on source
                        if row['source'] == 'owner':
                            await conn.execute("""
                                UPDATE users
                                SET server_mute_channel_ids = array_remove(server_mute_channel_ids, $2),
                                    updated_at = NOW()
                                WHERE user_id = $1
                            """, user_id, after_channel.id)
                        elif row['source'] == 'manual':
                            await conn.execute("""
                                UPDATE users
                                SET manual_mute_channels = array_remove(manual_mute_channels, $2),
                                    updated_at = NOW()
                                WHERE user_id = $1
                            """, user_id, after_channel.id)
                        else:
                            await conn.execute("""
                                UPDATE users
                                SET mute_channel_ids = array_remove(mute_channel_ids, $2),
                                    updated_at = NOW()
                                WHERE user_id = $1
                            """, user_id, after_channel.id)
    
                        row = None  # Clear the row since mute expired
    
                    if row:
                        print(f"DEBUG: Found active mute record, source: {row['source']}")
                        # Has active mute record - enforce based on source
                        if row['source'] in ('manual', 'bot', 'owner', 'bot_owner') and not after.mute:
                            try:
                                await member.edit(mute=True, reason=f"Enforcing mute in {after_channel.name}")
                                print(f"Enforced mute on {member.display_name} in {after_channel.name}")
                            except discord.Forbidden:
                                print(f"No permission to mute {member.display_name}")
                            except discord.HTTPException as e:
                                print(f"Failed to mute {member.display_name}: {e}")
                        elif row['source'] == 'unmuted' and after.mute:
                            try:
                                await member.edit(mute=False, reason=f"Enforcing unmute in {after_channel.name}")
                                print(f"Enforced unmute on {member.display_name} in {after_channel.name}")
                            except discord.Forbidden:
                                print(f"No permission to unmute {member.display_name}")
                            except discord.HTTPException as e:
                                print(f"Failed to unmute {member.display_name}: {e}")
                    else:
                        print(f"DEBUG: No active mute record found, checking user arrays")
                        # No active record - check user arrays for missing records
                        user_mute_data = await conn.fetchrow("""
                            SELECT mute_channel_ids, server_mute_channel_ids, manual_mute_channels
                            FROM users 
                            WHERE user_id = $1
                        """, user_id)
    
                        print(f"DEBUG: User mute data for {member.display_name}: {user_mute_data}")
    
                        should_be_muted = False
                        source = 'bot'
    
                        if user_mute_data:
                            mute_channels = user_mute_data['mute_channel_ids'] or []
                            server_mute_channels = user_mute_data['server_mute_channel_ids'] or []
                            manual_mute_channels = user_mute_data['manual_mute_channels'] or []
    
                            print(f"DEBUG: mute_channels: {mute_channels}")
                            print(f"DEBUG: server_mute_channels: {server_mute_channels}")
                            print(f"DEBUG: manual_mute_channels: {manual_mute_channels}")
                            print(f"DEBUG: Current channel ID: {after_channel.id}")
    
                            if after_channel.id in server_mute_channels:
                                should_be_muted = True
                                source = 'owner'
                                print(f"DEBUG: Found in server_mute_channels")
                            elif after_channel.id in manual_mute_channels:
                                should_be_muted = True
                                source = 'manual'
                                print(f"DEBUG: Found in manual_mute_channels")
                            elif after_channel.id in mute_channels:
                                should_be_muted = True
                                source = 'bot'
                                print(f"DEBUG: Found in mute_channels")
    
                        print(f"DEBUG: should_be_muted: {should_be_muted}")
    
                        if should_be_muted:
                            # Create missing active_mutes record
                            await conn.execute("""
                                INSERT INTO active_mutes (user_id, channel_id, source, issuer_id)
                                VALUES ($1, $2, $3, $4)
                                ON CONFLICT (user_id, channel_id) DO NOTHING
                            """, user_id, after_channel.id, source, member.guild.owner_id)
    
                            print(f"DEBUG: Created missing active_mutes record")
    
                            # ALWAYS enforce the mute regardless of current state
                            try:
                                await member.edit(mute=True, reason=f"Enforcing mute in {after_channel.name} (found in arrays)")
                                print(f"Enforced mute on {member.display_name} in {after_channel.name} (found in arrays)")
                            except discord.Forbidden:
                                print(f"No permission to mute {member.display_name}")
                            except discord.HTTPException as e:
                                print(f"Failed to mute {member.display_name}: {e}")
                        elif after.mute:
                            # User is muted but shouldn't be - unmute them
                            try:
                                await member.edit(mute=False, reason=f"Auto-unmuting in {after_channel.name} (no mute record)")
                                print(f"Auto-unmuted {member.display_name} in {after_channel.name} (no mute record)")
                            except discord.Forbidden:
                                print(f"No permission to unmute {member.display_name}")
                            except discord.HTTPException as e:
                                print(f"Failed to unmute {member.display_name}: {e}")
#            async with self.db_pool.acquire() as conn:
#                # Handle leaving a voice channel
#                if before_channel and not after_channel:
#                    row = await conn.fetchrow("""
#                        SELECT source FROM active_mutes
#                        WHERE user_id = $1 AND channel_id = $2
#                    """, user_id, before_channel.id)
#    
#                    if row and row['source'] == 'manual':
#                        await conn.execute("""
#                            DELETE FROM active_mutes
#                            WHERE user_id = $1 AND channel_id = $2
#                        """, user_id, before_channel.id)
#    
#                        # Update the correct column based on your schema
#                        await conn.execute("""
#                            UPDATE users
#                            SET mute_channel_ids = array_remove(mute_channel_ids, $2),
#                                updated_at = NOW()
#                            WHERE user_id = $1
#                        """, user_id, before_channel.id)
#    
#                # Handle joining or switching into a voice channel
#                if after_channel:
#                    # Check for any active mute in this channel
#                    row = await conn.fetchrow("""
#                        SELECT source, expires_at FROM active_mutes
#                        WHERE user_id = $1 AND channel_id = $2
#                    """, user_id, after_channel.id)
#    
#                    # Check if mute has expired
#                    if row and row['expires_at'] and row['expires_at'] <= datetime.utcnow():
#                        # Mute has expired, clean it up
#                        await conn.execute("""
#                            DELETE FROM active_mutes
#                            WHERE user_id = $1 AND channel_id = $2
#                        """, user_id, after_channel.id)
#    
#                        # Remove from user arrays based on source
#                        if row['source'] == 'owner':
#                            await conn.execute("""
#                                UPDATE users
#                                SET server_mute_channel_ids = array_remove(server_mute_channel_ids, $2),
#                                    updated_at = NOW()
#                                WHERE user_id = $1
#                            """, user_id, after_channel.id)
#                        else:
#                            await conn.execute("""
#                                UPDATE users
#                                SET mute_channel_ids = array_remove(mute_channel_ids, $2),
#                                    updated_at = NOW()
#                                WHERE user_id = $1
#                            """, user_id, after_channel.id)
#    
#                        row = None  # Clear the row since mute expired
#    
#                    # Case 1: The channel is muted for this user (and not expired)
#                    if row and row['source'] in ('manual', 'bot', 'owner', 'bot_owner'):
#                        # ALWAYS enforce the mute, regardless of current state
#                        try:
#                            await member.edit(mute=True, reason=f"Enforcing mute in {after_channel.name}")
#                            print(f"Enforced mute on {member.display_name} in {after_channel.name}")
#                        except discord.Forbidden:
#                            print(f"No permission to mute {member.display_name}")
#                        except discord.HTTPException as e:
#                            print(f"Failed to mute {member.display_name}: {e}")
#    
#                    # Case 2: The channel is explicitly unmuted for this user
#                    elif row and row['source'] == 'unmuted':
#                        if after.mute:
#                            # User joined muted but should be unmuted - unmute them
#                            try:
#                                await member.edit(mute=False, reason=f"Auto-unmuting in {after_channel.name} (explicitly unmuted)")
#                                print(f"Auto-unmuted {member.display_name} in {after_channel.name} (explicitly unmuted)")
#                            except discord.Forbidden:
#                                print(f"No permission to unmute {member.display_name}")
#                            except discord.HTTPException as e:
#                                print(f"Failed to unmute {member.display_name}: {e}")
#                        # Otherwise: already unmuted, do nothing
#    
#                    # Case 3: No mute state is stored → check if they should be unmuted
#                    elif not row:
#                        if after.mute:
#                            # Only unmute if they're NOT in any mute arrays in the users table
#                            user_mute_data = await conn.fetchrow("""
#                                SELECT mute_channel_ids, server_mute_channel_ids 
#                                FROM users 
#                                WHERE user_id = $1
#                            """, user_id)
#    
#                            should_unmute = True
#                            if user_mute_data:
#                                mute_channels = user_mute_data['mute_channel_ids'] or []
#                                server_mute_channels = user_mute_data['server_mute_channel_ids'] or []
#    
#                                # If this channel is in either mute array, don't unmute
#                                if after_channel.id in mute_channels or after_channel.id in server_mute_channels:
#                                    should_unmute = False
#    
#                            if should_unmute:
#                                try:
#                                    await member.edit(mute=False, reason=f"Auto-unmuting in {after_channel.name} (no mute record)")
#                                    print(f"Auto-unmuted {member.display_name} in {after_channel.name} (no mute record)")
#                                except discord.Forbidden:
#                                    print(f"No permission to unmute {member.display_name}")
#                                except discord.HTTPException as e:
#                                    print(f"Failed to unmute {member.display_name}: {e}")
#                                # They should be muted but have no active_mutes record - create one and enforce mute
#                                # Determine source based on which array contains the channel
#                                if user_mute_data and user_mute_data['server_mute_channel_ids'] and after_channel.id in user_mute_data['server_mute_channel_ids']:
#                                    source = 'owner'
#                                else:
#                                    source = 'bot'  # or whatever default you prefer
#    
#                                await conn.execute("""
#                                    INSERT INTO active_mutes (user_id, channel_id, source, issuer_id)
#                                    VALUES ($1, $2, $3, $4)
#                                    ON CONFLICT (user_id, channel_id) DO NOTHING
#                                """, user_id, after_channel.id, source, member.guild.owner_id)
#    
#                                # ALWAYS enforce the mute
#                                try:
#                                    await member.edit(mute=True, reason=f"Enforcing mute in {after_channel.name} (found in mute arrays)")
#                                    print(f"Enforced mute on {member.display_name} in {after_channel.name} (found in mute arrays)")
#                                except discord.Forbidden:
#                                    print(f"No permission to mute {member.display_name}")
#                                except discord.HTTPException as e:
#                                    print(f"Failed to mute {member.display_name}: {e}")
                    # If not muted and no record, that's the expected state - do nothing
                    # If not muted and no record, that's the expected state - do nothing
    
                try:
                    async with self.bot.db_pool.acquire() as conn:
                        is_flagged = await conn.fetchval(
                            '''
                            SELECT 1 FROM users
                            WHERE user_id = $1
                              AND $2 = ANY(flagged_channel_ids)
                            ''',
                            member.id, after_channel.id
                        )
                        if not is_flagged:
                            return
                        if isinstance(channel, discord.VoiceChannel):
                            await channel.send(
                            f'⚠️ <@{member.id}> has joined voice channel <#{channel.id}> and is flagged.',
                            allowed_mentions=discord.AllowedMentions.none()
                        )
                except Exception as e:
                    print(f"🔥 Error in on_voice_state_update: {e}")

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
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

    
async def setup(bot: commands.Bot):
    await bot.add_cog(EventListeners(bot))

