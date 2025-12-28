
''' aliases.py A discord.py cog containing command aliases for the Vyrtuous bot.

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
from datetime import datetime, timedelta, timezone
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.ban import Ban
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.duration import Duration
from vyrtuous.utils.flag import Flag
from vyrtuous.utils.reason import Reason
from vyrtuous.utils.statistics import Statistics
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.state import State
from vyrtuous.utils.invincibility import Invincibility
from vyrtuous.utils.voice_mute import VoiceMute
from vyrtuous.utils.emojis import Emojis

import time

class Aliases(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.alias_help = {
            'ban': [
                "**member** (Optional): Tag a member or include their snowflake ID",
                "**duration** (Optional): (+|-)duration(m|h|d)\n0 = permanent / 24h = default\n`+` to append, `-` to delete, `=` to overwrite reason",
                "**reason** (Optional): Optional reason (required for 7 days or more)"
            ],
            'cow': [
                "**member** (Optional): Tag a member or include their snowflake ID"
            ],
            'uncow': [
                "**member** (Optional): Tag a member or include their snowflake ID"
            ],
            'unban': [
                "**member** (Optional): Tag a member or include their snowflake ID"
            ],
            'flag': [
                "**member** (Optional): Tag a member or include their snowflake ID",
                "**reason** (Optional): Optional reason for flagging the user"
            ],
            'unflag': [
                "**member** (Optional): Tag a member or include their snowflake ID"
            ],
            'mute': [
                "**member** (Optional): Tag a member or include their snowflake ID",
                "**duration** (Optional): (+|-)duration(m|h|d)\n0 = permanent / 24h = default\n`+` to append, `-` to delete, `=` to overwrite reason",
                "**reason** (Optional): Optional reason (required for 7 days or more)"
            ],
            'unmute': [
                "**member** (Optional): Tag a member or include their snowflake ID"
            ],
            'tmute': [
                "**member** (Optional): Tag a member or include their snowflake ID",
                "**duration** (Optional): (+|-)duration(m|h|d)\n0 = permanent / 24h = default\n`+` to append, `-` to delete, `=` to overwrite reason",
                "**reason** (Optional): Optional reason (required for 7 days or more)"
            ],
            'untmute': [
                "**member** (Optional): Tag a member or include their snowflake ID"
            ],
            'role': [
                "**member** (Optional): Tag a member or include their snowflake ID",
                "**role** (Optional): Role to assign"
            ],
            'unrole': [
                "**member** (Optional): Tag a member or include their snowflake ID",
                "**role** (Optional): Role to remove"
            ]
        }
        self.alias_type_to_description = {
            'ban': 'Bans a user from the server.',
            'cow': 'Verifies a user as going vegan.',
            'uncow': 'Unverifies a user as going vegan.',
            'unban': 'Unbans a user from the server.',
            'flag': 'Flags a user for moderation review.',
            'unflag': 'Removes a flag from a user.',
            'mute': 'Mutes a user in voice channels.',
            'unmute': 'Unmutes a user in voice channels.',
            'tmute': 'Mutes a user in text channels.',
            'untmute': 'Unmutes a user in text channels.',
            'role': 'Assigns a role to a user.',
            'unrole': 'Removes a role from a user.'
        }
        self.alias_type_to_permission_level = {
            'ban': 'Moderator',
            'cow': 'Moderator',
            'uncow': 'Moderator',
            'unban': 'Moderator',
            'mute': 'Moderator',
            'unmute': 'Moderator',
            'tmute': 'Moderator',
            'untmute': 'Moderator',
            'flag': 'Moderator',
            'unflag': 'Moderator',
            'role': 'Coordinator',
            'unrole': 'Coordinator'
        }
        self.bot = bot
        self.emoji = Emojis()
        self.channel_service = ChannelService()
        self.member_service = MemberService()
        self.invincible_members = Invincibility.get_invincible_members()
    
    async def handle_ban_alias(self, alias, args, channel_obj, executor_role, existing_moderation, is_duration_modification, is_reason_modification, member_obj, message, state):
        embed = discord.Embed(title=f"{self.emoji.get_random_emoji()} Ban", color=discord.Color.orange())
        if is_reason_modification and existing_moderation:
            modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
            match is_reason_modification:
                case '+':
                    updated_reason = existing_moderation.reason + modified_reason
                case '=' | '-':
                    updated_reason = modified_reason
            await Alias.update_reason(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, updated_reason=updated_reason)
            await Statistics.send_statistic(message, alias.alias_type, member_obj, channel_obj, str(updated_reason), existing_moderation.expires_at, alias.alias_name, None, True, executor_role)
            try:
                return await state.end(success=f"{self.emoji.get_random_emoji()} {member_obj.display_name} ban has been updated with reason {updated_reason}.")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        if is_duration_modification:
            if existing_moderation:
                duration = Duration(args[1] if len(args) > 1 else '24h')
                match str(duration)[0]:
                    case '+':
                        updated_expires_at = existing_moderation.expires_at + duration.to_timedelta()
                    case '=':
                        updated_expires_at = datetime.now() + duration.to_timedelta()
                    case '-':
                        updated_expires_at = existing_moderation.expires_at - duration.to_timedelta()
                await Alias.update_duration(channel_snowflake=channel_obj.id, expires_at=updated_expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)
            await Statistics.send_statistic(message, alias.alias_type, member_obj, channel_obj, str(updated_reason), updated_expires_at, alias.alias_name, None, True, executor_role)
            try:
                return await state.end(success=f"{self.emoji.get_random_emoji()} {member_obj.display_name} ban has been updated with duration {updated_expires_at}.")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        duration = Duration(args[1] if len(args) > 1 else '24h')
        reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'
        try:
            await channel_obj.set_permissions(member_obj, view_channel=False, reason=reason)
        except discord.Forbidden:
            return await state.end(error=f'\U0001F3C6 {e}')
        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            try:
                await member_obj.move_to(None, reason=reason)
            except discord.Forbidden:
                return await state.end(error=f'\U0001F3C6 {e}')
            is_in_channel = True
        else:
            is_in_channel = False
        ban = Ban(channel_snowflake=channel_obj.id, expires_at=duration.expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason)
        await ban.create()
        await Statistics.send_statistic(message, 'ban', member_obj, channel_obj, str(duration), str(reason), duration.expires_at, alias.alias_name, is_in_channel, False, executor_role)
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} has been banned",
            description=f"**User:** {member_obj.mention}\n**Channel:** {channel_obj.mention}\n**Duration:** {str(duration)}\n**Reason:** {reason}",
            color=discord.Color.orange()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
    
    # DONE
    async def handle_cow_alias(self, message: discord.Message, alias: Alias, args):
        highest_role = await is_owner_developer_administrator_coordinator_moderator(message)
        if highest_role == 'Everyone':
            return await message.reply(content='\U000026A0\U0000FE0F You are not permitted to cow users.')
        member = args[0] if len(args) > 0 else None
        if not message.guild:
            return await message.reply(content='\U000026A0\U0000FE0F This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj:
            return await message.reply(content=f'\U000026A0\U0000FE0F Could not resolve a valid member from input: {member}.')
        if member_obj.id == message.guild.me.id:
            return await message.reply(content='\U000026A0\U0000FE0F You cannot cow the bot.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel_snowflake)
        allowed, highest_role = await has_equal_or_higher_role(message, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, sender_snowflake=message.author.id)
        if not allowed:
            return await message.reply(content=f'\U000026A0\U0000FE0F You are not allowed to cow this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        select_sql = '''
            SELECT 1
            FROM active_cows
            WHERE discord_snowflake = $1
            AND channel_id = $2
        '''
        insert_cow_sql = '''
            INSERT INTO active_cows (guild_id, discord_snowflake, channel_id, created_at)
            VALUES ($1, $2, $3, $4)
            ON CONFLICT (guild_id, discord_snowflake, channel_id) DO NOTHING
        '''
        insert_log_sql = '''
            INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
            VALUES ($1, $2, $3, $4, $5, $6)
            RETURNING created_at
        '''
        async with self.bot.db_pool.acquire() as conn:
            already_cowed = await conn.fetchval(select_sql, member_obj.id, channel_obj.id)
            if already_cowed:
                return await message.reply(content=f'\U000026A0\U0000FE0F {member_obj.mention} is already going vegan.')
            created_at = await conn.fetchval(insert_log_sql, 'cow', member_obj.id, message.author.id, message.guild.id, channel_obj.id, 'Cowed a user')
            await conn.execute(insert_cow_sql, message.guild.id, member_obj.id, channel_obj.id, created_at)
            await message.reply(content=f'\U0001F525 {member_obj.mention} is going vegan!!! \U0001F525')
    
    # DONE
    async def handle_flag_alias(self, alias, args, channel_obj, executor_role, existing_moderation, is_duration_modification, is_reason_modification, member_obj, message, state):
        embed = discord.Embed(title=f"{self.emoji.get_random_emoji()} Flag", color=discord.Color.orange())
        if is_reason_modification and existing_moderation:
            modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
            match is_reason_modification:
                case '+':
                    updated_reason = existing_moderation.reason + modified_reason
                case '=' | '-':
                    updated_reason = modified_reason
            await Alias.update_reason(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, updated_reason=updated_reason)
            await Statistics.send_statistic(message, alias.alias_type, member_obj, channel_obj, str(updated_reason), existing_moderation.expires_at, alias.alias_name, None, True, executor_role)
            try:
                return await state.end(success=f"{self.emoji.get_random_emoji()} {member_obj.display_name} flag has been updated with reason {updated_reason}.")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        if is_duration_modification:
            if existing_moderation:
                duration = Duration(args[1] if len(args) > 1 else '24h')
                match str(duration)[0]:
                    case '+':
                        updated_expires_at = existing_moderation.expires_at + duration.to_timedelta()
                    case '=':
                        updated_expires_at = datetime.now() + duration.to_timedelta()
                    case '-':
                        updated_expires_at = existing_moderation.expires_at - duration.to_timedelta()
                await Alias.update_duration(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, expires_at=updated_expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)
            await Statistics.send_statistic(message, alias.alias_type, member_obj, channel_obj, str(updated_reason), updated_expires_at, alias.alias_name, None, True, executor_role)
            try:
                return await state.end(success=f"{self.emoji.get_random_emoji()} {member_obj.display_name} flag has been updated with duration {updated_expires_at}.")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        duration = Duration(args[1] if len(args) > 1 else '24h')
        reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'
        try:
            await channel_obj.set_permissions(member_obj, view_channel=False, reason=reason)
        except discord.Forbidden:
            return await state.end(error=f'\U0001F3C6 {e}')
        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            try:
                await member_obj.move_to(None, reason=reason)
            except discord.Forbidden:
                return await state.end(error=f'\U0001F3C6 {e}')
            is_in_channel = True
        else:
            is_in_channel = False
        ban = Ban(channel_snowflake=channel_obj.id, expires_at=duration.expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason)
        await ban.create()
        await Statistics.send_statistic(message, 'flag', member_obj, channel_obj, str(duration), str(reason), duration.expires_at, alias.alias_name, is_in_channel, False, executor_role)
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} has been flagged",
            description=f"**User:** {member_obj.mention}\n**Channel:** {channel_obj.mention}\n**Duration:** {str(duration)}\n**Reason:** {reason}",
            color=discord.Color.orange()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
   
    # DONE
    async def handle_role_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_moderation,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        role_obj = message.guild.get_role(alias.role_snowflake)
        if not role_obj:
            return await state.end(warning=f"\U000026A0\U0000FE0F Could not resolve role with ID `{alias.role_snowflake}`.")
        if role_obj in member_obj.roles:
            return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} already has {role_obj.mention}.')
        try:
            await member_obj.add_roles(role_obj, reason='Added role')
        except discord.Forbidden:
            return await state.end(error=f'\U0001F3C6 {member_obj.mention} was not successfully roled.')
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} Role Added",
            description=(
                f"**User:** {member_obj.mention}\n"
                f"**Role:** {role_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')

    # DONE
    async def handle_text_mute_alias(self, alias, args, channel_obj, executor_role, existing_moderation, is_duration_modification, is_reason_modification, member_obj, message, state):
        embed = discord.Embed(title=f"{self.emoji.get_random_emoji()} Flag", color=discord.Color.orange())
        if is_reason_modification and existing_moderation:
            modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
            match is_reason_modification:
                case '+':
                    updated_reason = existing_moderation.reason + modified_reason
                case '=' | '-':
                    updated_reason = modified_reason
            await Alias.update_reason(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, updated_reason=updated_reason)
            await Statistics.send_statistic(message, alias.alias_type, member_obj, channel_obj, str(updated_reason), existing_moderation.expires_at, alias.alias_name, None, True, executor_role)
            try:
                return await state.end(success=f"{self.emoji.get_random_emoji()} {member_obj.display_name} text-mute has been updated with reason {updated_reason}.")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        if is_duration_modification:
            if existing_moderation:
                duration = Duration(args[1] if len(args) > 1 else '24h')
                match str(duration)[0]:
                    case '+':
                        updated_expires_at = existing_moderation.expires_at + duration.to_timedelta()
                    case '=':
                        updated_expires_at = datetime.now() + duration.to_timedelta()
                    case '-':
                        updated_expires_at = existing_moderation.expires_at - duration.to_timedelta()
                await Alias.update_duration(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, expires_at=updated_expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)
            await Statistics.send_statistic(message, alias.alias_type, member_obj, channel_obj, str(updated_reason), updated_expires_at, alias.alias_name, None, True, executor_role)
            try:
                return await state.end(success=f"{self.emoji.get_random_emoji()} {member_obj.display_name} text-mute has been updated with duration {updated_expires_at}.")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        duration = Duration(args[1] if len(args) > 1 else '24h')
        reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'
        try:
            await channel_obj.set_permissions(member_obj, view_channel=False, reason=reason)
        except discord.Forbidden:
            return await state.end(error=f'\U0001F3C6 {e}')
        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            try:
                await member_obj.move_to(None, reason=reason)
            except discord.Forbidden:
                return await state.end(error=f'\U0001F3C6 {e}')
            is_in_channel = True
        else:
            is_in_channel = False
        ban = Ban(channel_snowflake=channel_obj.id, expires_at=duration.expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason)
        await ban.create()
        await Statistics.send_statistic(message, 'tmute', member_obj, channel_obj, str(duration), str(reason), duration.expires_at, alias.alias_name, is_in_channel, False, executor_role)
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} has been text-muted",
            description=f"**User:** {member_obj.mention}\n**Channel:** {channel_obj.mention}\n**Duration:** {str(duration)}\n**Reason:** {reason}",
            color=discord.Color.orange()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
        
    # DONE
    async def handle_voice_mute_alias(self, alias, args, channel_obj, executor_role, existing_moderation, is_duration_modification, is_reason_modification, member_obj, message, state):
        embed = discord.Embed(title=f"{self.emoji.get_random_emoji()} Flag", color=discord.Color.orange())
        if is_reason_modification and existing_moderation:
            modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
            match is_reason_modification:
                case '+':
                    updated_reason = existing_moderation.reason + modified_reason
                case '=' | '-':
                    updated_reason = modified_reason
            await Alias.update_reason(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, updated_reason=updated_reason)
            await Statistics.send_statistic(message, alias.alias_type, member_obj, channel_obj, str(updated_reason), existing_moderation.expires_at, alias.alias_name, None, True, executor_role)
            try:
                return await state.end(success=f"{self.emoji.get_random_emoji()} {member_obj.display_name} voice-mute has been updated with reason {updated_reason}.")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        if is_duration_modification:
            if existing_moderation:
                duration = Duration(args[1] if len(args) > 1 else '24h')
                match str(duration)[0]:
                    case '+':
                        updated_expires_at = existing_moderation.expires_at + duration.to_timedelta()
                    case '=':
                        updated_expires_at = datetime.now() + duration.to_timedelta()
                    case '-':
                        updated_expires_at = existing_moderation.expires_at - duration.to_timedelta()
                await Alias.update_duration(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, expires_at=updated_expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)
            await Statistics.send_statistic(message, alias.alias_type, member_obj, channel_obj, str(updated_reason), updated_expires_at, alias.alias_name, None, True, executor_role)
            try:
                return await state.end(success=f"{self.emoji.get_random_emoji()} {member_obj.display_name} voice-mute has been updated with duration {updated_expires_at}.")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        duration = Duration(args[1] if len(args) > 1 else '24h')
        reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'
        try:
            await channel_obj.set_permissions(member_obj, view_channel=False, reason=reason)
        except discord.Forbidden:
            return await state.end(error=f'\U0001F3C6 {e}')
        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            try:
                await member_obj.move_to(None, reason=reason)
            except discord.Forbidden:
                return await state.end(error=f'\U0001F3C6 {e}')
            is_in_channel = True
        else:
            is_in_channel = False
        ban = Ban(channel_snowflake=channel_obj.id, expires_at=duration.expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason)
        await ban.create()
        await Statistics.send_statistic(message, 'vmute', member_obj, channel_obj, str(duration), str(reason), duration.expires_at, alias.alias_name, is_in_channel, False, executor_role)
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} has been voice-muted",
            description=f"**User:** {member_obj.mention}\n**Channel:** {channel_obj.mention}\n**Duration:** {str(duration)}\n**Reason:** {reason}",
            color=discord.Color.orange()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
        
    # DONE
    async def handle_unban_alias(self, alias, args, channel_obj, executor_role, existing_moderation, is_duration_modification, is_reason_modification, member_obj, message, state
    ):
        bans = await Ban.fetch_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        if not bans:
            return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently banned in {channel_obj.mention}.')
        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            try:
                await member_obj.move_to(None, reason="Unbanned")
            except discord.Forbidden:
                return await state.end(error=f'\U0001F3C6 Could not move {member_obj.mention} from the voice channel.')
        for ban in bans:
            if ban.expires_at is None and executor_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator'):
                return await state.end(warning='\U000026A0\U0000FE0F Only coordinators and above can undo permanent bans.')
            try:
                await channel_obj.set_permissions(member_obj, overwrite=None)
            except discord.Forbidden:
                return await state.end(error=f'\U0001F3C6 {member_obj.mention} was not successfully unbanned.')
        await Ban.delete_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} has been unbanned",
            description=(
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Unbanned By:** {message.author.mention}"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
    
    # DONE
    async def handle_uncow_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_moderation,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        select_sql = '''
            SELECT 1
            FROM active_cows
            WHERE channel_id = $1 AND discord_snowflake = $2 AND guild_id = $3
        '''
        delete_sql = '''
            DELETE FROM active_cows
            WHERE channel_id = $1 AND discord_snowflake = $2 AND guild_id = $3
        '''
        async with self.bot.db_pool.acquire() as conn:
            is_flagged = await conn.fetchval(select_sql, channel_obj.id, member_obj.id, message.guild.id)
            if not is_flagged:
                return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} has no active cow record in {channel_obj.mention}.')
            await conn.execute(delete_sql, channel_obj.id, member_obj.id, message.guild.id)
            await conn.execute(
                '''
                INSERT INTO moderation_logs (
                    action_type,
                    target_discord_snowflake,
                    executor_discord_snowflake,
                    guild_id,
                    channel_id,
                    reason
                )
                VALUES ($1, $2, $3, $4, $5, $6)
                ''',
                'uncow',
                member_obj.id,
                message.author.id,
                message.guild.id,
                channel_obj.id,
                'Uncowed a user'
            )
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} Cow Removed",
            description=(
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Action:** No longer going vegan"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
        
    # DONE
    async def handle_unflag_alias(self, alias, args, channel_obj, executor_role, existing_moderation, is_duration_modification, is_reason_modification, member_obj, message, state
    ):
        flags = await Flag.fetch_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        if not flags:
            return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently flagged in {channel_obj.mention}.')
        for flag in flags:
            if flag.expires_at is None and executor_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator'):
                return await state.end(warning='\U000026A0\U0000FE0F Only coordinators and above can undo permanent flags.')
            try:
                await channel_obj.set_permissions(member_obj, overwrite=None)
            except discord.Forbidden:
                return await state.end(error=f'\U0001F3C6 {member_obj.mention} was not successfully unflagged.')
        await Flag.delete_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} has been unflagged",
            description=(
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Unflagged By:** {message.author.mention}"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')

    # DONE
    async def handle_unmute_alias(self, alias, args, channel_obj, executor_role, existing_moderation, is_duration_modification, is_reason_modification, member_obj, message, state):
        voice_mutes = await VoiceMute.fetch_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        if not voice_mutes:
            return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently voice-muted in {channel_obj.mention}.')
        for voice_mute in voice_mutes:
            if voice_mute.expires_at is None and executor_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator'):
                return await state.end(warning='\U000026A0\U0000FE0F Only coordinators and above can undo permanent voice-mutes.')
            if member_obj.voice and member_obj.voice.channel:
                try:
                    await member_obj.edit(mute=False)
                except discord.Forbidden :
                    return await state.end(error=f'\U0001F3C6 {member_obj.mention}\'s voice-mute was not successfuly undone.')
        await VoiceMute.delete_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name}\'s voice-mute was undone.",
            description=(
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Voice-Mute Undone By:** {message.author.mention}"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')

    # DONE
    async def handle_unrole_alias(self, message: discord.Message, alias: Alias, args):
        highest_role = await is_owner_developer_administrator_coordinator_moderator(message)
        if highest_role == 'Everyone':
            return await message.reply(content='\U000026A0\U0000FE0F You are not permitted to unrole users.')
        member = args[0] if len(args) > 0 else None
        if not message.guild:
            return await message.reply(content='\U000026A0\U0000FE0F This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj:
            return await message.reply(content=f'\U000026A0\U0000FE0F Could not resolve a valid member from input: {member}.')
        if member_obj.id == message.guild.me.id:
            return await message.reply(content='\U000026A0\U0000FE0F You cannot unrole the bot.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel_snowflake)
        allowed, highest_role = await has_equal_or_higher_role(message, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)
        if not allowed:
            return await message.reply(content=f'\U000026A0\U0000FE0F You are not allowed to unrole this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        role_obj = message.guild.get_role(alias.role_snowflake)
        if not role_obj:
            return await message.reply(content=f"\U000026A0\U0000FE0F Could not resolve role with ID `{alias.role_snowflake}`.")
        if role_obj not in member_obj.roles:
            return await message.reply(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} does not have {role_obj.mention}.')
        try:
            await member_obj.remove_roles(role_obj)
        except discord.Forbidden:
            return await message.reply(content=f'\U000026A0\U0000FE0F {member_obj.mention} was not successfully unroled.')
        return await message.reply(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} had {role_obj.mention} removed.')
    
    # DONE
    async def handle_untextmute_alias(self, alias, args, channel_obj, executor_role, existing_moderation, is_duration_modification, is_reason_modification, member_obj, message, state):
        text_mutes = await TextMute.fetch_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        if not text_mutes:
            return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently text-muted in {channel_obj.mention}.')
        for text_mute in text_mutes:
            if text_mute.expires_at is None and executor_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator'):
                return await state.end(warning='\U000026A0\U0000FE0F Only coordinators and above can undo permanent text-mutes.')
            try:
                await channel_obj.set_permissions(member_obj, send_messages=None)
            except discord.Forbidden :
                return await state.end(error=f'\U0001F3C6 {member_obj.mention}\'s text-mute was not successfuly undone.')
        await VoiceMute.delete_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name}\'s text-mute was undone.",
            description=(
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Text-Mute Undone By:** {message.author.mention}"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
        
async def setup(bot: DiscordBot):
    cog = Aliases(bot)
    await bot.add_cog(cog)
