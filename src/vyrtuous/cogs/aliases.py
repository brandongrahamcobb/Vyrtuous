
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
from vyrtuous.utils.duration import Duration, DurationObject
from vyrtuous.utils.flag import Flag
from vyrtuous.utils.reason import Reason
from vyrtuous.utils.statistics import Statistics
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.state import State
from vyrtuous.utils.invincibility import Invincibility
from vyrtuous.utils.vegan import Vegan
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
            'vegan': [
                "**member** (Optional): Tag a member or include their snowflake ID"
            ],
            'carnist': [
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
            'voice_mute': [
                "**member** (Optional): Tag a member or include their snowflake ID",
                "**duration** (Optional): (+|-)duration(m|h|d)\n0 = permanent / 24h = default\n`+` to append, `-` to delete, `=` to overwrite reason",
                "**reason** (Optional): Optional reason (required for 7 days or more)"
            ],
            'unvoice_mute': [
                "**member** (Optional): Tag a member or include their snowflake ID"
            ],
            'text_mute': [
                "**member** (Optional): Tag a member or include their snowflake ID",
                "**duration** (Optional): (+|-)duration(m|h|d)\n0 = permanent / 24h = default\n`+` to append, `-` to delete, `=` to overwrite reason",
                "**reason** (Optional): Optional reason (required for 7 days or more)"
            ],
            'untext_mute': [
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
            'vegan': 'Verifies a user as going vegan.',
            'carnist': 'Unverifies a user as going vegan.',
            'unban': 'Unbans a user from the server.',
            'flag': 'Flags a user for moderation review.',
            'unflag': 'Removes a flag from a user.',
            'voice_mute': 'Mutes a user in voice channels.',
            'unvoice_mute': 'Unmutes a user in voice channels.',
            'text_mute': 'Mutes a user in text channels.',
            'untext_mute': 'Unmutes a user in text channels.',
            'role': 'Assigns a role to a user.',
            'unrole': 'Removes a role from a user.'
        }
        self.alias_type_to_permission_level = {
            'ban': 'Moderator',
            'vegan': 'Moderator',
            'carnist': 'Moderator',
            'unban': 'Moderator',
            'voice_mute': 'Moderator',
            'unvoice_mute': 'Moderator',
            'text_mute': 'Moderator',
            'untext_mute': 'Moderator',
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
    
    async def handle_ban_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        is_channel_scope = False
        is_modification = False
        if is_reason_modification and existing_guestroom_alias_event:
            is_modification = True
            duration = DurationObject.from_expires_at(existing_guestroom_alias_event.expires_at)
            modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
            match is_reason_modification:
                case '+':
                    updated_reason = existing_guestroom_alias_event.reason + modified_reason
                case '=' | '-':
                    updated_reason = modified_reason
            await Alias.update_reason(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, updated_reason=updated_reason)
        elif is_duration_modification and existing_guestroom_alias_event:
            is_modification = True
            duration = DurationObject(args[1] if len(args) > 1 else '24h')
            match duration.prefix:
                case '+':
                    updated_expires_at = existing_guestroom_alias_event.expires_at + duration.to_timedelta()
                case '=':
                    updated_expires_at = datetime.now(timezone.utc) + duration.to_timedelta()
                case '-':
                    updated_expires_at = existing_guestroom_alias_event.expires_at - duration.to_timedelta()
            duration = DurationObject.from_expires_at(updated_expires_at)
            await Alias.update_duration(channel_snowflake=channel_obj.id, expires_at=updated_expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=Ban)
        else:
            duration = DurationObject(args[1] if len(args) > 1 else '24h')
            reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'

        try:
            await channel_obj.set_permissions(member_obj, view_channel=False, reason=reason)
        except discord.Forbidden:
            return await state.end(error=f'\U0001F3C6 {e}')
        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            is_channel_scope = True
            try:
                await member_obj.move_to(None, reason=reason)
            except discord.Forbidden:
                try:
                    return await state.end(error=f'\U0001F3C6 Unsuccessfully banned {member_obj.mention}.')
                except:
                    return await state.end(error=f'\U0001F3C6 {e}')
                
        ban = Ban(channel_snowflake=channel_obj.id, expires_at=duration.expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason)
        await ban.create()

        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Banned",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Duration:** {duration}\n"
                f"**Reason:** {reason}"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
    
    # DONE
    async def handle_vegan_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message, 
        state
    ):
        duration = None
        is_channel_scope = False
        is_modification = False
        reason = None

        vegan = Vegan(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)
        await vegan.create()

        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Vegan",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"Celebrate!** \U0001F525 {member_obj.mention} is going vegan!!! \U0001F525"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
    
    # DONE
    async def handle_flag_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        duration = None
        is_channel_scope = False
        is_modification = False

        if is_reason_modification and existing_guestroom_alias_event:
            is_modification = True
            modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
            match is_reason_modification:
                case '+':
                    reason = existing_guestroom_alias_event.reason + modified_reason
                case '=' | '-':
                    reason = modified_reason
            await Alias.update_reason(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, updated_reason=reason)
        else:
            reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'
        
        flag = Flag(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason)
        await flag.create()

        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)

        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Flagged",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Duration:** {duration}\n"
                f"**Reason:** {reason}"
            ),
            color=discord.Color.green()
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
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        duration = None
        is_channel_scope = False
        is_modification = False
        reason = None

        try:
            role = await self.role_service.resolve_role(alias.role_snowflake)
        except Exception as e:
            try:
                return await state.end(warning=f'\U0001F3C6 Role `{alias.rle_snowflake}` was not found.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        try:
            await member_obj.add_roles(role, reason='Added role')
        except discord.Forbidden:
            try:
                return await state.end(error=f'\U0001F3C6 {member_obj.mention} was not successfully roled.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
            
        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)

        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Roled",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')

    # DONE
    async def handle_text_mute_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        is_channel_scope = False
        is_modification = False

        if is_reason_modification and existing_guestroom_alias_event:
            is_modification = True
            modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
            match is_reason_modification:
                case '+':
                    updated_reason = existing_guestroom_alias_event.reason + modified_reason
                case '=' | '-':
                    updated_reason = modified_reason
            await Alias.update_reason(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, updated_reason=updated_reason)
        elif is_duration_modification:
            is_modification = True
            if existing_guestroom_alias_event:
                duration = DurationObject(args[1] if len(args) > 1 else '24h')
                match duration.prefix:
                    case '+':
                        updated_expires_at = existing_guestroom_alias_event.expires_at + duration.to_timedelta()
                    case '=':
                        updated_expires_at = datetime.now(timezone.utc) + duration.to_timedelta()
                    case '-':
                        updated_expires_at = existing_guestroom_alias_event.expires_at - duration.to_timedelta()
                duration = DurationObject.from_expires_at(updated_expires_at)
                await Alias.update_duration(channel_snowflake=channel_obj.id, expires_at=updated_expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=TextMute)
        else:
            duration = DurationObject(args[1] if len(args) > 1 else '24h')
            reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'

        try:
            await channel_obj.set_permissions(member_obj, send_messages=False, add_reactions=False, reason=reason)
        except discord.Forbidden:
            try:
                return await state.end(error=f'\U0001F3C6 {member_obj.mention} was not successfully text-muted.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
            
        text_mute = TextMute(channel_snowflake=channel_obj.id, expires_at=duration.expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason)
        await text_mute.create()

        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)

        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Muted",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Duration:** {duration}\n"
                f"**Reason:** {reason}"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
        
    # DONE
    async def handle_voice_mute_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message, 
        state
    ):
        is_channel_scope = False
        is_modification = False
        if is_reason_modification and existing_guestroom_alias_event:
            is_modification = True
            duration = DurationObject.from_expires_at(existing_guestroom_alias_event.expires_at)
            modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
            match is_reason_modification:
                case '+':
                    reason = existing_guestroom_alias_event.reason + modified_reason
                case '=' | '-':
                    reason = modified_reason
            await Alias.update_reason(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, updated_reason=reason)
        elif is_duration_modification:
            is_modification = True
            if existing_guestroom_alias_event:
                duration = DurationObject(args[1] if len(args) > 1 else '24h')
                match duration.prefix:
                    case '+':
                        expires_at = existing_guestroom_alias_event.expires_at + duration.to_timedelta()
                    case '=':
                        expires_at = datetime.now(timezone.utc) + duration.to_timedelta()
                    case '-':
                        expires_at = existing_guestroom_alias_event.expires_at - duration.to_timedelta()
                duration = DurationObject.from_expires_at(expires_at)
                await Alias.update_duration(channel_snowflake=channel_obj.id, expires_at=expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=VoiceMute)
        else:
            duration = DurationObject(args[1] if len(args) > 1 else '24h')
            reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'

        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            is_channel_scope = True
            try:
                await member_obj.move_to(None, reason=reason)
            except discord.Forbidden:
                return await state.end(error=f'\U0001F3C6 {e}')

        voice_mute = VoiceMute(channel_snowflake=channel_obj.id, expires_at=duration.expires_at, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason, target="user")
        await voice_mute.create()

        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)

        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Muted",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Duration:** {duration}\n"
                f"**Reason:** {reason}"
            ),
            color=discord.Color.green()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
        
    # DONE
    async def handle_unban_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        duration = None
        is_channel_scope = False
        is_modification = True
        reason = None

        ban = await Ban.fetch_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        if not ban:
            return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently banned in {channel_obj.mention}.')
        if ban.expires_at is None and executor_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator'):
            try:
                return await state.end(warning='\U000026A0\U0000FE0F Only coordinators and above can undo permanent bans.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        await Ban.delete_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )

        try:
            await channel_obj.set_permissions(member_obj, overwrite=None)
        except discord.Forbidden:
            try:
                return await state.end(error=f'\U0001F3C6 {member_obj.mention} was not successfully unbanned.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
            

        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            is_channel_scope = True
            try:
                await member_obj.move_to(None, reason="Unbanned")
            except discord.Forbidden:
                try:
                    return await state.end(error=f'\U0001F3C6 Could not move {member_obj.mention} from the voice channel.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}')

        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)

        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Unbanned",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}",
            ),
            color=discord.Color.yellow()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
    
    # DONE
    async def handle_carnist_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        duration = None
        is_channel_scope = False
        is_modification = True
        reason = None

        await Vegan.delete_by_channel_guild_and_member(channel_snowflake=channel_obj.id, member_snowflake=member_obj.id, guild_snowflake=message.guild.id)
        
        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)

        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Carnist",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}"
            ),
            color=discord.Color.yellow()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
        
    # DONE
    async def handle_unflag_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        duration = None
        is_channel_scope = False
        is_modification = True
        reason = None

        flag = await Flag.fetch_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        if not flag:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently flagged in {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        await Flag.delete_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )   

        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
                                        
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Unflagged",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}"
            ),
            color=discord.Color.red()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')

    # DONE
    async def handle_unmute_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        duration = None
        is_channel_scope = False
        is_modification = True
        reason = False

        voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
            target="user"
        )
        if not voice_mute:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently voice-muted in {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        if voice_mute.expires_at is None and executor_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator'):
            try:
                return await state.end(warning='\U000026A0\U0000FE0F Only coordinators and above can undo permanent voice-mutes.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        
        if member_obj.voice and member_obj.voice.channel:
            try:
                is_channel_scope = True
                await member_obj.edit(mute=False)
            except discord.Forbidden:
                try:
                    return await state.end(error=f'\U0001F3C6 {member_obj.mention}\'s voice-mute was not successfuly undone.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}')
        await VoiceMute.delete_by_channel_guild_member_and_target(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
            target="user"
        )

        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)

        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Unmuted",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}"
            ),
            color=discord.Color.yellow()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')

    # DONE
    async def handle_unrole_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        duration = None
        is_channel_scope = False
        is_modification = True
        reason = None

        try:
            role = await self.role_service.resolve_role(alias.role_snowflake)
        except Exception as e:
            try:
                return await state.end(warning=f'\U0001F3C6 Role `{alias.role_snowflake}` was not found.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        if role not in member_obj.roles:
            try:
                return await state.end(warning=f'{self.emoji.get_random_emoji()} {member_obj.mention} does not have {role.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        try:
            await member_obj.remove_roles(role)
        except discord.Forbidden:
            try:
                return await state.end(error=f'\U000026A0\U0000FE0F {member_obj.mention} was not successfully unroled.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
            
        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)

        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Unroled",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
    
    # DONE
    async def handle_untextmute_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state
    ):
        duration = None
        is_channel_scope = False
        is_modification = True
        reason = None

        text_mute = await TextMute.fetch_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
        if not text_mute:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently text-muted in {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        if text_mute.expires_at is None and executor_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator'):
            try:
                return await state.end(warning='\U000026A0\U0000FE0F Only coordinators and above can undo permanent text-mutes.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        await TextMute.delete_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id
        )
                    
        try:
            await channel_obj.set_permissions(member_obj, send_messages=None)
        except discord.Forbidden :
            try:
                return await state.end(error=f'\U0001F3C6 {member_obj.mention}\'s text-mute was not successfuly undone.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
            
        await Statistics.send_statistic(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)

        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} Unmuted",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
            ),
            color=discord.Color.yellow()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
        
async def setup(bot: DiscordBot):
    cog = Aliases(bot)
    await bot.add_cog(cog)
