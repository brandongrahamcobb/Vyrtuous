
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
from datetime import datetime, timezone
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.role_service import RoleService
from vyrtuous.utils.alias import Alias
from vyrtuous.moderation_action.ban import Ban
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.properties.duration import DurationObject
from vyrtuous.utils.emojis import Emojis
from vyrtuous.moderation_action.flag import Flag
from vyrtuous.utils.invincibility import Invincibility
from vyrtuous.utils.history import History
from vyrtuous.moderation_action.text_mute import TextMute
from vyrtuous.moderation_action.vegan import Vegan
from vyrtuous.moderation_action.voice_mute import VoiceMute

class Aliases(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.alias_help = {
            'ban': [
                '**member** (Required): Tag a member or include their ID',
                '**duration** (Optional): (+|-)duration(m|h|d)\n0 = permanent / 24h = default\n`+` to append, `-` to delete, `=` to overwrite reason',
                '**reason** (Optional): Reason (required for 7 days or more)'
            ],
            'vegan': [
                '**member** (Required): Tag a member or include their ID'
            ],
            'carnist': [
                '**member** (Required): Tag a member or include their ID'
            ],
            'unban': [
                '**member** (Required): Tag a member or include their ID'
            ],
            'flag': [
                '**member** (Required): Tag a member or include their ID',
                '**reason** (Optional): Reason for flagging the user'
            ],
            'unflag': [
                '**member** (Required): Tag a member or include their ID'
            ],
            'voice_mute': [
                '**member** (Required): Tag a member or include their ID',
                '**duration** (Optional): (+|-)duration(m|h|d)\n0 = permanent / 24h = default\n`+` to append, `-` to delete, `=` to overwrite reason',
                '**reason** (Optional): Reason (required for 7 days or more)'
            ],
            'unvoice_mute': [
                '**member** (Required): Tag a member or include their ID'
            ],
            'text_mute': [
                '**member** (Required): Tag a member or include their ID',
                '**duration** (Required): (+|-)duration(m|h|d)\n0 = permanent / 24h = default\n`+` to append, `-` to delete, `=` to overwrite reason',
                '**reason** (Required): Reason (required for 7 days or more)'
            ],
            'untext_mute': [
                '**member** (Required): Tag a member or include their ID'
            ],
            'role': [
                '**member** (Required): Tag a member or include their ID',
                '**role** (Required): Role to assign'
            ],
            'unrole': [
                '**member** (Required): Tag a member or include their ID',
                '**role** (Required): Role to remove'
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
        self.role_service = RoleService()
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
        try:
            ban = None
            is_channel_scope = False
            is_modification = False
            override = False
            reason = 'No reason provided.'

            cap = await Cap.fetch_by_channel_guild_and_moderation_type(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, moderation_type='ban')
            if not hasattr(cap, 'duration'):
                cap_duration = DurationObject('24h').to_seconds()
            else:
                cap_duration = cap.duration    
                            
            if is_reason_modification and existing_guestroom_alias_event:
                is_modification = True
                duration = DurationObject.from_expires_in(existing_guestroom_alias_event.expires_in)
                modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
                match args[1]:
                    case '+':
                        reason = existing_guestroom_alias_event.reason + modified_reason
                    case '=' | '-':
                        reason = modified_reason
                await Alias.update_reason(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=Ban, updated_reason=reason)
            elif is_duration_modification and existing_guestroom_alias_event:
                is_modification = True
                duration = DurationObject(args[1] if len(args) > 1 else '8h')
                match duration.prefix:
                    case '+':
                        updated_expires_in = existing_guestroom_alias_event.expires_in + duration.to_timedelta()
                    case '=':
                        updated_expires_in = datetime.now(timezone.utc) + duration.to_timedelta()
                    case '-':
                        updated_expires_in = existing_guestroom_alias_event.expires_in - duration.to_timedelta()
                duration = DurationObject.from_expires_in(updated_expires_in)
                delta = updated_expires_in - datetime.now(timezone.utc)
                if delta.total_seconds() > cap_duration and executor_role not in ('System Owner', 'Guild Owner', 'Administrator', 'Coordinator'):
                    duration = DurationObject.from_seconds(cap_duration)
                    try:
                        return await state.end(warning=f'\u274C Cannot extend the ban beyond {duration} as a {executor_role} in {channel_obj.mention}.')
                    except:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')  
                await Alias.update_duration(channel_snowflake=channel_obj.id, expires_in=updated_expires_in, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=Ban)
            else:
                duration = DurationObject(args[1] if len(args) > 1 else '8h')
                delta = duration.expires_in - datetime.now(timezone.utc)
                if delta.total_seconds() < 0 and duration.number != 0:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to decrease the duration below the current time.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}') 
                ban = await Ban.fetch_by_channel_guild_and_member(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)
                if duration.number == 0:
                    override = True
                if not override:
                    if ban:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F An existing ban already exists for {member_obj.mention}. Modify the ban by putting a '+', '-' or '=' in front of the duration (ex. +8h) or a single '+', '-' or '=' and a reason to update the reason.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                if duration.to_seconds() > cap_duration and executor_role not in ('System Owner', 'Guild Owner', 'Administrator', 'Coordinator'):
                    duration = DurationObject.from_seconds(cap_duration)
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F Cannot set the ban beyond {duration} as a {executor_role} in {channel_obj.mention}.')
                    except:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')  
                reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'
                if ban and override:
                    await Alias.update_duration(channel_snowflake=channel_obj.id, expires_in=None, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=Ban)
                else:
                    ban = Ban(channel_snowflake=channel_obj.id, expires_in=duration.expires_in, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason)
                    await ban.create()
    
            try:
                await channel_obj.set_permissions(member_obj, view_channel=False, reason=reason)
            except discord.Forbidden as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                is_channel_scope = True
                try:
                    await member_obj.move_to(None, reason=reason)
                except discord.Forbidden as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')

            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {member_obj.display_name} has been Banned',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}\n'
                    f'**Expires:** {duration}\n'
                    f'**Reason:** {reason}'
                ),
                color=discord.Color.blue()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise
    
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
        try:
            duration = None
            is_channel_scope = False
            is_modification = False
            reason = 'No reason provided.'
    
            vegan = Vegan(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)
            await vegan.create()
    
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
            embed = discord.Embed(
                title=f'\U0001F525\U0001F525 {member_obj.display_name} is going Vegan!!!\U0001F525\U0001F525',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}\n'
                    f'**Celebrate!** Stick around and do some activism with us!'
                ),
                color=discord.Color.green()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise
    
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
        try:
            duration = None
            is_channel_scope = False
            is_modification = False
            reason = 'No reason provided.'
    
            if is_reason_modification and existing_guestroom_alias_event:
                is_modification = True
                modified_reason = ' '.join(args[1:]) if len(args) > 1 else ''
                match is_reason_modification:
                    case '+':
                        reason = existing_guestroom_alias_event.reason + modified_reason
                    case '=' | '-':
                        reason = modified_reason
                await Alias.update_reason(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=Flag, updated_reason=reason)
            else:
                if not is_modification:
                    flag = await Flag.fetch_by_channel_guild_and_member(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)   
                    if flag:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F An existing flag already exists for {member_obj.mention}. Modify the flag by putting a single '+', '-' or '=' and a reason to update the reason.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                reason = ' '.join(args[1:]) if len(args) > 1 else 'No reason provided.'
            
            flag = Flag(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason)
            await flag.create()
    
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
    
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {member_obj.display_name} Flagged',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}\n'
                    f'**Reason:** {reason}'
                ),
                color=discord.Color.red()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)

            bot = DiscordBot.get_instance()
            cog = bot.get_cog('EventListeners')
            cog.flags.append(flag)

            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise
   
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
        try:
            duration = None
            is_channel_scope = False
            is_modification = False
            reason = 'No reason provided.'
    
            try:
                role = await self.role_service.resolve_role(message, alias.role_snowflake)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Role `{alias.role_snowflake}` was not found.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            try:
                await member_obj.add_roles(role, reason='Added role')
            except discord.Forbidden as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
                
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
    
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {member_obj.display_name} Roled',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}\n'
                    f'**Role:** {role.mention}'
                ),
                color=discord.Color.blurple()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise

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
        try:
            is_channel_scope = False
            is_modification = False
            override = False
            reason = 'No reason provided.'
            text_mute = None
    
            cap = await Cap.fetch_by_channel_guild_and_moderation_type(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, moderation_type='text_mute')
            if not hasattr(cap, 'duration'):
                cap_duration = DurationObject('24h').to_seconds()
            else:
                cap_duration = cap.duration
                    
            if is_reason_modification and existing_guestroom_alias_event:
                is_modification = True
                modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
                match args[1]:
                    case '+':
                        reason = existing_guestroom_alias_event.reason + modified_reason
                    case '=' | '-':
                        reason = modified_reason
                await Alias.update_reason(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=TextMute, updated_reason=reason)
            elif is_duration_modification and existing_guestroom_alias_event:
                is_modification = True
                duration = DurationObject(args[1] if len(args) > 1 else '8h')
                match duration.prefix:
                    case '+':
                        updated_expires_in = existing_guestroom_alias_event.expires_in + duration.to_timedelta()
                    case '=':
                        updated_expires_in = datetime.now(timezone.utc) + duration.to_timedelta()
                    case '-':
                        updated_expires_in = existing_guestroom_alias_event.expires_in - duration.to_timedelta()
                duration = DurationObject.from_expires_in(updated_expires_in)
                delta = updated_expires_in - datetime.now(timezone.utc)
                if delta.total_seconds() > cap_duration and executor_role not in ('System Owner', 'Guild Owner', 'Administrator', 'Coordinator'):
                    duration = DurationObject.from_seconds(cap_duration)
                    try:
                        return await state.end(warning=f'\u274C Cannot extend the ban beyond {duration} as a {executor_role} in {channel_obj.mention}.')
                    except:
                        return await state.end(error=f'\u274C {str(e).capitalize()}') 
                await Alias.update_duration(channel_snowflake=channel_obj.id, expires_in=updated_expires_in, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=TextMute)
            else:
                duration = DurationObject(args[1] if len(args) > 1 else '8h')
                delta = duration.expires_in - datetime.now(timezone.utc)
                if delta.total_seconds() < 0 and duration.number != 0:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to decrease the duration below the current time.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}') 
                text_mute = await TextMute.fetch_by_channel_guild_and_member(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id)   
                if duration.number == 0:
                    override = True
                if not override:
                    if text_mute:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F An existing text-mute already exists for {member_obj.mention}. Modify the text-mute by putting a '+', '-' or '=' in front of the duration (ex. +8h) or a single '+', '-' or '=' and a reason to update the reason.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                if duration.to_seconds() > cap_duration and executor_role not in ('System Owner', 'Guild Owner', 'Administrator', 'Coordinator'):
                    duration = DurationObject.from_seconds(cap_duration)
                    try:
                        return await state.end(warning=f'\u274C Cannot set the ban beyond {duration} as a {executor_role} in {channel_obj.mention}.')
                    except:
                        return await state.end(error=f'\u274C {str(e).capitalize()}') 
                reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'
                if text_mute and override:
                    await Alias.update_duration(channel_snowflake=channel_obj.id, expires_in=None, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=TextMute)
                else:
                    text_mute = TextMute(channel_snowflake=channel_obj.id, expires_in=duration.expires_in, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason)
                    await text_mute.create()
    
            try:
                await channel_obj.set_permissions(member_obj, send_messages=False, add_reactions=False, reason=reason)
            except discord.Forbidden as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
                
    
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
    
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {member_obj.display_name} Text Muted',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}\n'
                    f'**Expires:** {duration}\n'
                    f'**Reason:** {reason}'
                ),
                color=discord.Color.green()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise
        
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
        try:
            is_channel_scope = False
            is_modification = False
            override = False
            target = 'user'
            reason = 'No reason provided.'
            voice_mute = None
    
            cap = await Cap.fetch_by_channel_guild_and_moderation_type(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, moderation_type='voice_mute')
            if not hasattr(cap, 'duration'):
                cap_duration = DurationObject('24h').to_seconds()
            else:
                cap_duration = cap.duration
    
            if is_reason_modification and existing_guestroom_alias_event:
                is_modification = True
                duration = DurationObject.from_expires_in(existing_guestroom_alias_event.expires_in)
                modified_reason = ' '.join(args[2:]) if len(args) > 2 else ''
                match args[1]:
                    case '+':
                        reason = existing_guestroom_alias_event.reason + modified_reason
                    case '=' | '-':
                        reason = modified_reason
                await Alias.update_reason(alias_type=alias.alias_type, channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=VoiceMute, updated_reason=reason)
            elif is_duration_modification and existing_guestroom_alias_event:
                is_modification = True
                duration = DurationObject(args[1] if len(args) > 1 else '8h')
                match duration.prefix:
                    case '+':
                        updated_expires_in = existing_guestroom_alias_event.expires_in + duration.to_timedelta()
                    case '=':
                        updated_expires_in = datetime.now(timezone.utc) + duration.to_timedelta()
                    case '-':
                        updated_expires_in = existing_guestroom_alias_event.expires_in - duration.to_timedelta()
                duration = DurationObject.from_expires_in(updated_expires_in)
                delta = updated_expires_in - datetime.now(timezone.utc)
                if delta.total_seconds() > cap_duration and executor_role not in ('System Owner', 'Guild Owner', 'Administrator', 'Coordinator'):
                    duration = DurationObject.from_seconds(cap_duration)
                    try:
                        return await state.end(warning=f'\u274C Cannot extend the ban beyond {duration} as a {executor_role} in {channel_obj.mention}.')
                    except:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                await Alias.update_duration(channel_snowflake=channel_obj.id, expires_in=updated_expires_in, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=VoiceMute)
            else:
                duration = DurationObject(args[1] if len(args) > 1 else '8h')
                delta = duration.expires_in - datetime.now(timezone.utc)
                if delta.total_seconds() < 0 and duration.number != 0:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to decrease the duration below the current time.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}') 
                voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, target=target)    
                if duration.number == 0:
                    override = True
                if not override:
                    if voice_mute:
                        try:
                            return await state.end(warning=f"\U000026A0\U0000FE0F An existing voice-mute already exists for {member_obj.mention}. Modify the voice-mute by putting a '+', '-' or '=' in front of the duration (ex. +8h) or a single '+', '-' or '=' and a reason to update the reason.")
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                if duration.to_seconds() > cap_duration and executor_role not in ('System Owner', 'Guild Owner', 'Administrator', 'Coordinator'):
                    duration = DurationObject.from_seconds(cap_duration)
                    try:
                        return await state.end(warning=f'\u274C Cannot set the ban beyond {duration} as a {executor_role} in {channel_obj.mention}.')
                    except:
                        return await state.end(error=f'\u274C {str(e).capitalize()}') 
                reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'
                if voice_mute and override:
                    await Alias.update_duration(channel_snowflake=channel_obj.id, expires_in=None, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, moderation_type=VoiceMute)
                else:
                    voice_mute = VoiceMute(channel_snowflake=channel_obj.id, expires_in=duration.expires_in, guild_snowflake=message.guild.id, member_snowflake=member_obj.id, reason=reason, target=target)
                    await voice_mute.create()


            if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                is_channel_scope = True
                try:
                    await member_obj.edit(mute=True, reason=reason)
                except discord.Forbidden as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
    
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
    
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {member_obj.display_name} has been Voice Muted',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}\n'
                    f'**Expires:** {duration}\n'
                    f'**Reason:** {reason}'
                ),
                color=discord.Color.green()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise
        
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
        try:
            duration = None
            is_channel_scope = False
            is_modification = True
            reason = 'No reason provided.'
    
            ban = await Ban.fetch_by_channel_guild_and_member(
                channel_snowflake=channel_obj.id,
                guild_snowflake=message.guild.id,
                member_snowflake=member_obj.id
            )
            if not ban:
                return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently banned in {channel_obj.mention}.')
            if ban.expires_in is None and executor_role not in ('System Owner', 'Guild Owner', 'Administrator', 'Coordinator'):
                try:
                    return await state.end(warning='\U000026A0\U0000FE0F Only coordinators and above can undo permanent bans.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            await Ban.delete_by_channel_guild_and_member(
                channel_snowflake=channel_obj.id,
                guild_snowflake=message.guild.id,
                member_snowflake=member_obj.id
            )
    
            try:
                await channel_obj.set_permissions(member_obj, overwrite=None)
            except discord.Forbidden as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
                
    
            if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
                is_channel_scope = True
                try:
                    await member_obj.move_to(None, reason='Unbanned')
                except discord.Forbidden as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
    
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
    
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {member_obj.display_name} has been Unbanned',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}'
                ),
                color=discord.Color.yellow()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise
    
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
        try:
            duration = None
            is_channel_scope = False
            is_modification = True
            reason = 'No reason provided.'
    
            await Vegan.delete_by_channel_guild_and_member(channel_snowflake=channel_obj.id, member_snowflake=member_obj.id, guild_snowflake=message.guild.id)
            
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
    
            embed = discord.Embed(
                title=f'\U0001F44E\U0001F44E {member_obj.display_name} is a Carnist \U0001F44E\U0001F44E',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}'
                ),
                color=discord.Color.red()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise
        
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
        try:
            duration = None
            is_channel_scope = False
            is_modification = True
            reason = 'No reason provided.'
    
            flag = await Flag.fetch_by_channel_guild_and_member(
                channel_snowflake=channel_obj.id,
                guild_snowflake=message.guild.id,
                member_snowflake=member_obj.id
            )
            if not flag:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently flagged in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')

            bot = DiscordBot.get_instance()
            cog = bot.get_cog('EventListeners')
            for flag in cog.flags:
                if flag.channel_snowflake == channel_obj.id:
                    cog.flags.remove(flag)
                    break

            await Flag.delete_by_channel_guild_and_member(
                channel_snowflake=channel_obj.id,
                guild_snowflake=message.guild.id,
                member_snowflake=member_obj.id
            )   
    
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
                                            
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {member_obj.display_name} has been Unflagged',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}'
                ),
                color=discord.Color.yellow()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)

            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise

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
        try:
            duration = None
            is_channel_scope = False
            is_modification = True
            reason = 'No reason provided.'
    
            voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(
                channel_snowflake=channel_obj.id,
                guild_snowflake=message.guild.id,
                member_snowflake=member_obj.id,
                target='user'
            )
            if not voice_mute:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently voice-muted in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            if voice_mute.expires_in is None and executor_role not in ('System Owner', 'Guild Owner', 'Administrator', 'Coordinator'):
                try:
                    return await state.end(warning='\U000026A0\U0000FE0F Only coordinators and above can undo permanent voice-mutes.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            
            await VoiceMute.delete_by_channel_guild_member_and_target(
                channel_snowflake=channel_obj.id,
                guild_snowflake=message.guild.id,
                member_snowflake=member_obj.id,
                target='user'
            )
                    
            if member_obj.voice and member_obj.voice.channel:
                try:
                    is_channel_scope = True
                    await member_obj.edit(mute=False)
                except discord.Forbidden as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
    
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
    
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {member_obj.display_name} has been Unmuted',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}'
                ),
                color=discord.Color.yellow()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise

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
        try:
            duration = None
            is_channel_scope = False
            is_modification = True
            reason = 'No reason provided.'
    
            try:
                role = await self.role_service.resolve_role(message, alias.role_snowflake)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Role `{alias.role_snowflake}` was not found.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            if role not in member_obj.roles:
                try:
                    return await state.end(warning=f'{self.emoji.get_random_emoji()} {member_obj.mention} does not have {role.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            try:
                await member_obj.remove_roles(role)
            except discord.Forbidden as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
                
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
    
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {member_obj.display_name} has been Unroled',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}\n'
                    f'**Role:** {role.mention}'
                ),
                color=discord.Color.yellow()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise
    
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
        try:
            duration = None
            is_channel_scope = False
            is_modification = True
            reason = 'No reason provided.'
    
            text_mute = await TextMute.fetch_by_channel_guild_and_member(
                channel_snowflake=channel_obj.id,
                guild_snowflake=message.guild.id,
                member_snowflake=member_obj.id
            )
            if not text_mute:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not currently text-muted in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            if text_mute.expires_in is None and executor_role not in ('System Owner', 'Guild Owner', 'Administrator', 'Coordinator'):
                try:
                    return await state.end(warning='\U000026A0\U0000FE0F Only coordinators and above can undo permanent text-mutes.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
                
            await TextMute.delete_by_channel_guild_and_member(
                channel_snowflake=channel_obj.id,
                guild_snowflake=message.guild.id,
                member_snowflake=member_obj.id
            )
                        
            try:
                await channel_obj.set_permissions(member_obj, send_messages=None)
            except discord.Forbidden as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
                
            await History.send_entry(alias, channel_obj, duration, executor_role, is_channel_scope, is_modification, member_obj, message, reason)
    
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {member_obj.display_name} has been Unmuted',
                description=(
                    f'**By:** {message.author.mention}\n'
                    f'**User:** {member_obj.mention}\n'
                    f'**Channel:** {channel_obj.mention}'
                ),
                color=discord.Color.yellow()
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        except Exception as e:
            raise
        
async def setup(bot: DiscordBot):
    cog = Aliases(bot)
    await bot.add_cog(cog)
