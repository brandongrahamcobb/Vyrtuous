''' admin_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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
from typing import Optional

from discord import app_commands

from vyrtuous.inc.helpers import *
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.enhanced_members.administrator import AdministratorRole
from vyrtuous.enhanced_members.coordinator import Coordinator
from vyrtuous.enhanced_members.moderator import Moderator
from vyrtuous.moderation_action.ban import Ban
from vyrtuous.moderation_action.flag import Flag
from vyrtuous.moderation_action.server_mute import ServerMute
from vyrtuous.moderation_action.text_mute import TextMute
from vyrtuous.moderation_action.vegan import Vegan
from vyrtuous.moderation_action.voice_mute import VoiceMute
from vyrtuous.rooms.stage import Stage
from vyrtuous.rooms.temporary_room import TemporaryRoom
from vyrtuous.rooms.video_room import VideoRoom
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.check_service import *
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.role_service import RoleService
from vyrtuous.service.state_service import State
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.cancel_confirm import VerifyView
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.properties.duration import AppDuration, Duration, DurationObject
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.properties.moderation_type import AppModerationType, ModerationType
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.setup_logging import logger
from vyrtuous.utils.permission import TARGET_PERMISSIONS
from vyrtuous.utils.properties.snowflake import *
from vyrtuous.utils.history import History

class AdminCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.emoji = Emojis()
        self.message_service = MessageService(self.bot, self.bot.db_pool)
        self.member_service = MemberService()
        self.role_service = RoleService()
    
    # DONE
    @app_commands.command(name='alias', description='Alias creation.')
    @administrator_predicator()
    @app_commands.describe(
        alias_name='Alias/Pseudonym',
        moderation_type='One of: vegan, carnist, vmute, unvmute,' \
                        'ban, unban, flag, unflag, tmute, untmute, role, unrole',
        channel='Tag a channel or include its ID',
        role='Role ID (only for role/unrole)'
    )
    async def create_alias_app_command(
        self,
        interaction: discord.Interaction,
        moderation_type: AppModerationType,
        alias_name: Optional[str],
        channel: AppChannelSnowflake,
        role: AppRoleSnowflake = None
    ):
        state = State(interaction)
        channel_obj, role_obj = None, None
        try:
            channel_obj = \
                await self.channel_service.search(interaction, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        try:
            role_obj = await self.role_service.search(interaction, role)
            role_snowflake = role_obj.id
        except Exception as e:
            role_snowflake = None
        alias = Alias(alias_name=alias_name, alias_type=moderation_type,
            channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id,
            role_snowflake=role_snowflake
        )
        await alias.create()
        if role_snowflake:
            msg = f'Alias `{alias_name}` of type `{moderation_type}` ' \
                f'created successfully for channel {channel_obj.mention} ' \
                f'and role {role_obj.mention}.'
        else:
            msg = f'Alias `{alias_name}`  of type `{moderation_type}` ' \
                f'created successfully for channel {channel_obj.mention}.'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(
        name='alias',
        help='Set an alias for a vegan, carnist, mute, unmute, ban, ' \
             'unban, flag, unflag, tmute, untmute, role, or unrole action.'
    )
    @administrator_predicator()
    async def create_alias_text_command(
        self,
        ctx: commands.Context,
        moderation_type: ModerationType = commands.parameter(
            default=None,
            description='One of: `vegan`, `carnist`, `vmute`, ' \
                '`unvmute`, `ban`, `unban`, `flag`, ' \
                '`unflag`, `tmute`, `untmute`, `role`, `unrole`'
        ),
        alias_name: str = commands.parameter(
            default=None,
            description='Alias/Pseudonym'
        ),
        channel: ChannelSnowflake = commands.parameter(
            default=None,
            description='Tag a channel or include its ID'
        ), *,
        role: Optional[RoleSnowflake] = commands.parameter(
            default=None,
            description='Role ID (only for role/unrole)'
        )
    ):
        state = State(ctx)
        channel_obj, role_obj = None, None
        try:
            channel_obj = \
                await self.channel_service.search(ctx, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f' {str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        try:
            role_obj = await self.role_service.search(ctx, role)
            role_snowflake = role_obj.id
        except Exception as e:
            role_snowflake = None
        alias = Alias(alias_name=alias_name,  alias_type=moderation_type,
            channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id,
            role_snowflake=role_snowflake
        )
        await alias.create()
        if role_snowflake:
            msg = f'Alias `{alias_name}` of type `{moderation_type}` ' \
                f'created successfully for channel {channel_obj.mention} ' \
                f'and role {role_obj.mention}.'
        else:
            msg = f'Alias `{alias_name}`  of type `{moderation_type}` ' \
                f'created successfully for channel {channel_obj.mention}.'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    @app_commands.command(name='aroles', description='Administrator roles.')
    @app_commands.describe(scope="Specify one of: 'all', server ID or empty.")
    @administrator_predicator()
    async def list_administrator_roles_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        chunk_size, field_count, pages = 7, 0, []
        guild_obj = None
        is_at_home = False
        guild_dictionary, skipped_guilds, skipped_roles = {}, set(), {}
        title = f'{self.emoji.get_random_emoji()} Administrator Roles'

        highest_role = await permission_check(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('System Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F' \
                        'You are not authorized to list ' \
                        'administrator roles across all servers.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            administrator_roles = \
                await AdministratorRole.fetch_all_guilds_and_roles()
        elif scope:
            if highest_role not in (
                'System Owner',
                'Developer',
                'Guild Owner',
                'Administrator'
            ):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        'You are not authorized ' \
                        'to list all administrator roles in a specific server.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            guild_obj = self.bot.get_guild(int(scope))
            if not guild_obj:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        'Scope must be one of: ' \
                        f"'all', channel ID/mention, server ID or empty. " \
                        f'Received: {scope}.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            administrator_roles = \
                await AdministratorRole.fetch_by_guild(
                    guild_snowflake=guild_obj.id
                )
        else:
            administrator_roles = \
                await AdministratorRole.fetch_by_guild(
                    guild_snowflake=interaction.guild.id
                )
            guild_obj = interaction.guild
        
        if not administrator_roles:
            try:
                if scope:
                    msg = f'No administrator roles setup for scope: {scope}.'
                else:
                    msg = f'No administrator roles setup in {guild_obj.name}.'
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        for administrator_role in administrator_roles:
            guild_dictionary.setdefault(administrator_role.guild_snowflake, {})
            guild_dictionary[administrator_role.guild_snowflake].setdefault(
                'role_snowflake',
                []
            ).append(administrator_role.role_snowflake)
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = \
                dict(sorted(guild_dictionary[guild_snowflake].items()))
        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guilds.add(guild_snowflake)
                continue
            embed = discord.Embed(
                title=title,
                description=guild.name,
                color=discord.Color.blue()
            )
            for role_snowflake in guild_data.get('role_snowflake', []):
                role = guild.get_role(role_snowflake)
                if not role:
                    skipped_roles.setdefault(
                        guild_snowflake, 
                        []
                    ).append(role_snowflake)
                    continue
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(
                        title=f'{title} continued...',
                        description=guild.name,
                        color=discord.Color.blue()
                    )
                    field_count = 0
                embed.add_field(name=role.name, value=role.mention, inline=False)
                field_count += 1
            pages.append(embed)
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guilds:
                embed = discord.Embed(
                    title='Skipped Servers',
                    description='\u200b',
                    color=discord.Color.blue()
                )
                lines = []
                for guild_snowflake in skipped_guilds:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            title='Skipped Servers continued...',
                            color=discord.Color.red()
                        )
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_roles:
                for guild_snowflake, role_snowflakes in skipped_roles.items():
                    embed = discord.Embed(
                        title=f'Skipped Roles in {guild_snowflake}',
                        description='\u200b',
                        color=discord.Color.red()
                    )
                    field_count = 0
                    for role_snowflake in role_snowflakes:
                        if field_count >= chunk_size:
                            pages.append(embed)
                            embed = discord.Embed(
                                title=f'Skipped Roles in ' \
                                    f'{guild_snowflake} continued...',
                                description='\u200b',
                                color=discord.Color.red()
                            )
                            field_count = 0
                        embed.add_field(
                            name=str(role_snowflake),
                            value='\u200b',
                            inline=False
                        )
                        field_count += 1
                    pages.append(embed)
            
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        'Embed size is too large. Limit the scope.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    'No administrator roles found.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    @commands.command(name='aroles', help='Administrator roles.')
    @administrator_predicator()
    async def list_administrator_roles_text_command(
        self,
        ctx: commands.Context,
        *,
        scope: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: 'all', " \
                'channel ID/mention, server ID or empty.'
        )
    ):
        state = State(ctx)
        chunk_size, field_count, pages = 7, 0, []
        is_at_home = False
        guild_obj = None
        guild_dictionary, skipped_guilds, skipped_roles = {}, set(), {}
        title = f'{self.emoji.get_random_emoji()} Administrator Roles'

        highest_role = await permission_check(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('System Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        'You are not authorized to list ' \
                        'administrator roles across all servers.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            administrator_roles = await AdministratorRole.fetch_all()
        elif scope:
            if highest_role not in (
                'System Owner',
                'Developer',
                'Guild Owner',
                'Administrator'
            ):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        'You are not authorized to list ' \
                        'all administrator roles in a specific server.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            guild_obj = self.bot.get_guild(int(scope))
            if not guild_obj:
                try:
                    return await state.end(warning=f"\U000026A0\U0000FE0F " \
                        f"Scope must be one of: 'all', channel ID/mention, " \
                        f"server ID or empty. Received: {scope}."
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            administrator_roles = \
                await AdministratorRole.fetch_by_guild(guild_snowflake=guild_obj.id)
        else:
            administrator_roles = \
                await AdministratorRole.fetch_by_guild(guild_snowflake=ctx.guild.id)
            guild_obj = ctx.guild
        
        if not administrator_roles:
            try:
                if scope:
                    msg = f'No administrator roles setup for scope: {scope}.'
                else:
                    msg = f'No administrator roles setup in {guild_obj.name}.'
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        for administrator_role in administrator_roles:
            guild_dictionary.setdefault(administrator_role.guild_snowflake, {})
            guild_dictionary[administrator_role.guild_snowflake].setdefault(
                'role_snowflake',
                []).append(administrator_role.role_snowflake)
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = \
                dict(sorted(guild_dictionary[guild_snowflake].items()))
        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guilds.add(guild_snowflake)
                continue
            embed = discord.Embed(
                title=title,
                description=guild.name,
                color=discord.Color.blue()
            )
            for role_snowflake in guild_data.get('role_snowflake', []):
                role = guild.get_role(role_snowflake)
                if not role:
                    skipped_roles.setdefault(
                        guild_snowflake,
                        []
                    ).append(role_snowflake)
                    continue
                if field_count >= chunk_size:
                    pages.append(embed)
                    embed = discord.Embed(
                        title=f'{title} continued...',
                        description=guild.name,
                        color=discord.Color.blue()
                    )
                    field_count = 0
                embed.add_field(
                    name=role.name,
                    value=role.mention,
                    inline=False
                )
                field_count += 1
            pages.append(embed)
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guilds:
                embed = discord.Embed(
                    title='Skipped Servers',
                    description='\u200b',
                    color=discord.Color.blue()
                )
                lines = []
                for guild_snowflake in skipped_guilds:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            title='Skipped Servers continued...',
                            color=discord.Color.red()
                )
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_roles:
                for guild_snowflake, role_snowflakes in skipped_roles.items():
                    embed = discord.Embed(
                        title=f'Skipped Roles in {guild_snowflake}',
                        description='\u200b',
                        color=discord.Color.red()
                    )
                    field_count = 0
                    for role_snowflake in role_snowflakes:
                        if field_count >= chunk_size:
                            pages.append(embed)
                            embed = discord.Embed(
                                title=f'Skipped Roles in '\
                                    f'{guild_snowflake} continued...',
                                description='\u200b',
                                color=discord.Color.red
                            )
                            field_count = 0
                        embed.add_field(
                            name=str(role_snowflake),
                            value='\u200b',
                            inline=False
                        )
                        field_count += 1
                    pages.append(embed)
            
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        'Embed size is too large. Limit the scope.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    'No administrator roles found.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    @app_commands.command(name='cap', description='Cap alias duration for mods.')
    @administrator_predicator()
    @app_commands.describe(
        channel='Tag a channel or include its ID',
        moderation_type='One of: `mute`, `ban`, `tmute`',
        hours='(+|-)duration(m|h|d), 0=permanent, default=24h'
    )
    async def cap_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        moderation_type: AppModerationType,
        hours: int
    ):
        state = State(interaction)
        channel_obj = None
        seconds = int(hours) * 3600
        try:
            channel_obj = \
                await self.channel_service.search(interaction, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        cap = await Cap.fetch_by_channel_guild_and_moderation_type(
            channel_snowflake=channel_obj.id,
            guild_snowflake=interaction.guild.id,
            moderation_type=moderation_type
        )
        if cap and hours:
            await Cap.update_by_channel_and_duration(
                channel_snowflake=channel_obj.id,
                duration=seconds
            )
            msg = f'Cap `{moderation_type}` modified for {channel_obj.mention}.'
        elif cap:
            await Cap.delete_by_channel_guild_and_moderation_type(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                moderation_type=moderation_type
            )
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                    f'Cap of type {moderation_type} ' \
                    f'and channel {channel_obj.mention} deleted successfully.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            cap = Cap(
                channel_snowflake=channel_obj.id,
                duration=seconds,
                guild_snowflake=interaction.guild.id,
                moderation_type=moderation_type
            )
            await cap.create()
            msg = f'Cap `{moderation_type}` created for ' \
                f'{channel_obj.mention} successfully.'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='cap', help='Cap alias duration for mods.')
    @administrator_predicator()
    async def cap_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None,
            description='Tag a channel or include its ID'
        ),
        moderation_type: ModerationType = commands.parameter(
            default=None,
            description='One of: `mute`, `ban`, `tmute`'
        ), *,
        hours: int = commands.parameter(default=24, description='# of hours')
    ):
        state = State(ctx)
        channel_obj = None
        seconds = int(hours) * 3600
        try:
            channel_obj = \
                await self.channel_service.search(ctx, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        cap = await Cap.fetch_by_channel_guild_and_moderation_type(
            channel_snowflake=channel_obj.id,
            guild_snowflake=ctx.guild.id,
            moderation_type=moderation_type
        )
        if cap and hours:
            await Cap.update_by_channel_and_duration(
                channel_snowflake=channel_obj.id,
                duration=seconds
            )
            msg = f'Cap `{moderation_type}` modified for {channel_obj.mention}.'
        elif cap:
            await Cap.delete_by_channel_guild_and_moderation_type(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                moderation_type=moderation_type
            )
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                    f'Cap of type {moderation_type} ' \
                    f'and channel {channel_obj.mention} deleted successfully.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            cap = Cap(
                channel_snowflake=channel_obj.id,
                duration=seconds,
                guild_snowflake=ctx.guild.id,
                moderation_type=moderation_type
            )
            await cap.create()
            msg = f'Cap `{moderation_type}` created for ' \
                f'{channel_obj.mention} successfully.'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(
        name='chown',
        description='Change the owner of a temporary room.'
    )
    @app_commands.describe(
        member='Tag a user or provide their ID',
        channel='Tag a channel or provide it\'s ID'
    )
    @administrator_predicator()
    async def change_temp_room_owner_app_command(
        self,
        interaction,
        channel: AppChannelSnowflake,
        member: AppMemberSnowflake
    ):
        state = State(interaction)
        channel_obj, member_obj = None, None
        try:
            channel_obj = \
                await self.channel_service.search(interaction, channel)
            member_obj = \
                await self.member_service.search(interaction, member)
            check_not_self(interaction, member_snowflake=member_obj.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        await TemporaryRoom.update_owner(
            channel_snowflake=channel_obj.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=member_obj.id
        )
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                f'Temporary room {channel_obj.mention} ownership ' \
                f'transferred to {member_obj.mention}.'
            )
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(
        name='chown',
        help='Change the owner of a temporary room.',
        hidden=True
    )
    @administrator_predicator()
    async def change_temp_room_owner_text_command(
        self,
        ctx,
        channel: ChannelSnowflake = commands.parameter(
            default=None,
            description='Tag a channel or provide it\'s ID'
        ),
        member: MemberSnowflake = commands.parameter(
            default=None,
            description='Tag a user or provide their ID'
        )
    ):
        state = State(ctx)
        channel_obj, member_obj = None, None
        try:
            channel_obj = \
                await self.channel_service.search(ctx, channel)
            member_obj = \
                await self.member_service.search(ctx, member)
            check_not_self(ctx, member_snowflake=member_obj.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        await TemporaryRoom.update_owner(
            channel_snowflake=channel_obj.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=member_obj.id
        )
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                f'Temporary room {channel_obj.mention} ownership ' \
                f'transferred to {member_obj.mention}.'
            )
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @app_commands.command(name='clear', description='Reset channel/member.')
    @app_commands.describe(
        scope='Tag a channel/member or include the ID',
        action_type='Specify one of: `alias`, `all`, `ban`, ' \
            '`coord`, `flag`, `mod`, `temp`, `tmute`, `track`, `vegan`, `vmute` or `vr`.'
    )
    @administrator_predicator()
    async def clear_channel_access_app_command(
        self,
        interaction: discord.Interaction,
        scope: str,
        action_type: str
    ):
        state = State(interaction)
        channel_obj, member_obj = None, None
        is_modification = True
        target = 'user'
        if not action_type:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    "You must specify either `alias`, `all`, `ban`, `coord`, " \
                    "`flag`, `mod`, `temp`, `tmute`, `vegan` or `vmute` or `vr`."
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')  
        try:
            channel_obj = \
                await self.channel_service.search(interaction, scope)
        except Exception as e:
            try:
                member_obj = await self.member_service.search(interaction, scope)
                highest_role = await has_equal_or_higher_role(
                    interaction, channel_snowflake=interaction.channel.id,
                    guild_snowflake=interaction.guild.id,
                    member_snowflake=member_obj.id,
                    sender_snowflake=interaction.user.id
                )
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        highest_role = await permission_check(interaction)
        if channel_obj:
            if highest_role not in ('System Owner', 'Developer', 'Guild Owner'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        'You must be a system owner, a developer ' \
                        'or guild owner to delete channel associations.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            view = VerifyView(
                action_type=action_type,
                author_snowflake=interaction.user.id,
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id
            )
            embed = view.build_embed(
                action_type=view.action_type,
                target=view.target
            )
            await interaction.response.send_message(embed=embed, view=view)
            await view.wait()
            state = State(interaction)
            if view.result == True:
                match action_type.lower():
                    case 'all':
                        await Alias.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        await Ban.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=interaction,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        await Coordinator.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        await Flag.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=interaction,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        await History.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        await Moderator.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        await TemporaryRoom.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        await TextMute.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=interaction,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        await Vegan.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=interaction,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        await VoiceMute.clear_by_channel_guild_highest_role_modification_and_target(
                            ctx_interaction_or_message=interaction,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            s_modification=is_modification,
                            target=target
                        )
                        await VideoRoom.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        msg = f'Deleted all associated aliases, ' \
                            'moderation actions, roles, room setups ' \
                            f'and tracking for {channel_obj.mention}.'
                    case 'alias':
                        await Alias.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        msg = f'Deleted all associated aliases in {channel_obj.mention}.'
                    case 'ban':
                        await Ban.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=interaction,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        msg = f'Deleted all associated bans in {channel_obj.mention}.'
                    case 'coord':
                        await Coordinator.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        msg = f'Deleted all associated coordinators ' \
                            f'in {channel_obj.mention}.'
                    case 'flag':
                        await Flag.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=interaction,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        msg = f'Deleted all associated flags in {channel_obj.mention}.'
                    case 'mod':
                        await Moderator.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        msg = f'Deleted all associated moderators in {channel_obj.mention}.'
                    case 'temp':
                        await TemporaryRoom.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        msg = f'Deleted the associated temporary channel ' \
                            f'for {channel_obj.mention}.'
                    case 'tmute':
                        await TextMute.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=interaction,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        msg = f'Deleted all associated text-mutes in {channel_obj.mention}.'
                    case 'track':
                        await History.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        msg = f'Deleted all associated text-mutes in {channel_obj.mention}.'
                    case 'vegan':
                        await Vegan.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=interaction,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification)
                        msg = f'Deleted all associated new vegans in {channel_obj.mention}.'
                    case 'vmute':
                        await VoiceMute.clear_by_channel_guild_highest_role_modification_and_target(
                            ctx_interaction_or_message=interaction,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            target=target
                        )
                        msg = f'Deleted all associated voice-mutes in {channel_obj.mention}.'
                    case 'vr':
                        await VideoRoom.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=interaction.guild.id
                        )
                        msg = f'Deleted the associated video room in {channel_obj.mention}.'
                    case _:
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                                f'{action_type} is an unknown action type.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
        elif member_obj:
            view = VerifyView(
                action_type=action_type,
                author_snowflake=interaction.user.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id
            )
            embed = view.build_embed(
                action_type=view.action_type,
                target=view.target
            )
            await interaction.response.send_message(embed=embed, view=view)
            await view.wait()
            state = State(interaction)
            if view.result == True:
                match action_type.lower():
                    case 'all':
                        await Ban.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=interaction,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        await Coordinator.delete_by_guild_and_member(guild_snowflake=interaction.guild.id,
                            member_snowflake=member_obj.id
                        )
                        await Flag.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=interaction,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        await Moderator.delete_by_guild_and_member(
                            guild_snowflake=interaction.guild.id,
                            member_snowflake=member_obj.id
                        )
                        await TextMute.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=interaction,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        await Vegan.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=interaction,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        await VoiceMute.clear_by_member_guild_highest_role_modification_and_target(
                            ctx_interaction_or_message=interaction,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id,
                            target=target
                        )
                        msg = f'Deleted all associated moderation actions and ' \
                            f'roles for {member_obj.mention}.'
                    case 'ban':
                        await Ban.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=interaction,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated bans on {member_obj.mention}.'
                    case 'coord':
                        await Coordinator.delete_by_guild_and_member(
                            guild_snowflake=interaction.guild.id,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated coordinator channels for {member_obj.mention}.'
                    case 'flag':
                        await Flag.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=interaction,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated flags on {member_obj.mention}.'
                    case 'mod':
                        await Moderator.delete_by_guild_and_member(
                            guild_snowflake=interaction.guild.id,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated moderator channels for {member_obj.mention}.'
                    case 'tmute':
                        await TextMute.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=interaction,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated text-mutes on {member_obj.mention}.'
                    case 'vegan':
                        await Vegan.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=interaction,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated vegan channels on {member_obj.mention}.'
                    case 'vmute':
                        await VoiceMute.clear_by_guild_highest_role_member_modification_and_target(
                            ctx_interaction_or_message=interaction,
                            guild_snowflake=interaction.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id,
                            target=target
                        )
                        msg = f'Deleted all associated voice-mutes on {member_obj.mention}.'
                    case _:
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F {action_type} ' \
                                'is an unknown action type.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No associated records ' \
                    f'found for scope: {scope}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')    
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='clear', help='Reset channel/member.')
    @administrator_predicator()
    async def clear_channel_access_text_command(
        self,
        ctx: commands.Context,
        scope: str = commands.parameter(
            default=None,
            description='Tag a channel, a member or include the ID'
        ), *,
        action_type: Optional[str] = commands.parameter(
            default=None,
            description='Specify one of: `alias`, `all`, `ban`, `coord`, ' \
                'flag`, `mod`, `temp`, `tmute`, `track`, `vegan`, `vmute` or `vr`.'
        )
    ):
        state = State(ctx)
        channel_obj, member_obj = None, None
        is_modification = True
        target = 'user'
        if not action_type:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    'You must specify either `alias`, `all`, `ban`, `coord`, ' \
                    '`flag`, `mod`, `temp`, `tmute`, `vegan` or `vmute` or `vr`.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')   
        try:
            channel_obj = \
                await self.channel_service.search(ctx, scope)
        except Exception as e:
            try:
                member_obj = await self.member_service.search(ctx, scope)
                highest_role = await has_equal_or_higher_role(
                    ctx,
                    channel_snowflake=ctx.channel.id,
                    guild_snowflake=ctx.guild.id,
                    member_snowflake=member_obj.id,
                    sender_snowflake=ctx.author.id
                )
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')  
        highest_role = await permission_check(ctx)
        if channel_obj:
            if highest_role not in ('System Owner', 'Developer', 'Guild Owner'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        'You must be a system owner, a developer or guild owner ' \
                        'to delete channel associations.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            view = VerifyView(
                action_type=action_type,
                author_snowflake=ctx.author.id,
                guild_snowflake=ctx.guild.id,
                channel_snowflake=channel_obj.id
            )
            embed = view.build_embed(
                action_type=view.action_type,
                target=view.target
            )
            await ctx.send(embed=embed, view=view)
            await view.wait()
            state = State(ctx)
            if view.result == True:
                match action_type.lower():
                    case 'all':
                        await Alias.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        await Ban.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=ctx,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        await Coordinator.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        await Flag.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=ctx,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification)
                        await History.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        await Moderator.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        await TemporaryRoom.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        await TextMute.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=ctx,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        await Vegan.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=ctx,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        await VoiceMute.clear_by_channel_guild_highest_role_modification_and_target(
                            ctx_interaction_or_message=ctx,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            target=target
                        )
                        await VideoRoom.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        msg = f'Deleted all associated aliases, moderation actions, ' \
                            f'roles, room setups and tracking for {channel_obj.mention}.'
                    case 'alias':
                        await Alias.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        msg = f'Deleted all associated aliases in {channel_obj.mention}.'
                    case 'ban':
                        await Ban.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=ctx,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        msg = f'Deleted all associated bans in {channel_obj.mention}.'
                    case 'coord':
                        await Coordinator.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        msg = f'Deleted all associated coordinators in {channel_obj.mention}.'
                    case 'flag':
                        await Flag.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=ctx,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        msg = f'Deleted all associated flags in {channel_obj.mention}.'
                    case 'mod':
                        await Moderator.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        msg = f'Deleted all associated moderators in {channel_obj.mention}.'
                    case 'temp':
                        await TemporaryRoom.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        msg = f'Deleted the associated temporary channel for ' \
                            f'{channel_obj.mention}.'
                    case 'tmute':
                        await TextMute.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=ctx,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        msg = f'Deleted all associated text-mutes in {channel_obj.mention}.'
                    case 'track':
                        await History.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        msg = f'Deleted all associated tracking in {channel_obj.mention}.'
                    case 'vegan':
                        await Vegan.clear_by_channel_guild_highest_role_and_modification(
                            ctx_interaction_or_message=ctx,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification
                        )
                        msg = f'Deleted all associated new vegans in {channel_obj.mention}.'
                    case 'vmute':
                        await VoiceMute.clear_by_channel_guild_highest_role_modification_and_target(
                            ctx_interaction_or_message=ctx,
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            target=target
                        )
                        msg = f'Deleted all associated voice-mutes in {channel_obj.mention}.'
                    case 'vr':
                        await VideoRoom.delete_by_channel_and_guild(
                            channel_snowflake=channel_obj.id,
                            guild_snowflake=ctx.guild.id
                        )
                        msg = f'Deleted the associated video room in {channel_obj.mention}.'
                    case _:
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                                f'{action_type} is an unknown action type.'
                            )
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
        elif member_obj:
            view = VerifyView(
                action_type=action_type,
                author_snowflake=ctx.author.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id
            )
            embed = view.build_embed(
                action_type=view.action_type,
                target=view.target
            )
            await ctx.send(embed=embed, view=view)
            await view.wait()
            state = State(ctx)
            if view.result == True:
                match action_type.lower():
                    case 'all':
                        await Ban.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=ctx,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        await Coordinator.delete_by_guild_and_member(
                            guild_snowflake=ctx.guild.id,
                            member_snowflake=member_obj.id
                        )
                        await Flag.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=ctx,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        await Moderator.delete_by_guild_and_member(
                            guild_snowflake=ctx.guild.id,
                            member_snowflake=member_obj.id
                        )
                        await TextMute.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=ctx,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        await Vegan.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=ctx,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        await VoiceMute.clear_by_guild_highest_role_member_modification_and_target(
                            ctx_interaction_or_message=ctx,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id,
                            target=target
                        )
                        msg = f'Deleted all associated moderation actions and roles ' \
                            f'for {member_obj.mention}.'
                    case 'ban':
                        await Ban.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=ctx,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated bans on {member_obj.mention}.'
                    case 'coord':
                        await Coordinator.delete_by_guild_and_member(
                            guild_snowflake=ctx.guild.id,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated coordinator channels ' \
                            f'for {member_obj.mention}.'
                    case 'flag':
                        await Flag.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=ctx,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated flags on {member_obj.mention}.'
                    case 'mod':
                        await Moderator.delete_by_guild_and_member(
                            guild_snowflake=ctx.guild.id,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated moderator channels for {member_obj.mention}.'
                    case 'tmute':
                        await TextMute.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=ctx,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated text-mutes on {member_obj.mention}.'
                    case 'vegan':
                        await Vegan.clear_by_guild_highest_role_member_and_modification(
                            ctx_interaction_or_message=ctx,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id
                        )
                        msg = f'Deleted all associated vegan channels on {member_obj.mention}.'
                    case 'vmute':
                        await VoiceMute.clear_by_guild_highest_role_member_modification_and_target(
                            ctx_interaction_or_message=ctx,
                            guild_snowflake=ctx.guild.id,
                            highest_role=highest_role,
                            is_modification=is_modification,
                            member_snowflake=member_obj.id,
                            target=target
                        )
                        msg = f'Deleted all associated voice-mutes on {member_obj.mention}.'
                    case _:
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                                f'{action_type} is an unknown action type.'
                            )
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No ' \
                    f'associated records found for scope: {scope}.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @app_commands.command(
        name='coord',
        description='Grant/revoke coords.'
    )
    @app_commands.describe(
        member='Tag a member or include their ID',
        channel='Tag a channel or include its ID'
    )
    @administrator_predicator()
    async def create_coordinator_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj, member_obj = None, None
        try:
            channel_obj = \
                await self.channel_service.search(interaction, channel)
        except Exception as e:
            channel_obj = interaction.channel
            await self.message_service.send_message(
                interaction,
                content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.'
            )
        try:
            member_obj = await self.member_service.search(interaction, member)
            check_not_self(interaction, member_snowflake=member_obj.id)
            highest_role = await has_equal_or_higher_role(
                message_ctx_interaction=interaction,
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id,
                sender_snowflake=interaction.user.id
            )
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        coordinator = await Coordinator.fetch_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=member_obj.id
        )
        if coordinator:
            await Coordinator.delete_by_channel_guild_and_member(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.channel.id,
                member_snowflake=member_obj.id
            )
            action = 'revoked'
        else:
            coordinator = Coordinator(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id
            )
            await coordinator.grant()
            action = 'granted'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                f'Coordinator access has been {action} for {member_obj.mention} ' \
                f'in {channel_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='coord', help='Grant/revoke coords.')
    @administrator_predicator()
    async def create_coordinator_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            default=None,
            description='Tag a member or include their ID'
        ),
        channel: ChannelSnowflake = commands.parameter(
            default=None,
            description='Tag a channel or include its ID'
        )
    ):
        state = State(ctx)
        channel_obj, member_obj = None, None
        try:
            channel_obj = \
                await self.channel_service.search(ctx, channel)
        except Exception as e:
            channel_obj = ctx.channel
            await self.message_service.send_message(
                ctx,
                content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.'
            )
        try:
            member_obj = await self.member_service.search(ctx, member)
            check_not_self(ctx, member_snowflake=member_obj.id)
            highest_role = await has_equal_or_higher_role(
                message_ctx_interaction=ctx,
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id,
                sender_snowflake=ctx.author.id
            )
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        coordinator = await Coordinator.fetch_by_channel_guild_and_member(
            channel_snowflake=channel_obj.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=member_obj.id
        )
        if coordinator:
            await Coordinator.delete_by_channel_guild_and_member(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id
            )
            action = 'revoked'
        else:
            coordinator = Coordinator(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id
            )
            await coordinator.grant()
            action = 'granted'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                f'Coordinator access has been {action} for {member_obj.mention} ' \
                f'in {channel_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # @app_commands.command(name='pc', description='View permissions.')
    # @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, server ID or empty.")

    @commands.command(name='pc', help='View permissions.')
    @administrator_predicator()
    async def list_permissions_text_command(
        self,
        ctx: commands.Context,
        scope: str = commands.parameter(
            default=None,
            description='Specify one of: `all`, channel ' \
                'ID/mention, server ID or empty.'
            )
    ):
        state = State(ctx)
        chunk_size, field_count, lines, pages = 7, 0, [], []
        channel_obj, guild_obj = None, None
        is_at_home = False
        guild_dictionary, skipped_channels, skipped_guilds = {}, set()
        title = f'{self.emoji.get_random_emoji()} Permissions'

        highest_role = await permission_check(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('System Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'You are not authorized to list missing permissions ' \
                        f'across all servers.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        elif scope:
            try:
                channel_obj = \
                    await self.channel_service.search(ctx, scope)
                channels = [channel_obj] 
            except Exception as e:
                if highest_role not in (
                    'System Owner',
                    'Developer',
                    'Guild Owner',
                    'Administrator'
                ):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            'You are not authorized to list all missing permissions in a specific server.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            'Scope must be one of: `all`, channel ID/mention, ' \
                            f'server ID or empty. Received: {scope}.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                channels = guild_obj.channels
        else:
            channels = [ctx.channel]
            guild_obj = ctx.guild
        
        if not channels:
            try:
                if scope:
                    msg = f'No channels exist for scope: {scope}.'
                else:
                    msg = f'No channels exist in {guild_obj.name}.'
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        for channel in channels:
            permissions = channel.permissions_for(guild_obj.me)
            missing = []
            for permission in TARGET_PERMISSIONS:
                if not getattr(permissions, permission):
                    missing.append(permission)
            if not missing:
                continue
            guild_dictionary.setdefault(channel.guild.id, {})[channel.id] = [
                {'Missing Permissions': missing}
            ]
                            
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = \
                dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            current_channel = None
            lines = []
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guilds.add(guild_snowflake)
                continue
            embed = discord.Embed(
                title=title,
                description=guild.name,
                color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channels.setdefault(
                        guild_snowflake,
                        []
                    ).append(channel_snowflake)
                    continue
                channel_lines = []
                if channel_data:
                    for entry in channel_data:
                        for section_name, permissions in entry.items():
                            channel_lines.append(section_name)
                            for permission in permissions:
                                channel_lines.append(f'  ↳ {permission}')
                if not channel_lines:
                    lines.append('')
                    current_channel = channel
                else:
                    i = 0
                    while i < len(channel_lines):
                        remaining_space = chunk_size
                        chunk = channel_lines[i:i + remaining_space]
                        if not lines:
                            current_channel = channel
                        lines.extend(chunk)
                        i += remaining_space
                if len(lines) >= chunk_size:
                    embed.add_field(
                        name=f'Channel: {current_channel.mention}',
                        value='\n'.join(lines),
                        inline=False
                    )
                    pages.append(embed)
                    embed = discord.Embed(
                        title=title,
                        description=f'{guild.name} continued...',
                        color=discord.Color.blue()    
                    )
                    lines = []
                    current_channel = None
            if lines:
                embed.add_field(
                    name=f'Channel: {current_channel.mention}',
                    value='\n'.join(lines),
                    inline=False
                )
            pages.append(embed)
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guilds:
                embed = discord.Embed(
                    title='Skipped Servers',
                    description='\u200b',
                    color=discord.Color.blue()
                )
                lines = []
                for guild_snowflake in skipped_guilds:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            title='Skipped Servers continued...',
                            color=discord.Color.red()
                        )
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_channels:
                for guild_snowflake, channel_list in skipped_channels.items():
                    embed = discord.Embed(
                        color=discord.Color.red(),
                        title=f'Skipped Channels in Server ({guild_snowflake})'
                    )
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(
                                color=discord.Color.red(),
                                title=f'Skipped Channels in Server ' \
                                    f'({guild_snowflake}) continued...'
                            )
                            field_count = 0
                            lines = []
                        lines.append(str(channel_snowflake))
                        field_count += 1
                    embed.description = '\n'.join(lines)
                    pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'Embed size is too large. Limit the scope.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'No permissions found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(name='rmv', description='VC move.')
    @app_commands.describe(
        source_channel='Tag the source channel or include its ID',
        target_channel='Tag the target channel or include its ID'
    )
    @administrator_predicator()
    async def room_move_all_app_command(
        self,
        interaction: discord.Interaction,
        source_channel: AppChannelSnowflake,
        target_channel: AppChannelSnowflake
    ):
        state = State(interaction)
        failed, moved = [], []
        try:
            source_channel_obj = \
                await self.channel_service.search(interaction, source_channel)
            target_channel_obj = \
                await self.channel_service.search(interaction, target_channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
        for member in source_channel_obj.members:
            try:
                await member.move_to(target_channel_obj)
                moved.append(member)
            except discord.Forbidden as e:
                failed.append(member)
                logger.warning(f'Unable to move member ' \
                    f'{member.display_name} ({member.id}) from channel ' \
                    f'{source_channel_obj.name} ({source_channel}) to channel ' \
                    f'{target_channel_obj.name} ({target_channel}) in guild ' \
                    f'{interaction.guild.name} ({interaction.guild.id}).'
                )
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} ' \
                f'Moved {source_channel_obj.mention} to ' \
                f'{target_channel_obj.mention}',
            color=discord.Color.green()
        )
        if moved:
            embed.add_field(
                name=f'Successfully Moved (`{len(moved)}`)',
                value='\n'.join(member.mention for member in moved),
                inline=False
            )
        else:
            embed.add_field(
                name='Successfully Moved',
                value='None',
                inline=False
            )
        if failed:
            embed.add_field(
                name=f'Failed to Move ({len(failed)})',
                value='\n'.join(member.mention for member in failed),
                inline=False
            )
        embed.set_footer(
            text=f'Moved from {source_channel_obj.mention} ' \
                f'to {target_channel_obj.mention}'
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='rmv', help='VC move.')
    @administrator_predicator()
    async def room_move_all_text_command(
        self,
        ctx: commands.Context,
        source_channel: ChannelSnowflake = commands.parameter(
            default=None,
            description='Tag a channel or include its ID'
        ),
        target_channel: ChannelSnowflake = commands.parameter(
            default=None,
            description='Tag a channel or include its ID'
        )
    ):
        state = State(ctx)
        failed, moved = [], []
        try:
            source_channel_obj = \
                await self.channel_service.search(ctx, source_channel)
            target_channel_obj = \
                await self.channel_service.search(ctx, target_channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                f'{str(e).capitalize()}'
            )
        for member in source_channel_obj.members:
            try:
                await member.move_to(target_channel_obj)
                moved.append(member)
            except discord.Forbidden as e:
                failed.append(member)
                logger.warning(f'Unable to move member ' \
                    f'{member.display_name} ({member.id}) from channel ' \
                    f'{source_channel_obj.name} ({source_channel}) to channel ' \
                    f'{target_channel_obj.name} ({target_channel}) in guild ' \
                    f'{ctx.guild.name} ({ctx.guild.id}).'
                )
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} Moved ' \
                f'{source_channel_obj.mention} to ' \
                f'{target_channel_obj.mention}.',
            color=discord.Color.green()
        )
        if moved:
            embed.add_field(
                name=f'Successfully Moved (`{len(moved)}`)',
                value='\n'.join(member.mention for member in moved),
                inline=False
            )
        else:
            embed.add_field(
                name='Successfully Moved',
                value='None',
                inline=False
            )
        if failed:
            embed.add_field(
                name=f'Failed to Move ({len(failed)})',
                value='\n'.join(member.mention for member in failed),
                inline=False
            )
        embed.set_footer(
            text=f'Moved from {source_channel_obj.mention} ' \
                f'to {target_channel_obj.mention}'
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @app_commands.command(
        name='smute',
        description='Server mute/server unmute.'
    )
    @app_commands.describe(
        member='Tag a member or include their ID',
        reason='Optional reason (required for 7 days or more)'
    )
    @administrator_predicator()
    async def toggle_server_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        reason: Optional[str] = 'No reason provided'
    ):
        state = State(interaction)
        member_obj = None
        try:
            member_obj = await self.member_service.search(interaction, member)
            check_not_self(interaction, member_snowflake=member_obj.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        server_mute = await ServerMute.fetch_by_member(
            member_snowflake=member_obj.id
        )
        if not server_mute:
            server_mute = ServerMute(
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id,
                reason=reason
            )
            await server_mute.create()
            action = 'muted'
            should_be_muted = True
        else:
            await ServerMute.delete_by_guild_and_member(
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id
            )
            action = 'unmuted'
            should_be_muted = False
        if member_obj.voice and member_obj.voice.channel:
            try:
                await member_obj.edit(mute=should_be_muted)
            except discord.Forbidden as e:
                return await state.end(error=f'\u274C ' \
                    f'{str(e).capitalize()}'
                )
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                f'Successfully server {action} {member_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(name='smute', help='Server mute/server unmute.')
    @administrator_predicator()
    async def toggle_server_mute_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            default=None,
            description='Tag a member or include their ID'
        ), *,
        reason: Optional[str] = commands.parameter(
            default='No reason provided',
            description='Optional reason (required for 7 days or more)'
        )
    ):
        state = State(ctx)
        member_obj = None
        try:
            member_obj = \
                await self.member_service.search(ctx, member)
            check_not_self(ctx, member_snowflake=member_obj.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        server_mute = await ServerMute.fetch_by_member(
            member_snowflake=member_obj.id
        )
        if not server_mute:
            server_mute = ServerMute(
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id,
                reason=reason
            )
            await server_mute.create()
            action = 'muted'
            should_be_muted = True
        else:
            await ServerMute.delete_by_guild_and_member(
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id
            )
            action = 'unmuted'
            should_be_muted = False
        if member_obj.voice and member_obj.voice.channel:
            try:
                await member_obj.edit(mute=should_be_muted)
            except discord.Forbidden as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                f'Successfully server {action} {member_obj.mention}.'
            )
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @app_commands.command(
        name='stage',
        description='Start/stop stage.'
    )
    @app_commands.describe(
        channel='Tag a voice/stage channel',
        duration='Duration of the stage (e.g., 1h, 30m)'
    )
    @administrator_predicator()
    async def stage_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        duration: AppDuration = DurationObject('1h')
    ):
        state = State(interaction)
        channel_obj = None
        is_modification = False
        failed, skipped, succeeded = [], [], []
        pages = []
        if not isinstance(duration, DurationObject):
            duration = DurationObject(duration)
        if duration.is_modification:
            is_modification = True
        try:
            channel_obj = \
                await self.channel_service.search(interaction, channel)
        except Exception as e:
            channel_obj = interaction.channel
        stage = await Stage.fetch_by_channel_and_guild(
            channel_snowflake=channel_obj.id,
            guild_snowflake=interaction.guild.id
        )
        if is_modification and stage:
            delta = duration.expires_in - datetime.now(timezone.utc)
            if delta.total_seconds() < 0:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'You are not authorized to decrease the duration ' \
                        f'below the current time.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')    
            if stage:
                await Stage.update_duration(
                    channel_snowflake=channel_obj.id,
                    expires_in=duration.expires_in,
                    guild_snowflake=interaction.guild.id
                )
                description_lines = [
                    f'**Channel:** {channel_obj.mention}',
                    f'**Expires:** {duration}'
                ]
                embed = discord.Embed(
                    description='\n'.join(description_lines),
                    title=f'{self.emoji.get_random_emoji()} Stage Modified',
                    color=discord.Color.blurple()
                )
        elif stage:
            title = f'{self.emoji.get_random_emoji()} ' \
                f'Stage Ended in {channel_obj.mention}'
            await Stage.delete_by_channel_and_guild(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id
            )
            for member in channel_obj.members:
                await VoiceMute.delete_by_channel_guild_member_and_target(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=interaction.guild.id,
                    member_snowflake=member.id,
                    target='room'
                )
                voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=interaction.guild.id,
                    member_snowflake=member.id,
                    target='user'
                )
                if not voice_mute and member.voice and member.voice.mute:
                    try:
                        await member.edit(mute=False, reason='Stage ended ' \
                            f'— no user-specific mute found')
                        succeeded.append(member)
                    except discord.Forbidden as e:
                        logger.warning(f'Unable to undo voice-mute ' \
                            f'for member {member.display_name} ({member.id}) in ' \
                            f'channel {channel_obj.name} ({channel_obj.id}) in ' \
                            f'guild {interaction.guild.name} ({interaction.guild.id}).'
                        )
                        failed.append(member)
            description_lines = [
                f'**Channel:** {channel_obj.mention}',
                f'**Unmuted:** {len(succeeded)} users'
            ]
            if failed:
                description_lines.append(f'**Failed:** {len(failed)}')
            embed = discord.Embed(
                description='\n'.join(description_lines),
                title=title,
                color=discord.Color.blurple()
            )
            pages.append(embed)
        else:
            stage = Stage(
                channel_snowflake=channel_obj.id,
                expires_in=duration.expires_in,
                guild_snowflake=interaction.guild.id,
                member_snowflake=interaction.user.id
            )
            await stage.create()
            for member in channel_obj.members:
                if await permission_check_specific(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=interaction.guild.id,
                    member_snowflake=member.id
                ) or member.id == interaction.user.id:
                    skipped.append(member)
                    continue
                voice_mute = await VoiceMute(
                    channel_snowflake=channel_obj.id,
                    expires_in=duration.expires_in,
                    guild_snowflake=interaction.guild.id,
                    member_snowflake=member.id,
                    target='room',
                    reason='Stage mute'
                )
                await voice_mute.create()
                try:
                    if member.voice and member.voice.channel.id == channel_obj.id:
                        await member.edit(mute=True)
                    succeeded.append(member)
                except Exception as e:
                    logger.warning(f'Unable to voice-mute ' \
                        f'member {member.display_name} ({member.id}) ' \
                        f'in channel {channel_obj.name} ({channel_obj.id}) ' \
                        f'in guild {interaction.guild.name} ({interaction.guild.id}).')
                    failed.append(member)
            description_lines = [
                f'**Channel:** {channel_obj.mention}',
                f'**Expires:** {duration}',
                f'**Muted:** {len(succeeded)} users',
                f'**Skipped:** {len(skipped)}'
            ]
            if failed:
                description_lines.append(f'**Failed:** {len(failed)}')
            embed = discord.Embed(
                description='\n'.join(description_lines),
                title=f'{self.emoji.get_random_emoji()} Stage Created in {channel_obj.name}',
                color=discord.Color.blurple()
            )
            pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'Embed size is too large. Limit the scope.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No stages found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(name='stage', help='Start/stop stage')
    @administrator_predicator()
    async def stage_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None,
            description='Tag a channel or include its ID'
        ), *,
        duration: Duration = commands.parameter(
            default=DurationObject('1h'),
            description='Options: (+|-)duration(m|h|d) ' \
                f'0 - permanent / 24h - default')
    ):
        state = State(ctx)
        channel_obj = None
        is_modification = False
        failed, skipped, succeeded = [], [], []
        pages = []
        if not isinstance(duration, DurationObject):
            duration = DurationObject(duration)
        if duration.is_modification:
            is_modification = True
        try:
            channel_obj = \
                await self.channel_service.search(ctx, channel)
        except Exception as e:
            channel_obj = ctx.channel
        stage = await Stage.fetch_by_channel_and_guild(
            channel_snowflake=channel_obj.id,
            guild_snowflake=ctx.guild.id
        )
        if is_modification and stage:
            delta = duration.expires_in - datetime.now(timezone.utc)
            if delta.total_seconds() < 0:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'You are not authorized to decrease the duration ' \
                        f'below the current time.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}') 
            if stage:
                await Stage.update_duration(
                    channel_snowflake=channel_obj.id,
                    expires_in=duration.expires_in,
                    guild_snowflake=ctx.guild.id
                )
                description_lines = [
                    f'**Channel:** {channel_obj.mention}'
                    f'**Expires:** {duration}',
                ]
                embed = discord.Embed(
                    description='\n'.join(description_lines),
                    title=f'{self.emoji.get_random_emoji()} Stage Modified',
                    color=discord.Color.blurple()
                )
        elif stage:
            title = f'{self.emoji.get_random_emoji()} ' \
                f'Stage Ended in {channel_obj.mention}'
            await Stage.delete_by_channel_and_guild(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id
            )
            for member in channel_obj.members:
                await VoiceMute.delete_by_channel_guild_member_and_target(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=ctx.guild.id,
                    member_snowflake=member.id,
                    target='room'
                )
                voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=ctx.guild.id,
                    member_snowflake=member.id,
                    target='user'
                )
                if not voice_mute and member.voice and member.voice.mute:
                    try:
                        await member.edit(mute=False, reason=\
                            f'Stage ended — no user-specific mute found'
                        )
                        succeeded.append(member)
                    except discord.Forbidden as e:
                        logger.warning(f'Unable to undo voice-mute for ' \
                            f'member {member.display_name} ({member.id}) ' \
                            f'in channel {channel_obj.name} ({channel_obj.id}) ' \
                            f'in guild {ctx.guild.name} ({ctx.guild.id}).'
                        )
                        failed.append(member)
            description_lines = [
                f'**Channel:** {channel_obj.mention}',
                f'**Unmuted:** {len(succeeded)} users'
            ]
            if failed:
                description_lines.append(f'**Failed:** {len(failed)}')
            embed = discord.Embed(
                description='\n'.join(description_lines),
                title=title,
                color=discord.Color.blurple()
            )
            pages.append(embed)
        else: 
            stage = Stage(
                channel_snowflake=channel_obj.id,
                expires_in=duration.expires_in,
                guild_snowflake=ctx.guild.id,
                member_snowflake=ctx.author.id
            )
            await stage.create()
            for member in channel_obj.members:
                if await permission_check_specific(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=ctx.guild.id,
                    member_snowflake=member.id
                ) or member.id == ctx.author.id:
                    skipped.append(member)
                    continue
                voice_mute = await VoiceMute(
                    channel_snowflake=channel_obj.id,
                    expires_in=duration.expires_in,
                    guild_snowflake=ctx.guild.id,
                    member_snowflake=member.id,
                    target='room',
                    reason='Stage mute'
                )
                await voice_mute.create()
                try:
                    if member.voice and member.voice.channel.id == channel_obj.id:
                        await member.edit(mute=True)
                    succeeded.append(member)
                except Exception as e:
                    logger.warning(f'Unable to voice-mute ' \
                        f'for member {member.display_name} ({member.id}) ' \
                        f'in channel {channel_obj.name} ({channel_obj.id}) ' \
                        f'in guild {ctx.guild.name} ({ctx.guild.id}).'
                    )
                    failed.append(member)
            description_lines = [
                f'**Channel:** {channel_obj.mention}',
                f'**Expires:** {duration}',
                f'**Muted:** {len(succeeded)} users',
                f'**Skipped:** {len(skipped)}'
            ]
            if failed:
                description_lines.append(f'**Failed:** {len(failed)}')
            embed = discord.Embed(
                description='\n'.join(description_lines),
                title=f'{self.emoji.get_random_emoji()} ' \
                    f'Stage Created in {channel_obj.name}',
                color=discord.Color.blurple()
            )
            pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'Embed size is too large. Limit the scope.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No stages found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @app_commands.command(
        name='temp',
        description='Toggle a temporary room and assign an owner.'
    )
    @app_commands.describe(
        channel='Tag a channel or include its ID',
        owner='Tag a member or include their ID'
    )
    @administrator_predicator()
    async def toggle_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        owner: AppMemberSnowflake
    ):
        state = State(interaction)
        action = None
        channel_obj, member_obj = None, None
        try:
            channel_obj = \
                await self.channel_service.search(interaction, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        temporary_room = await TemporaryRoom.fetch_by_channel_and_guild(
            channel_snowflake=channel_obj.id,
            guild_snowflake=interaction.guild.id
        )
        if temporary_room:
            if temporary_room.member_snowflake:
                await Moderator.delete_by_channel_guild_and_member(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=interaction.guild.id,
                    member_snowflake=temporary_room.member_snowflake
                )
            # await Alias.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            await TemporaryRoom.delete_by_channel_and_guild(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id
            )
            action = 'removed'
        else:
            try:
                member_obj = \
                    await self.member_service.search(interaction, owner)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'{str(e).capitalize()}'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            moderator = Moderator(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id
            )
            await moderator.grant()
            temporary_room = TemporaryRoom(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id,
                room_name=channel_obj.name
            )
            await temporary_room.create()
            action = f'created'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                f'Temporary room {action} in {channel_obj.mention}.'
            )
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(
        name='temp',
        help='Toggle a temporary room and assign an owner.',
        hidden=True
    )
    @administrator_predicator()
    async def toggle_temp_room_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None,
            description='Tag a channel or include its ID'
        ),
        owner: MemberSnowflake = commands.parameter(
            default=None,
            description='Tag a member or include their Discord ID'
        )
    ):
        state = State(ctx)
        action = None
        channel_obj, member_obj = None, None
        try:
            channel_obj = await self.channel_service.search(ctx, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        temporary_room = await TemporaryRoom.fetch_by_channel_and_guild(
            channel_snowflake=channel_obj.id,
            guild_snowflake=ctx.guild.id
        )
        if temporary_room:
            if temporary_room.member_snowflake:
                await Moderator.delete_by_channel_guild_and_member(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=ctx.guild.id,
                    member_snowflake=temporary_room.member_snowflake
                )
            # await Alias.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            await TemporaryRoom.delete_by_channel_and_guild(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id
            )
            action = 'removed' 
        else:
            try:
                member_obj = await self.member_service.search(ctx, owner)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'{str(e).capitalize()}'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            moderator = Moderator(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id
            )
            await moderator.grant()
            temporary_room = TemporaryRoom(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id,
                room_name=channel_obj.name)
            await temporary_room.create()
            action = 'created'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                f'Temporary room has been {action} in {channel_obj.mention}.'
            )
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @app_commands.command(
        name='temps',
        description='List temporary rooms with matching command aliases.'
    )
    @administrator_predicator()
    async def list_temp_rooms_app_command(
        self,
        interaction: discord.Interaction,
        scope: str = None
    ):
        state = State(interaction)
        aliases, temporary_rooms = [], []
        is_at_home = False
        channel_obj, guild_obj = None, None
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary, skipped_channels, skipped_guilds = {}, {}, set()
        title = f'{self.emoji.get_random_emoji()} Temporary Rooms'
        
        highest_role = await permission_check(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('System Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'You are not authorized to list temporary rooms ' \
                        f'across all servers.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            aliases = await Alias.fetch_all()
            temporary_rooms = await TemporaryRoom.fetch_all()
        elif scope:
            try:
                channel_obj = \
                    await self.channel_service.search(interaction, scope) 
                aliases = await Alias.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=interaction.guild.id)
                temporary_room = await TemporaryRoom.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=interaction.guild.id
                )
                temporary_rooms = [temporary_room] if temporary_room else []
            except Exception as e:
                if highest_role not in (
                    'System Owner',
                    'Developer',
                    'Guild Owner',
                    'Administrator'
                ):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'You are not authorized to list temporary rooms ' \
                            f'for specific servers.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'Scope must be one of: `all`, channel ID/mention, ' \
                            f'server ID or empty. Received: {scope}.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                aliases = await Alias.fetch_by_guild(
                    guild_snowflake=int(scope)
                )
                temporary_rooms = await TemporaryRoom.fetch_by_guild(
                    guild_snowflake=int(scope)
                )
        else:
            aliases = await Alias.fetch_by_channel_and_guild(
                channel_snowflake=interaction.channel.id,
                guild_snowflake=interaction.guild.id
            )
            temporary_room = await TemporaryRoom.fetch_by_channel_and_guild(
                channel_snowflake=interaction.channel.id,
                guild_snowflake=interaction.guild.id
            )
            temporary_rooms = [temporary_room] if temporary_room else []
            channel_obj = interaction.channel
            guild_obj = interaction.guild
        
        if not temporary_rooms:
            try:
                if scope:
                    msg = f'No temporary rooms setup for scope: {scope}.'
                else:
                    msg = f'No temporary room setup for {channel_obj.mention} in {guild_obj.name}.'
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        for temporary_room in temporary_rooms:
            guild_id = temporary_room.guild_snowflake
            channel_id = temporary_room.channel_snowflake
            guild_dictionary.setdefault(guild_id, {})
            guild_dictionary[guild_id].setdefault(channel_id, [])
            entry = {}
            guild_dictionary[guild_id][channel_id].append(entry)
            if aliases:
                for alias in aliases:
                    if alias.guild_snowflake == guild_id \
                        and alias.channel_snowflake == channel_id \
                    :
                        entry.setdefault(alias.alias_type, [])
                        entry[alias.alias_type].append(alias.alias_name)
                            
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = \
                dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            current_channel = None
            lines = []
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guilds.add(guild_snowflake)
                continue
            embed = discord.Embed(
                title=title,
                description=guild.name,
                color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channels.setdefault(
                        guild_snowflake,
                        []
                    ).append(channel_snowflake)
                    continue
                channel_lines = []
                if channel_data:
                    for entry in channel_data:
                        for alias_type, alias_names in entry.items():
                            channel_lines.append(f'{alias_type}')
                            for name in alias_names:
                                channel_lines.append(f'  ↳ {name}')
                if not channel_lines:
                    lines.append('')
                    current_channel = channel
                else:
                    i = 0
                    while i < len(channel_lines):
                        remaining_space = chunk_size
                        chunk = channel_lines[i:i + remaining_space]
                        if not lines:
                            current_channel = channel
                        lines.extend(chunk)
                        i += remaining_space
                if len(lines) >= chunk_size:
                    embed.add_field(
                        name=f'Channel: {current_channel.mention}',
                        value='\n'.join(lines),
                        inline=False)
                    pages.append(embed)
                    embed = discord.Embed(
                        title=title,
                        description=f'{guild.name} continued...',
                        color=discord.Color.blue()
                    )
                    lines = []
                    current_channel = None
            if lines:
                embed.add_field(
                    name=f'Channel: {current_channel.mention}',
                    value='\n'.join(lines),
                    inline=False
                )
            pages.append(embed)
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guilds:
                embed = discord.Embed(
                    title='Skipped Servers',
                    description='\u200b',
                    color=discord.Color.blue()
                )
                lines = []
                for guild_snowflake in skipped_guilds:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            title='Skipped Servers continued...',
                            color=discord.Color.red()
                        )
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_channels:
                for guild_snowflake, channel_list in skipped_channels.items():
                    embed = discord.Embed(
                        color=discord.Color.red(),
                        title=f'Skipped Channels in Server ({guild_snowflake})'
                    )
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(
                                color=discord.Color.red(),
                                title=f'Skipped Channels in Server ' \
                                    f'({guild_snowflake}) continued...'
                            )
                            field_count = 0
                            lines = []
                        lines.append(str(channel_snowflake))
                        field_count += 1
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
                
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'Embed size is too large. Limit the scope.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'No temporary rooms found.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(
        name='temps',
        help='List temporary rooms with matching command aliases.',
        hidden=True
    )
    @administrator_predicator()
    async def list_temp_rooms_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(
            default=None,
            description='Specify one of: `all`, channel ID/mention, ' \
                f'server ID or empty.')
    ):
        state = State(ctx)
        aliases, temporary_rooms = [], [], [], []
        is_at_home = False
        channel_obj, guild_obj = None
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary, skipped_channels, skipped_guilds = {}, {}, set()
        title = f'{self.emoji.get_random_emoji()} Temporary Rooms'

        highest_role = await permission_check(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('System Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'You are not authorized to list temporary rooms ' \
                        f'across all servers.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            aliases = await Alias.fetch_all()
            temporary_rooms = await TemporaryRoom.fetch_all()
        elif scope:
            try:
                channel_obj = \
                    await self.channel_service.search(ctx, scope)
                aliases = await Alias.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=ctx.guild.id
                )
                temporary_room = await TemporaryRoom.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=ctx.guild.id
                )
                temporary_rooms = [temporary_room] if temporary_room else []
            except Exception as e:
                if highest_role not in (
                    'System Owner',
                    'Developer',
                    'Guild Owner',
                    'Administrator'
                ):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'You are not authorized to list temporary rooms ' \
                            f'for specific servers.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'Scope must be one of: `all`, channel ID/mention, ' \
                            f'server ID or empty. Received: {scope}.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                aliases = await Alias.fetch_by_guild(guild_snowflake=int(scope))
                temporary_rooms = await TemporaryRoom.fetch_by_guild(
                    guild_snowflake=int(scope)
                )
        else:
            aliases = await Alias.fetch_by_channel_and_guild(
                channel_snowflake=ctx.channel.id,
                guild_snowflake=ctx.guild.id
            )
            temporary_room = await TemporaryRoom.fetch_by_channel_and_guild(
                channel_snowflake=ctx.channel.id,
                guild_snowflake=ctx.guild.id
            )
            temporary_rooms = [temporary_room] if temporary_room else []
            channel_obj = ctx.channel
            guild_obj = ctx.guild
        
        if not temporary_rooms:
            try:
                if scope:
                    msg = f'No temporary rooms setup for scope: {scope}.'
                else:
                    msg = f'No temporary room setup for {channel_obj.mention} in {guild_obj.name}.'
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        for temporary_room in temporary_rooms:
            guild_id = temporary_room.guild_snowflake
            channel_id = temporary_room.channel_snowflake
            guild_dictionary.setdefault(guild_id, {})
            guild_dictionary[guild_id].setdefault(channel_id, [])
            entry = {}
            guild_dictionary[guild_id][channel_id].append(entry)
            if aliases:
                for alias in aliases:
                    if alias.guild_snowflake == guild_id \
                        and alias.channel_snowflake == channel_id \
                    :
                        entry.setdefault(alias.alias_type, [])
                        entry[alias.alias_type].append(alias.alias_name)
                            
        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = \
                dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            current_channel = None
            lines = []
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guilds.add(guild_snowflake)
                continue
            embed = discord.Embed(
                title=title,
                description=guild.name,
                color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channels.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                channel_lines = []
                if channel_data:
                    for entry in channel_data:
                        for alias_type, alias_names in entry.items():
                            channel_lines.append(f'{alias_type}')
                            for name in alias_names:
                                channel_lines.append(f'  ↳ {name}')
                if not channel_lines:
                    lines.append('')
                    current_channel = channel
                else:
                    i = 0
                    while i < len(channel_lines):
                        remaining_space = chunk_size
                        chunk = channel_lines[i:i + remaining_space]
                        if not lines:
                            current_channel = channel
                        lines.extend(chunk)
                        i += remaining_space
                if len(lines) >= chunk_size:
                    embed.add_field(
                        name=f'Channel: {current_channel.mention}',
                        value='\n'.join(lines),
                        inline=False
                    )
                    pages.append(embed)
                    embed = discord.Embed(
                        title=title,
                        description=f'{guild.name} continued...',
                        color=discord.Color.blue()
                    )
                    lines = []
                    current_channel = None
            if lines:
                embed.add_field(
                    name=f'Channel: {current_channel.mention}',
                    value='\n'.join(lines),
                    inline=False
                )
            pages.append(embed)
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guilds:
                embed = discord.Embed(
                    title='Skipped Servers',
                    description='\u200b',
                    color=discord.Color.blue()
                )
                lines = []
                for guild_snowflake in skipped_guilds:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            title='Skipped Servers continued...',
                            color=discord.Color.red()
                        )
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_channels:
                for guild_snowflake, channel_list in skipped_channels.items():
                    embed = discord.Embed(
                        color=discord.Color.red(),
                        title=f'Skipped Channels in Server ({guild_snowflake})'
                    )
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(
                                color=discord.Color.red(),
                                title=f'Skipped Channels in Server ({guild_snowflake}) continued...'
                            )
                            field_count = 0
                            lines = []
                        lines.append(str(channel_snowflake))
                        field_count += 1
                    embed.description = '\n'.join(lines)
                    pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'Embed size is too large. Limit the scope.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'No temporary rooms found.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(name='track', description='Setup tracking.')
    @app_commands.describe(
        channel='Tag a channel or include its ID where the messages will be sent.',
        scope='create | modify | delete.',
        entry_type='all | channel | general.',
        snowflakes='Optional list of channel/member IDs to be tracked.'
    )
    @administrator_predicator()
    async def modify_tracking_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        scope: Optional[str] = None,
        entry_type: Optional[str] = None,
        snowflakes: Optional[str] = None
    ):
        state = State(interaction)
        action = None
        channel_obj = None
        channel_mentions, failed_snowflakes, member_mentions = [], [], []
        enabled = True
        try:
            channel_obj = \
                await self.channel_service.search(interaction, channel) 
        except Exception as e:
            channel_obj = interaction.channel
        if scope is None and entry_type is None:
            history = await History.fetch_by_channel_and_guild(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id
            )
            if history:
                for entry in history:
                    if not entry.enabled:
                        action = 'enabled'
                        enabled = True
                    else:
                        action = 'disabled'
                        enabled = False
                    await History.update_by_channel_enabled_and_guild(
                        channel_snowflake=channel_obj.id,
                        enabled=enabled,
                        guild_snowflake=interaction.guild.id
                    )
                    try:
                        return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                            f'Tracking has been {action} in {channel_obj.mention}.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
            else:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'No tracking exists for {channel_obj.mention}.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}') 
        elif scope and entry_type:
            match scope.lower():
                case 'create':
                    match entry_type:
                        case 'channel':
                            for snowflake in snowflakes:
                                try:
                                    channel = \
                                        await self.channel_service.search(interaction, snowflake)
                                    channel_mentions.append(channel.mention)
                                except Exception as e:
                                    failed_snowflakes.append(snowflake)
                                    continue
                        case 'member':
                            for snowflake in snowflakes:
                                try:
                                    member = \
                                        await self.member_service.search(interaction, snowflake)
                                    member_mentions.append(member.mention)
                                except Exception as e:
                                    failed_snowflakes.append(snowflake)
                                    continue
                    old_history = await History.fetch_by_channel_and_guild(
                        channel_snowflake=channel_obj.id,
                        guild_snowflake=interaction.guild.id
                    )
                    if old_history:
                        enabled = old_history[0].enabled
                    history = History(
                        channel_snowflake=channel_obj.id,
                        enabled=enabled,
                        guild_snowflake=interaction.guild.id,
                        entry_type=entry_type,
                        snowflakes=snowflakes
                    )
                    await history.create()
                    scope = 'created'
                case 'delete':
                    await History.delete_by_channel_and_guild(
                        channel_snowflake=channel_obj.id,
                        guild_snowflake=interaction.guild.id
                    )
                    scope = 'deleted'
                case 'modify':
                    match entry_type:
                        case 'channel':
                            for snowflake in snowflakes:
                                try:
                                    channel = \
                                        await self.channel_service.search(interaction, snowflake)
                                    channel_mentions.append(channel.mention)
                                except Exception as e:
                                    failed_snowflakes.append(snowflake)
                                    continue
                        case 'member':
                            for snowflake in snowflakes:
                                try:
                                    member = \
                                        await self.member_service.search(interaction, snowflake)
                                    member_mentions.append(member.mention)
                                except Exception as e:
                                    failed_snowflakes.append(snowflake)
                                    continue
                    await History.update_by_channel_guild_and_type(
                        channel_snowflake=channel_obj.id,
                        entry_type=entry_type,
                        guild_snowflake=interaction.guild.id,
                        snowflakes=snowflakes
                    )
                    scope = 'modified'
                case _:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'Scope must be one of `create`, `delete` or `modify`.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}') 
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'Scope must be one of `create`, `delete` or `modify`.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')        
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} ' \
                f'Tracking {scope.capitalize()} for {channel_obj.mention}',
                color=0x00FF00
            )
        if channel_mentions:
            embed.add_field(
                name='Processed Channels',
                value=', '.join(channel_mention for channel_mention in channel_mentions),
                inline=False
            )
        if member_mentions:
            embed.add_field(
                name='Processed Members',
                value=', '.join(member_mention for member_mention in member_mentions),
                inline=False
            )
        if failed_snowflakes:
            embed.add_field(
                name='Failed IDs',
                value=', '.join(str(s) for s in failed_snowflakes),
                inline=False
            )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='track', help='Setup tracking.')
    @administrator_predicator()
    async def modify_tracking_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None,
            description='Tag a channel or include its ID where the messages will be sent.'
        ),
        scope: Optional[str] = commands.parameter(
            default=None,
            description='create | modify | delete.'
        ),
        entry_type: Optional[str] = commands.parameter(
            default=None,
            description='all | channel | member.'
        ),
        *snowflakes: Optional[int]
    ):
        state = State(ctx)
        action = None
        channel_obj = None
        channel_mentions, failed_snowflakes, member_mentions = [], [], []
        enabled = True
        try:
            channel_obj = await self.channel_service.search(ctx, channel) 
        except Exception as e:
            channel_obj = ctx.channel
        if scope is None and entry_type is None:
            history = await History.fetch_by_channel_and_guild(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id
            )
            if history:
                for entry in history:
                    if not entry.enabled:
                        action = 'enabled'
                        enabled = True
                    else:
                        action = 'disabled'
                        enabled = False
                    await History.update_by_channel_enabled_and_guild(
                        channel_snowflake=channel_obj.id,
                        enabled=enabled,
                        guild_snowflake=ctx.guild.id
                    )
                    try:
                        return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                            f'Tracking has been {action} in {channel_obj.mention}.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
            else:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'No tracking exists for {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}') 
        elif scope and entry_type:
            match scope.lower():
                case 'create':
                    match entry_type:
                        case 'channel':
                            for snowflake in snowflakes:
                                try:
                                    channel = \
                                        await self.channel_service.search(ctx, snowflake)
                                    channel_mentions.append(channel.mention)
                                except Exception as e:
                                    failed_snowflakes.append(snowflake)
                                    continue
                        case 'member':
                            for snowflake in snowflakes:
                                try:
                                    member = \
                                        await self.member_service.search(ctx, snowflake)
                                    member_mentions.append(member.mention)
                                except Exception as e:
                                    failed_snowflakes.append(snowflake)
                                    continue
                    old_history = await History.fetch_by_channel_and_guild(
                        channel_snowflake=channel_obj.id,
                        guild_snowflake=ctx.guild.id
                    )
                    if old_history:
                        enabled = old_history[0].enabled
                    history = History(
                        channel_snowflake=channel_obj.id,
                        enabled=enabled,
                        guild_snowflake=ctx.guild.id,
                        entry_type=entry_type,
                        snowflakes=snowflakes
                    )
                    await history.create()
                    scope = 'created'
                case 'delete':
                    await History.delete_by_channel_and_guild(
                        channel_snowflake=channel_obj.id,
                        guild_snowflake=ctx.guild.id
                    )
                    scope = 'deleted'
                case 'modify':
                    match entry_type:
                        case 'channel':
                            for snowflake in snowflakes:
                                try:
                                    channel = \
                                        await self.channel_service.search(ctx, snowflake)
                                    channel_mentions.append(channel.mention)
                                except Exception as e:
                                    failed_snowflakes.append(snowflake)
                                    continue
                        case 'member':
                            for snowflake in snowflakes:
                                try:
                                    member = \
                                        await self.member_service.search(ctx, snowflake)
                                    member_mentions.append(member.mention)
                                except Exception as e:
                                    failed_snowflakes.append(snowflake)
                                    continue
                    await History.update_by_channel_guild_and_type(
                        channel_snowflake=channel_obj.id,
                        entry_type=entry_type,
                        guild_snowflake=ctx.guild.id,
                        snowflakes=snowflakes
                    )
                    scope = 'modified'
                case _:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'Scope must be one of `create``, `delete` or `modify`.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}') 
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'Scope must be one of `create`, `delete` or `modify`.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} ' \
                f'Tracking {scope.capitalize()} for ' \
                f'{channel_obj.mention}',
            color=0x00FF00
        )
        if channel_mentions:
            embed.add_field(
                name='Processed Channels',
                value=', '.join(channel_mention for channel_mention in channel_mentions),
                inline=False
            )
        if member_mentions:
            embed.add_field(
                name='Processed Members',
                value=', '.join(member_mention for member_mention in member_mentions),
                inline=False
            )
        if failed_snowflakes:
            embed.add_field(
                name='Failed IDs',
                value=', '.join(str(s) for s in failed_snowflakes),
                inline=False
            )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(name='tracks', description='List history channels.')
    @app_commands.describe(scope="Specify one of: 'all', channel ID/mention, server ID or empty.")
    @administrator_predicator()
    async def list_tracking_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        is_at_home = False
        channel_obj, guild_obj = None, None
        chunk_size, field_count, lines, pages = 7, 0, [], []
        history = []
        guild_dictionary, skipped_channels, skipped_guilds, skipped_snowflakes = {}, {}, set(), []
        title = f'{self.emoji.get_random_emoji()} Logging Routes'

        highest_role = await permission_check(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('System Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'You are not authorized to list logging routes ' \
                        f'across all servers.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            history = await History.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.search(interaction, scope) 
                history = await History.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=interaction.guild.id
                )
            except Exception as e:
                if highest_role not in (
                    'System Owner',
                    'Developer',
                    'Guild Owner',
                    'Administrator'
                ):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'You are not authorized to list logging ' \
                            f'routes for specific servers.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                try:
                    guild_obj = self.bot.get_guild(int(scope))
                    if not guild_obj:
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                                f'Scope must be one of: `all`, channel ID/mention, ' \
                                f'server ID or empty. Received: {scope}.'
                            )
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
                    history = await History.fetch_by_guild(
                        guild_snowflake=int(scope)
                    )
                except:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'Scope must be one of: `all`, channel ID/mention, ' \
                            f'server ID or empty. Received: {scope}.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            history = await History.fetch_by_guild(
                interaction.guild.id
            )
            channel_obj = interaction.channel
            guild_obj = interaction.guild

        if not history:
            try:
                if scope:
                    msg = f'No logging routes setup for scope: {scope}.'
                else:
                    msg = f'No logging routes setup ' \
                        f'for {channel_obj.mention} in {guild_obj.name}.'
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

        for entry in history:
            channel_id = entry.channel_snowflake
            guild_id = entry.guild_snowflake
            guild_dictionary.setdefault(guild_id, {})
            guild_dictionary[guild_id].setdefault(channel_id, [])
            entry_dict = {
                'channels': [],
                'enabled': entry.enabled,
                'entry_snowflakes': entry.snowflakes,
                'entry_type': entry.entry_type,
                'members': []
            }
            guild_dictionary[guild_id][channel_id].append(entry_dict)
            if entry_dict['entry_snowflakes']:
                guild = self.bot.get_guild(guild_id)
                if not guild:
                    continue
                for snowflake in entry.snowflakes:
                    channel = guild.get_channel(snowflake)
                    if channel:
                        entry_dict['channels'].append(channel.mention)
                    else:
                        member = guild.get_member(snowflake)
                        if member:
                            entry_dict['members'].append(member.mention)
                        else:
                            skipped_snowflakes.append(snowflake)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guilds.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channels.setdefault(
                        guild_snowflake,
                        []
                    ).append(channel_snowflake)
                    continue
                lines = []
                for entry_data in channel_data:
                    status = '\u2705' if entry_data['enabled'] else '\u26D4'
                    match entry_data['entry_type']:
                        case 'all':
                            lines.append(f'{status} ' \
                                f'All Actions: {channel.mention}'
                            )
                        case 'channel':
                            channels_from = ', '.join(channel_mention for channel_mention in entry_data.get(
                                'channels',
                                []
                            ))
                            lines.append(f'{status} ' \
                                f'Channel-specific Actions: ' \
                                f'{channel.mention} \n' \
                                f'From: {channels_from}'
                            )
                        case 'member':
                            members_from = ', '.join(member_mention for member_mention in entry_data.get(
                                'members',
                                []
                            ))
                            lines.append(f'{status} ' \
                                f'Member-specific Actions: ' \
                                f'{channel.mention} \n' \
                                f'From: {members_from}'
                            )
                    field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name='Channels',
                        value='\n'.join(lines),
                        inline=False
                    )
                    pages.append(embed)
                    embed = discord.Embed(
                        title=title,
                        description=f'{guild.name} continued...',
                        color=discord.Color.blue()
                    )
                    lines = []
                    field_count = 0
            if lines:
                embed.add_field(
                    name='Channels',
                    value='\n'.join(lines),
                    inline=False)
            pages.append(embed)
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guilds:
                embed = discord.Embed(
                    title='Skipped Servers',
                    description='\u200b',
                    color=discord.Color.blue()
                )
                lines = []
                for guild_snowflake in skipped_guilds:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            title='Skipped Servers continued...',
                            color=discord.Color.red()
                        )
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_channels:
                for guild_snowflake, channel_list in skipped_channels.items():
                    embed = discord.Embed(
                        color=discord.Color.red(),
                        title=f'Skipped Channels in Server ({guild_snowflake})'
                    )
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(
                                color=discord.Color.red(),
                                title=f'Skipped Channels in ' \
                                    f'Server ({guild_snowflake}) continued...'
                                )
                            field_count = 0
                            lines = []
                        lines.append(str(channel_snowflake))
                        field_count += 1
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
            if skipped_snowflakes:
                embed = discord.Embed(
                    color=discord.Color.blue(),
                    title='Skipped Snowflakes for Logging'
                )
                field_count = 0
                lines = []
                for skipped_snowflake in skipped_snowflakes:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            color=discord.Color.blue(),
                            title='Skipped Snowflakes for Logging continued...'
                        )
                        field_count = 0
                        lines = []
                    lines.append(str(skipped_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F '\
                        f'Embed size is too large. Limit the scope.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'No logging routes found.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
               
    # DONE
    @commands.command(name='tracks', help='List history channels.')
    @administrator_predicator()
    async def list_tracking_text_command(
        self,
        ctx: commands.Context,
        *,
        scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', channel ID/mention, server ID or empty.")
    ):
        state = State(ctx)
        is_at_home = False
        channel_obj, guild_obj = None, None
        chunk_size, field_count, lines, pages = 7, 0, [], []
        history = []
        guild_dictionary, skipped_channels, skipped_guilds, skipped_snowflakes = {}, {}, set(), []
        title = f'{self.emoji.get_random_emoji()} Logging Routes'

        highest_role = await permission_check(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('System Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'You are not authorized to list logging ' \
                        f'routes across all servers.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            history = await History.fetch_all()
        elif scope:
            try:
                channel_obj = await self.channel_service.search(ctx, scope) 
                history = await History.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=ctx.guild.id
                )
            except Exception as e:
                if highest_role not in (
                    'System Owner',
                    'Developer',
                    'Guild Owner',
                    'Administrator'
                ):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'You are not authorized to list ' \
                            f'logging routes for specific servers.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'Scope must be one of: `all`, channel ID/mention, ' \
                            f'server ID or empty. Received: {scope}.'
                        )
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                history = await History.fetch_by_guild(
                    guild_snowflake=int(scope)
                )
        else:
            history = await History.fetch_by_guild(ctx.guild.id)
            channel_obj = ctx.channel
            guild_obj = ctx.guild
        
        if not history:
            try:
                if scope:
                    msg = f'No logging routes setup for scope: {scope}.'
                else:
                    msg = f'No logging routes setup for ' \
                        f'{channel_obj.mention} in {guild_obj.name}.'
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
            
        for entry in history:
            channel_id = entry.channel_snowflake
            guild_id = entry.guild_snowflake
            guild_dictionary.setdefault(guild_id, {})
            guild_dictionary[guild_id].setdefault(channel_id, [])
            entry_dict = {
                'channels': [],
                'enabled': entry.enabled,
                'entry_snowflakes': entry.snowflakes,
                'entry_type': entry.entry_type,
                'members': []
            }
            guild_dictionary[guild_id][channel_id].append(entry_dict)
            if entry_dict['entry_snowflakes']:
                guild = self.bot.get_guild(guild_id)
                if not guild:
                    continue
                for snowflake in entry.snowflakes:
                    channel = guild.get_channel(snowflake)
                    if channel:
                        entry_dict['channels'].append(channel.mention)
                    else:
                        member = guild.get_member(snowflake)
                        if member:
                            entry_dict['members'].append(member.mention)
                        else:
                            skipped_snowflakes.append(snowflake)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guilds.add(guild_snowflake)
                continue
            embed = discord.Embed(
                title=title,
                description=guild.name,
                color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channels.setdefault(
                        guild_snowflake,
                        []
                    ).append(channel_snowflake)
                    continue
                lines = []
                for entry_data in channel_data:
                    status = '\u2705' if entry_data['enabled'] else '\u26D4'
                    match entry_data['entry_type']:
                        case 'all':
                            lines.append(f'{status} ' \
                                f'All Actions: {channel.mention}')
                        case 'channel':
                            channels_from = ', '.join(channel_mention for channel_mention in entry_data.get(
                                'channels',
                                []
                            ))
                            lines.append(f'{status} ' \
                                f'Channel-specific Actions: ' \
                                f'{channel.mention} \n' \
                                f'From: {channels_from}'
                            )
                        case 'member':
                            members_from = ', '.join(member_mention for member_mention in entry_data.get(
                                'members',
                                []
                            ))
                            lines.append(f'{status} ' \
                                f'Member-specific Actions: ' \
                                f'{channel.mention} \n' \
                                f'From: {members_from}')
                    field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name='Channels',
                        value='\n'.join(lines),
                        inline=False
                    )
                    pages.append(embed)
                    embed = discord.Embed(
                        title=title,
                        description=f'{guild.name} continued...',
                        color=discord.Color.blue()
                    )
                    lines = []
                    field_count = 0
            if lines:
                embed.add_field(
                    name='Channels',
                    value='\n'.join(lines),
                    inline=False
                )
            pages.append(embed)
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guilds:
                embed = discord.Embed(
                    title='Skipped Servers',
                    description='\u200b',
                    color=discord.Color.blue()
                )
                lines = []
                for guild_snowflake in skipped_guilds:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            title='Skipped Servers continued...',
                            color=discord.Color.red()
                        )
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_channels:
                for guild_snowflake, channel_list in skipped_channels.items():
                    embed = discord.Embed(
                        color=discord.Color.red(),
                        title=f'Skipped Channels in Server ({guild_snowflake})'
                    )
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(
                                color=discord.Color.red(),
                                title=f'Skipped Channels in Server ' \
                                    f'({guild_snowflake}) continued...')
                            field_count = 0
                            lines = []
                        lines.append(str(channel_snowflake))
                        field_count += 1
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
            if skipped_snowflakes:
                embed = discord.Embed(
                    color=discord.Color.blue(),
                    title='Skipped Snowflakes for Logging'
                )
                field_count = 0
                lines = []
                for skipped_snowflake in skipped_snowflakes:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            color=discord.Color.blue(),
                            title='Skipped Snowflakes for ' \
                                f'Logging continued...')
                        field_count = 0
                        lines = []
                    lines.append(str(skipped_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'Embed size is too large. Limit the scope.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'No logging routes found.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(
        name='vr',
        description='Start/stop video-only room.'
    )
    @app_commands.describe(channel='Tag a channel or include the ID')
    @administrator_predicator()
    async def toggle_video_room_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        action = None
        channel_obj = None
        try:
            channel_obj = \
                await self.channel_service.search(interaction, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        video_room = await VideoRoom.fetch_by_channel_and_guild(
            channel_snowflake=channel_obj.id,
            guild_snowflake=interaction.guild.id
        )
        if video_room:
            action = 'removed'
            VideoRoom.video_rooms = [
                vr
                for vr in VideoRoom.video_rooms
                if vr.channel_snowflake != video_room.channel_snowflake
            ]
            await VideoRoom.delete_by_channel_and_guild(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id
            )
        else:
            video_room = VideoRoom(
                channel_snowflake=channel_obj.id,
                guild_snowflake=interaction.guild.id
            )
            await video_room.create()
            VideoRoom.video_rooms.append(video_room)
            action = f'created'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                f'Video-only room {action} in {channel_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    @commands.command(name='vr', help='Start/stop video-only room.')
    @administrator_predicator()
    async def toggle_video_room_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(
            default=None,
            description='Tag a channel or include the ID'
        )
    ):
        state = State(ctx)
        action = None
        channel_obj = None
        try:
            channel_obj = \
                await self.channel_service.search(ctx, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        video_room = await VideoRoom.fetch_by_channel_and_guild(
            channel_snowflake=channel_obj.id,
            guild_snowflake=ctx.guild.id
        )
        if video_room:
            action = 'removed'
            VideoRoom.video_rooms = [
                vr
                for vr in VideoRoom.video_rooms
                if vr.channel_snowflake != video_room.channel_snowflake
            ]
            await VideoRoom.delete_by_channel_and_guild(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id
            )
        else:
            video_room = VideoRoom(
                channel_snowflake=channel_obj.id,
                guild_snowflake=ctx.guild.id
            )
            await video_room.create()
            VideoRoom.video_rooms.append(video_room)
            action = f'created'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} ' \
                f'Video-only room {action} in {channel_obj.mention}.'
            )
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @app_commands.command(
        name='vrs',
        description='List video rooms.'
    )
    @app_commands.describe(
        scope='Specify one of: `all`, channel ID/mention, ' \
            f'server ID or empty.')
    @administrator_predicator()
    async def list_video_rooms_app_command(
        self,
        interaction: discord.Interaction,
        scope: str = None
    ):
        state = State(interaction)
        aliases, video_rooms = [], [], [], []
        is_at_home = False
        channel_obj, guild_obj = None, None
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary, skipped_channels, skipped_guilds = {}, {}, set()
        title = f'{self.emoji.get_random_emoji()} Video Rooms'

        highest_role = await permission_check(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('System Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'You are not authorized to list video rooms ' \
                        f'across all servers.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            video_rooms = await VideoRoom.fetch_all()
        elif scope:
            try:
                channel_obj = \
                    await self.channel_service.search(interaction, scope) 
                aliases = await Alias.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=interaction.guild.id
                )
                video_room = await VideoRoom.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=interaction.guild.id
                )
                video_rooms = [video_room] if video_room else []
            except Exception as e:
                if highest_role not in (
                    'System Owner',
                    'Developer',
                    'Guild Owner',
                    'Administrator'
                ):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'You are not authorized to list ' \
                            f'all video rooms in a specific server.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'Scope must be one of: `all`, channel ID/mention, ' \
                            f'server ID or empty. Received: {scope}.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                aliases = await Alias.fetch_by_guild(guild_snowflake=int(scope))
                video_rooms = await VideoRoom.fetch_by_guild(guild_snowflake=int(scope))
        else:
            aliases = await Alias.fetch_by_channel_and_guild(
                channel_snowflake=interaction.channel.id,
                guild_snowflake=interaction.guild.id
            )
            video_room = await VideoRoom.fetch_by_channel_and_guild(
                channel_snowflake=interaction.channel.id,
                guild_snowflake=interaction.guild.id
            )
            video_rooms = [video_room] if video_room else []
            channel_obj = interaction.channel
            guild_obj = interaction.guild
        
        if not video_rooms:
            try:
                if scope:
                    msg = f'No video rooms setup for scope: {scope}.'
                else:
                    msg = f'No video room setup for {channel_obj.mention} '\
                        f'in {guild_obj.name}.'
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        for video_room in video_rooms:
            guild_id = video_room.guild_snowflake
            channel_id = video_room.channel_snowflake
            guild_dictionary.setdefault(guild_id, {})
            guild_dictionary[guild_id].setdefault(channel_id, [])
            entry = {}
            guild_dictionary[guild_id][channel_id].append(entry)
            if aliases:
                for alias in aliases:
                    if alias.guild_snowflake == guild_id \
                        and alias.channel_snowflake == channel_id \
                    :
                        entry.setdefault(alias.alias_type, [])
                        entry[alias.alias_type].append(alias.alias_name)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = \
                dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            current_channel = None
            lines = []
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guilds.add(guild_snowflake)
                continue
            embed = discord.Embed(
                title=title,
                description=guild.name,
                color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channels.setdefault(
                        guild_snowflake,
                        []
                    ).append(channel_snowflake)
                    continue
                channel_lines = []
                if channel_data:
                    for entry in channel_data:
                        for alias_type, alias_names in entry.items():
                            channel_lines.append(f'{alias_type}')
                            for name in alias_names:
                                channel_lines.append(f'  ↳ {name}')
                if not channel_lines:
                    lines.append('')
                    current_channel = channel
                else:
                    i = 0
                    while i < len(channel_lines):
                        remaining_space = chunk_size
                        chunk = channel_lines[i:i + remaining_space]
                        if not lines:
                            current_channel = channel
                        lines.extend(chunk)
                        i += remaining_space
                if len(lines) >= chunk_size:
                    embed.add_field(
                        name=f'Channel: {current_channel.mention}',
                        value='\n'.join(lines),
                        inline=False
                    )
                    pages.append(embed)
                    embed = discord.Embed(
                        title=title,
                        description=f'{guild.name} continued...',
                        color=discord.Color.blue()
                    )
                    lines = []
                    current_channel = None
            if lines:
                embed.add_field(
                    name=f'Channel: {current_channel.mention}',
                    value='\n'.join(lines),
                    inline=False
                )
            pages.append(embed)
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guilds:
                embed = discord.Embed(
                    title='Skipped Servers',
                    description='\u200b',
                    color=discord.Color.blue()
                )
                lines = []
                for guild_snowflake in skipped_guilds:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            title='Skipped Servers continued...',
                            color=discord.Color.red()
                        )
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_channels:
                for guild_snowflake, channel_list in skipped_channels.items():
                    embed = discord.Embed(
                        color=discord.Color.red(),
                        title=f'Skipped Channels in ' \
                            f'Server ({guild_snowflake})'
                    )
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(
                                color=discord.Color.red(),
                                title=f'Skipped Channels in Server ' \
                                    f'({guild_snowflake}) continued...'
                                )
                            field_count = 0
                            lines = []
                        lines.append(str(channel_snowflake))
                        field_count += 1
                    embed.description = '\n'.join(lines)
                    pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'Embed size is too large. Limit the scope.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'No video rooms found.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @commands.command(name='vrs', help='List video rooms.')
    @administrator_predicator()
    async def list_video_rooms_text_command(
        self,
        ctx: commands.Context,
        *,
        scope: Optional[str] = commands.parameter(
            default=None,
            description='Include `all`, channel or server ID.'
        )
    ):
        state = State(ctx)
        aliases, video_rooms = [], [], [], []
        is_at_home = False
        channel_obj, guild_obj = None, None
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary, skipped_channels, skipped_guilds = {}, {}, set()
        title = f'{self.emoji.get_random_emoji()} Video Rooms'

        highest_role = await permission_check(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('System Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'You are not authorized to list video rooms ' \
                        f'across all servers.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            video_rooms = await VideoRoom.fetch_all()
        elif scope:
            try:
                channel_obj = \
                    await self.channel_service.search(ctx, scope) 
                aliases = await Alias.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=ctx.guild.id
                )
                video_room = await VideoRoom.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=ctx.guild.id
                )
                video_rooms = [video_room] if video_room else []
            except Exception as e:
                if highest_role not in (
                    'System Owner',
                    'Developer',
                    'Guild Owner',
                    'Administrator'
                ):
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'You are not authorized to list ' \
                            f'all video rooms in a specific server.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                guild_obj = self.bot.get_guild(int(scope))
                if not guild_obj:
                    try:
                        return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                            f'Scope must be one of: `all`, channel ID/mention, ' \
                            f'server ID or empty. Received: {scope}.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
                aliases = await Alias.fetch_by_guild(guild_snowflake=int(scope))
                video_rooms = await VideoRoom.fetch_by_guild(guild_snowflake=int(scope))
        else:
            aliases = await Alias.fetch_by_channel_and_guild(
                channel_snowflake=ctx.channel.id,
                guild_snowflake=ctx.guild.id
            )
            video_room = await VideoRoom.fetch_by_channel_and_guild(
                channel_snowflake=ctx.channel.id,
                guild_snowflake=ctx.guild.id
            )
            video_rooms = [video_room] if video_room else []
            channel_obj = ctx.channel
            guild_obj = ctx.guild
        
        if not video_rooms:
            try:
                if scope:
                    msg = f'No video rooms setup for scope: {scope}.'
                else:
                    msg = f'No video room setup for {channel_obj.mention} in {guild_obj.name}.'
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        for video_room in video_rooms:
            guild_id = video_room.guild_snowflake
            channel_id = video_room.channel_snowflake
            guild_dictionary.setdefault(guild_id, {})
            guild_dictionary[guild_id].setdefault(channel_id, [])
            entry = {}
            guild_dictionary[guild_id][channel_id].append(entry)
            if aliases:
                for alias in aliases:
                    if alias.guild_snowflake == guild_id \
                        and alias.channel_snowflake == channel_id \
                    :
                        entry.setdefault(alias.alias_type, [])
                        entry[alias.alias_type].append(alias.alias_name)

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            current_channel = None
            lines = []
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guilds.add(guild_snowflake)
                continue
            embed = discord.Embed(
                title=title,
                description=guild.name,
                color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channels.setdefault(
                        guild_snowflake,
                        []
                    ).append(channel_snowflake)
                    continue
                channel_lines = []
                if channel_data:
                    for entry in channel_data:
                        for alias_type, alias_names in entry.items():
                            channel_lines.append(f'{alias_type}')
                            for name in alias_names:
                                channel_lines.append(f'  ↳ {name}')
                if not channel_lines:
                    lines.append('')
                    current_channel = channel
                else:
                    i = 0
                    while i < len(channel_lines):
                        remaining_space = chunk_size
                        chunk = channel_lines[i:i + remaining_space]
                        if not lines:
                            current_channel = channel
                        lines.extend(chunk)
                        i += remaining_space
                if len(lines) >= chunk_size:
                    embed.add_field(
                        name=f'Channel: {current_channel.mention}',
                        value='\n'.join(lines),
                        inline=False
                    )
                    pages.append(embed)
                    embed = discord.Embed(
                        title=title,
                        description=f'{guild.name} continued...',
                        color=discord.Color.blue()
                    )
                    lines = []
                    current_channel = None
            if lines:
                embed.add_field(
                    name=f'Channel: {current_channel.mention}',
                    value='\n'.join(lines),
                    inline=False
                )
            pages.append(embed)
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guilds:
                embed = discord.Embed(
                    title='Skipped Servers',
                    description='\u200b',
                    color=discord.Color.blue()
                )
                lines = []
                for guild_snowflake in skipped_guilds:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(
                            title='Skipped Servers continued...',
                            color=discord.Color.red()
                        )
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_channels:
                for guild_snowflake, channel_list in skipped_channels.items():
                    embed = discord.Embed(
                        color=discord.Color.red(),
                        title=f'Skipped Channels in ' \
                            f'Server ({guild_snowflake})'
                    )
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(
                                color=discord.Color.red(),
                                title=f'Skipped Channels in Server ' \
                                    f'({guild_snowflake}) continued...'
                                )
                            field_count = 0
                            lines = []
                        lines.append(str(channel_snowflake))
                        field_count += 1
                    embed.description = '\n'.join(lines)
                    pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'Embed size is too large. Limit the scope.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'No video rooms found.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @app_commands.command(name='xalias', description='Delete alias.')
    @administrator_predicator()
    @app_commands.describe(alias_name='Include an alias name')
    async def delete_alias_app_command(
        self,
        interaction: discord.Interaction,
        alias_name: str
    ):
        state = State(interaction)
        channel_obj = None
        alias = await Alias.fetch_by_guild_and_name(
            alias_name=alias_name,
            guild_snowflake=interaction.guild.i
        )
        if not alias:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'No aliases found for `{alias_name}`.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        await Alias.delete_by_guild_and_name(
            alias_name=alias.alias_name,
            guild_snowflake=interaction.guild.id
        )
        try:
            channel_obj = \
                await self.channel_service.search(interaction, alias.channel_snowflake)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        if alias.role_snowflake:
            msg = f'Alias `{alias.alias_name}` of type ' \
                f'`{alias.alias_type}` for channel {channel_obj.mention} ' \
                f' and role {alias.role_mention} deleted successfully.'
        else:
            msg = f'Alias `{alias.alias_name}` of type ' \
                f'`{alias.alias_type}` for channel {channel_obj.mention} ' \
                f'deleted successfully.'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @commands.command(name='xalias', help='Delete alias.')
    @administrator_predicator()
    async def delete_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: str = commands.parameter(
            default=None,
            description='Include an alias name'
        )
    ):
        state = State(ctx)
        channel_obj = None
        alias = await Alias.fetch_by_guild_and_name(
            alias_name=alias_name,
            guild_snowflake=ctx.guild.id
        )
        if not alias:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'No aliases found for `{alias_name}`.'
                )
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        await Alias.delete_by_guild_and_name(
            alias_name=alias.alias_name,
            guild_snowflake=ctx.guild.id
        )
        try:
            channel_obj = \
                await self.channel_service.search(ctx, alias.channel_snowflake)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                    f'{str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        if alias.role_snowflake:
            msg = f'Alias `{alias.alias_name}` of type `{alias.alias_type}` ' \
                f'for channel {channel_obj.mention} ' \
                f'and role {alias.role_mention} deleted successfully.'
        else:
            msg = f'Alias `{alias.alias_name}` of type `{alias.alias_type}` ' \
                f'for channel {channel_obj.mention} deleted successfully.'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
async def setup(bot: DiscordBot):
    await bot.add_cog(AdminCommands(bot))
