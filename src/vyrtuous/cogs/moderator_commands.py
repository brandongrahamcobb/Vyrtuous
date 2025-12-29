''' moderator_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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
from collections import defaultdict
from typing import Optional
from discord import app_commands
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.message_service import MessageService
from vyrtuous.utils.paginator import Paginator
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.all import All
from vyrtuous.utils.ban import Ban
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.duration import Duration, DurationObject
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.flag import Flag
from vyrtuous.utils.setup_logging import logger
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.stage import Stage
from vyrtuous.utils.state import State
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.vegan import Vegan
from vyrtuous.utils.voice_mute import VoiceMute
   
class ModeratorCommands(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.emoji = Emojis()
        self.handler = MessageService(self.bot, self.bot.db_pool)
        self.channel_service = ChannelService()
        self.member_service = MemberService()
    
    # DONE
    @app_commands.command(name='bans', description='Lists ban statistics.')
    @app_commands.describe(scope='"all", channel name/ID/mention, or user mention/ID')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_bans_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, scope)
            if member_obj.id == interaction.guild.me.id:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You cannot determine if {interaction.guild.me.mention} is banned.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            except:
                channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all bans in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            bans = await Ban.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active bans found in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=interaction.guild.id, moderations=bans, moderation_type=Ban)
        elif member_obj:
            bans = await Ban.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active bans found for user {member_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, moderations=bans, moderation_type=Ban)
        elif channel_obj:
            bans = await Ban.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active bans found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=interaction.guild.id, moderations=bans, moderation_type=Ban)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
        
    # DONE
    @commands.command(name='bans', description='Lists ban statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ):
        state = State(ctx)
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, scope)
            if member_obj.id == ctx.guild.me.id:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You cannot determine if {ctx.guild.me.mention} is banned.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            except:
                channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all bans in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            bans = await Ban.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active bans found in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=bans, moderation_type=Ban)
        elif member_obj:
            bans = await Ban.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not banned in any channels.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, moderations=bans, moderation_type=Ban)
        elif channel_obj:
            bans = await Ban.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not bans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active bans found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=ctx.guild.id, moderations=bans, moderation_type=Ban)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    
    # DONE
    @app_commands.command(name='caps', description='List active caps for a channel or all channels if "all" is provided.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    @app_commands.describe(scope='"all", channel name/ID/mention')
    async def list_caps_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, scope)
        except:
            channel_obj = interaction.channel
            await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all caps in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            caps = await Cap.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not caps:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No caps found in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            lines = []
            for cap in caps:
                ch = interaction.guild.get_channel(cap.channel_snowflake)
                ch_name = ch.mention if ch else f'Channel ID `{cap.channel_snowflake}`'
                lines.append(f'**{cap.moderation_type} in {ch_name}** → `{Duration.convert_timedelta_seconds(cap.duration)}`')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                embed = discord.Embed(
                    title="{self.emoji.get_random_emoji()} All Active Caps in Server",
                    description="\n".join(lines[i:i+chunk_size]),
                    color=discord.Color.red()
                )
                pages.append(embed)
        caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
        if not caps:
                return await state.end(warning=f'\U000026A0\U0000FE0F No caps found for {channel_obj.mention}.')
        lines = []
        for cap in caps:
            lines.append(f'**{cap.moderation_type} in {channel_obj.mention}** → `{Duration.convert_timedelta_seconds(cap.duration)}`')
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} Active Caps for {channel_obj.mention}",
            description="\n".join(lines),
            color=discord.Color.red()
        )
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')

    # DONE
    @commands.command(name='caps', help='List active caps for a channel or all channels if "all" is provided.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention')
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, scope)
        except:
            channel_obj = ctx.channel
            await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all caps in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            caps = await Cap.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not caps:
                return await state.end(warning=f'\U000026A0\U0000FE0F No caps found in {ctx.guild.name}.')
            lines = []
            for cap in caps:
                ch = ctx.guild.get_channel(cap.channel_snowflake)
                ch_name = ch.mention if ch else f'Channel ID `{cap.channel_snowflake}`'
                lines.append(f'**{cap.moderation_type} in {ch_name}** → `{Duration.convert_timedelta_seconds(cap.duration)}`')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                embed = discord.Embed(
                    title="{self.emoji.get_random_emoji()} All Active Caps in Server",
                    description="\n".join(lines[i:i+chunk_size]),
                    color=discord.Color.red()
                )
                pages.append(embed)
        caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
        if not caps:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No caps found in {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        lines = []
        for cap in caps:
            lines.append(f'**{cap.moderation_type} in {channel_obj.mention}** → `{Duration.convert_timedelta_seconds(cap.duration)}`')
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} Active Caps for {channel_obj.mention}',
            description='\n'.join(lines),
            color=discord.Color.red()
        )
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
    
 
    # DONE
    @app_commands.command(name='cmds', description='List command aliases routed to a specific channel, temp room, or all channels if "all" is provided.')
    @app_commands.describe(scope='"all", channel name/ID/mention')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_commands_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, scope)
        except:
            channel_obj = interaction.channel
        lines = []
        pages = []
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await state.end(warning=f'\U000026A0\U0000FE0F Only owners, developers or administrators can list all aliases in {interaction.guild.name}.')
            aliases = await Alias.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not aliases:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No aliases found in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            grouped = defaultdict(list)
            for alias in aliases:
                grouped[alias.channel_snowflake].append(alias)
            pages = []
            for channel_id, channel_aliases in grouped.items():
                lines = Alias.format_aliases(channel_aliases)
                channel_obj = interaction.guild.get_channel(channel_id)
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Aliases for {channel_obj.mention}',
                    description='\n'.join(lines),
                    color=discord.Color.blue()
                )
                pages.append(embed)
        else:
            aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not aliases:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No aliases found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            lines.extend(Alias.format_aliases(aliases))
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Aliases in {channel_obj.mention}',
                description='\n'.join(lines),
                color=0x57F287
            )
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        
    # DONE
    @commands.command(name='cmds', help='List command aliases routed to a specific channel, temp room, or all channels if "all" is provided.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or temp room name')
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, scope)
        except:
            channel_obj = ctx.channel
        lines = []
        pages = []
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await state.end(warning=f'\U000026A0\U0000FE0F Only owners, developers or administrators can list all aliases in {ctx.guild.name}.')
            aliases = await Alias.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not aliases:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No aliases found in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            grouped = defaultdict(list)
            for alias in aliases:
                grouped[alias.channel_snowflake].append(alias)
            pages = []
            for channel_id, channel_aliases in grouped.items():
                lines = Alias.format_aliases(channel_aliases)
                channel_obj = ctx.guild.get_channel(channel_id)
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Aliases for {channel_obj.mention}',
                    description='\n'.join(lines),
                    color=discord.Color.blue()
                )
                pages.append(embed)
        else:
            aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not aliases:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No aliases found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            lines.extend(Alias.format_aliases(aliases))
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Aliases in {channel_obj.mention}',
                description='\n'.join(lines),
                color=0x57F287
            )
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
 
    # DONE
    @app_commands.command(name='del', description='Delete a message by ID (only if you are coordinator/moderator of that temp room).')
    @app_commands.describe(
        message='Message snowflake ID',
        channel='Tag a channel or include its snowflake ID'
    )
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def delete_message_app_command(
        self,
        interaction: discord.Interaction,
        message: AppMessageSnowflake,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
            await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        member_permission_role = await is_owner_developer_administrator_coordinator_moderator_via_channel_member(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=interaction.user.id)
        if member_permission_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator', 'Moderator'):
            return await state.end(warning='\U000026A0\U0000FE0F You are not permitted to delete messages in {channel_obj.mention}.')
        msg = await channel_obj.fetch_message(message)
        if not msg:
            return await state.end(warning=f'\U000026A0\U0000FE0F Message `{message}` does not exist.')
        try:
            await msg.delete()
        except discord.Forbidden:
            return await state.end(warning='\U000026A0\U0000FE0F Missing permissions to delete the message.')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Message `{message}` deleted successfully.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')

    # DONE
    @commands.command(name='del', help='Delete a message by ID (only if you are coordinator/moderator of that temp room).')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        message: MessageSnowflake = commands.parameter(default=None, description='Message snowflake'),
        *,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Channel or snowflake')
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
            await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        member_permission_role = await is_owner_developer_administrator_coordinator_moderator_via_channel_member(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=ctx.author.id)
        if member_permission_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator', 'Moderator'):
            return await state.end(warning='\U000026A0\U0000FE0F You are not permitted to delete messages in {channel_obj.mention}.')
        msg = await channel_obj.fetch_message(message)
        if not msg:
            return await state.end(warning=f'\U000026A0\U0000FE0F Message `{message}` does not exist.')
        try:
            await msg.delete()
        except discord.Forbidden:
            return await state.end(warning='\U000026A0\U0000FE0F Missing permissions to delete the message.')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Message `{message}` deleted successfully.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
        
    # DONE
    @app_commands.command(name='flags', description='List flag statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_flags_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, scope)
            if member_obj.id == interaction.guild.me.id:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You cannot determine if {interaction.guild.me.mention} is flagged.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            except:
                channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all flags in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            flags = await Flag.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active flags found in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=interaction.guild.id, moderations=flags, moderation_type=Flag)
        elif member_obj:
            flags = await Flag.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not flagged in any channels.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, moderations=flags, moderation_type=Flag)
        elif channel_obj:
            flags = await Flag.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active flags found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=interaction.guild.id, moderations=flags, moderation_type=Flag)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    
    # DONE
    @commands.command(name='flags', help='List flag statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ):
        state = State(ctx)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, scope)
            if member_obj.id == ctx.guild.me.id:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You cannot determine if {ctx.guild.me.mention} is flagged.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            except:
                channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all flags in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            flags = await Flag.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active flags found in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=flags, moderation_type=Flag)
        elif member_obj:
            flags = await Flag.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not flagged in any channels.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, moderations=flags, moderation_type=Flag)
        elif channel_obj:
            flags = await Flag.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not flags:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active flags found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=ctx.guild.id, moderations=flags, moderation_type=Flag)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
        
    # DONE
    @app_commands.command(name='ls', description='List users veganed as going vegan in this guild.')
    @app_commands.describe(scope='A member or channel snowflake ID/mention')
    async def list_vegans_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, scope)
            if member_obj.id == interaction.guild.me.id:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You cannot determine if {interaction.guild.me.mention} is a new vegan.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            except:
                channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all new vegans in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            vegans = await Vegan.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No new vegans found in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=interaction.guild.id, moderations=vegans, moderation_type=Vegan)
        elif member_obj:
            vegans = await Vegan.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No new vegans found for user {member_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, moderations=vegans, moderation_type=Vegan)
        elif channel_obj:
            vegans = await Vegan.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active vegans found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=interaction.guild.id, moderations=vegans, moderation_type=Vegan)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
                            
    # DONE
    @commands.command(name='ls', help='List users veganed as going vegan in this guild.')
    async def list_members_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, scope)
            if member_obj.id == ctx.guild.me.id:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You cannot determine if {ctx.guild.me.mention} is a new vegan.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            except:
                channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all new vegans in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            vegans = await Vegan.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No new vegans found in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=vegans, moderation_type=Vegan)
        elif member_obj:
            vegans = await Vegan.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No new vegans found for user {member_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, moderations=vegans, moderation_type=Vegan)
        elif channel_obj:
            vegans = await Vegan.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not vegans:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active vegans found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=ctx.guild.id, moderations=vegans, moderation_type=Vegan)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')

 # DONE
    @app_commands.command(name='migrate', description='Migrate a temporary room to a new channel.')
    @app_commands.describe(
        old_name='Old temporary room name',
        channel='New channel to migrate to'
    )
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def migrate_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        old_name: str,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        old_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=interaction.guild.id, room_name=old_name)
        if old_room:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, channel)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            is_owner = old_room.member_snowflake == interaction.user.id
            highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
            if highest_role not in ('Owner', 'Developer', 'Administrator') or not is_owner:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners, developers and administrators can migrate rooms.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            await TemporaryRoom.update_by_source_and_target(guild_snowflake=interaction.guild.id, room_name=channel_obj.id, source_channel_snowflake=old_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            new_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=interaction.guild.id, room_name=channel_obj.name)
            await Alias.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Ban.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Cap.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Coordinator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Moderator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Stage.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await TextMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await VoiceMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Temporary room {old_name} migrated to {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        try:
            return await state.end(warning=f'\U000026A0\U0000FE0F No temporary rooms found called {old_name} in {interaction.guild.name}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    
    # DONE
    @commands.command(name='migrate', help='Migrate a temporary room to a new channel by snowflake.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def migrate_temp_room_text_command(
        self,
        ctx: commands.Context,
        old_name: Optional[str] = commands.parameter(default=None, description='Provide a channel name'),
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        channel_obj = None
        old_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=ctx.guild.id, room_name=old_name)
        if old_room:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, channel)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            is_owner = old_room.member_snowflake == ctx.author.id
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
            if highest_role not in ('Owner', 'Developer', 'Administrator') or not is_owner:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners, developers and administrators can migrate rooms.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            await TemporaryRoom.update_by_source_and_target(guild_snowflake=ctx.guild.id, room_name=channel_obj.id, source_channel_snowflake=old_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            new_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=ctx.guild.id, room_name=channel_obj.name)
            await Alias.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Ban.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Cap.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Coordinator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Moderator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Stage.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await TextMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await VoiceMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Temporary room {old_name} migrated to {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        try:
            return await state.end(warning=f'\U000026A0\U0000FE0F No temporary rooms found called {old_name} in {ctx.guild.name}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    # DONE
    @app_commands.command(name='mutes', description='Lists mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_mutes_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, scope)
            if member_obj.id == interaction.guild.me.id:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You cannot determine if {interaction.guild.me.mention} is muted.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            except:
                channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all voice mutes in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=interaction.guild.id, target="user")
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=interaction.guild.id, moderations=voice_mutes, moderation_type=VoiceMute)
        elif member_obj:
            voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, target="user")
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found for user {member_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, moderations=voice_mutes, moderation_type=VoiceMute)
        elif channel_obj:
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, target="user")
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=interaction.guild.id, moderations=voice_mutes, moderation_type=VoiceMute)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
            
    # DONE
    @commands.command(name='mutes', help='Lists mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ):
        state = State(ctx)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, scope)
            if member_obj.id == ctx.guild.me.id:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You cannot list flags on the bot.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            except:
                channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners or developers can list all voice mutes in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=ctx.guild.id, target="user")
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=voice_mutes, moderation_type=VoiceMute)
        elif member_obj:
            voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, target="user")
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found for user {member_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, moderations=voice_mutes, moderation_type=VoiceMute)
        elif channel_obj:
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, target="user")
            if not voice_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active voice mutes found in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(guild_snowflake=ctx.guild.id, moderations=voice_mutes, moderation_type=VoiceMute)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    
    # DONE
    @app_commands.command(name='mstage', description='Mute/unmute a member in the active stage.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    async def stage_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
            await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        try:
            member_obj = await self.member_service.resolve_member(interaction, member)
            await has_equal_or_higher_role(interaction, channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, sender_snowflake=interaction.user.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
        if not stage:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No active stage found.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        try:
            await member_obj.edit(mute=not member_obj.voice.mute)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been {"muted" if member_obj.voice.mute else "unmuted"}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
                
    # DONE
    @commands.command(name='mstage', help='Mute/unmute a member in the active stage.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def stage_mute_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: ChannelSnowflake = commands.parameter(default=None, description="Tag a channel or include it's snowflake ID")
    ):
        state = State(ctx)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
            await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        try:
            member_obj = await self.member_service.resolve_member(ctx, member)
            await has_equal_or_higher_role(ctx, channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, sender_snowflake=ctx.author.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
        if not stage:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No active stage found.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        try:
            await member_obj.edit(mute=not member_obj.voice.mute)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been {"muted" if member_obj.voice.mute else "unmuted"}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    
    # DONE
    @app_commands.command(name='stages', description='Lists stage mute statistics.')
    @app_commands.describe(scope='"all", channel name/ID/mention')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_stages_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, scope)
        except:
            channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners/devs can view all stages in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            stages = await Stage.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not stages:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active stages in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages, chunk_size = [], 8
            for i in range(0, len(stages), chunk_size):
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Active Stages in {interaction.guild.name}',
                    color=discord.Color.purple()
                )
                for s in stages[i:i+chunk_size]:
                    ch = interaction.guild.get_channel(s.channel_snowflake)
                    ch_name = ch.mention
                    voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=interaction.guild.id, target="room")
                    for voice_mute in voice_mutes:
                        embed.add_field(name=ch_name, value=f'Active stage mutes: {voice_mute.member_snowflake}', inline=False)
                pages.append(embed)
        else:
            stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not stage:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active stages in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, target="room")
            initiator = interaction.guild.get_member(stage.member_snowflake)
            initiator_name = initiator.mention
            expires = DurationObject.from_expires_at(stage.expires_at) if stage.expires_at else 'No expiration'
            lines = []
            for m in voice_mutes:
                user = interaction.guild.get_member(m.member_snowflake)
                duration_str = DurationObject.from_expires_at(m.expires_at) if m.expires_at else 'No expiration'
                reason = m.reason or 'No reason provided'
                lines.append(f'• {user.mention} — {reason} — {duration_str}')
            description = f'Initiator: {initiator_name}\nStage Expires: {expires}\nActive stage mutes: {len(lines)}\n\n'+'\n'.join(lines)
            pages, chunk_size = [], 18
            for i in range(0, len(description.splitlines()), chunk_size):
                embed=discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Stage info for {channel_obj.mention}',
                    description='\n'.join(description.splitlines()[i:i+chunk_size]),
                    color=discord.Color.purple()
                )
                pages.append(embed)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')          
    # DONE
    @commands.command(name='stages', help='Lists stage mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_stages_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention')
    ):
        state = State(ctx)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, scope)
        except:
            channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners/devs can view all stages in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            stages = await Stage.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not stages:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active stages in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages, chunk_size = [], 8
            for i in range(0, len(stages), chunk_size):
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Active Stages in {ctx.guild.name}',
                    color=discord.Color.purple()
                )
                for s in stages[i:i+chunk_size]:
                    ch = ctx.guild.get_channel(s.channel_snowflake)
                    ch_name = ch.mention
                    voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=ctx.guild.id, target="room")
                    for voice_mute in voice_mutes:
                        embed.add_field(name=ch_name, value=f'Active stage mutes: {voice_mute.member_snowflake}', inline=False)
                pages.append(embed)
        else:
            stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not stage:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No active stages in {channel_obj.mention}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, target="room")
            initiator = ctx.guild.get_member(stage.member_snowflake)
            initiator_name = initiator.mention
            expires = DurationObject.from_expires_at(stage.expires_at) if stage.expires_at else 'No expiration'
            lines = []
            for m in voice_mutes:
                user = ctx.guild.get_member(m.member_snowflake)
                duration_str = DurationObject.from_expires_at(m.expires_at) if m.expires_at else 'No expiration'
                reason = m.reason or 'No reason provided'
                lines.append(f'• {user.mention} — {reason} — {duration_str}')
            description = f'Initiator: {initiator_name}\nStage Expires: {expires}\nActive stage mutes: {len(lines)}\n\n'+'\n'.join(lines)
            pages, chunk_size = [], 18
            for i in range(0, len(description.splitlines()), chunk_size):
                embed=discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Stage info for {channel_obj.mention}',
                    description='\n'.join(description.splitlines()[i:i+chunk_size]),
                    color=discord.Color.purple()
                )
                pages.append(embed)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
            
    # DONE
    @app_commands.command(name='tmutes', description='Lists text-mute statistics.')
    @app_commands.describe(scope='"all", channel name/ID/mention, or user mention/ID')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_text_mutes_app_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str] = None
    ):
        state = State(interaction)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, scope)
            if member_obj.id == interaction.guild.me.id:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You cannot determine if {interaction.guild.me.mention} is flagged.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            except:
                channel_obj = interaction.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners, developers and administrators can list all text-mutes in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            text_mutes = await TextMute.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No users are currently text-muted in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=text_mutes, moderation_type=TextMute)
        elif member_obj:
            text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not text-muted in any channels.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=text_mutes, moderation_type=TextMute)
        elif channel_obj:
            text_mutes = await TextMute.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No text-muted users currently in {channel_obj.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, moderations=text_mutes, moderation_type=TextMute)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
            
    # DONE
    @commands.command(name='tmutes', help='Lists text-mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        scope: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ):
        state = State(ctx)
        pages = []
        channel_obj = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, scope)
            if member_obj.id == ctx.guild.me.id:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F You cannot determine if {ctx.guild.me.mention} is text-muted.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            except:
                channel_obj = ctx.channel
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if scope and scope.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Only owners, developers and administrators can list all text-mutes in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            text_mutes = await TextMute.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No users are currently text-muted in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild(guild_snowflake=ctx.guild.id, moderations=text_mutes, moderation_type=TextMute)
        elif member_obj:
            text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.mention} is not text-muted in any channels.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, moderations=text_mutes, moderation_type=TextMute)
        elif channel_obj:
            text_mutes = await TextMute.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not text_mutes:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No text-muted users currently in {channel_obj.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_from_moderations_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, moderations=text_mutes, moderation_type=TextMute)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
            
async def setup(bot: DiscordBot):
    cog = ModeratorCommands(bot)
    await bot.add_cog(cog)

