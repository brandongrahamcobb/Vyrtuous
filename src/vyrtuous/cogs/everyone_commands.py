''' commands.py

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
from discord import app_commands
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.message_service import MessageService
from vyrtuous.utils.paginator import Paginator
from vyrtuous.utils.administrator import Administrator
from vyrtuous.utils.all import All
from vyrtuous.utils.developer import Developer
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.moderator import Moderator
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.state import State
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.utils.setup_logging import logger
from vyrtuous.utils.time_to_complete import TimeToComplete
   
import discord
import time

class EveryoneCommands(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.emoji = Emojis()
        self.handler = MessageService(self.bot, self.bot.db_pool)
        self.channel_service = ChannelService()
        self.member_service = MemberService()
        
    # DONE
    @app_commands.command(name='admins', description='Lists all members with server mute privileges in this guild.')
    async def list_administrators_app_command(
        self,
        interaction: discord.Interaction,
        target: str
    ):
        state = State(interaction)
        pages = []
        administrators = await Administrator.fetch_members_by_guild(guild_snowflake=interaction.guild.id)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all admins across all guilds.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_to_show_guilds_by_members(members=administrators, member_type=Administrator)
        else:
            pages = await All.create_pages_to_show_members_by_guild(guild_snowflake=interaction.guild.id, members=administrators, member_type=Administrator)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No admins found in {interaction.guild.name}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
            
    # DONE
    @commands.command(name='admins', help='Lists all members with server mute privileges in this guild.')
    async def list_administrators_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(default=None, description='"all", or user mention/ID')
    ):
        state = State(ctx)
        pages = []
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, target)
        except:
            pass
        administrators = await Administrator.fetch_members_by_guild(guild_snowflake=ctx.guild.id)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all admins across all guilds.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            pages = await All.create_pages_to_show_guilds_by_members(members=administrators, member_type=Administrator)
        elif member_obj:
            guild_snowflakes = await Administrator.fetch_guilds_by_member(member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_guilds_by_member(guild_snowflakes=guild_snowflakes, member_snowflake=member_obj.id, member_type=Developer)
        else:
            pages = await All.create_pages_to_show_members_by_guild(guild_snowflake=ctx.guild.id, members=administrators, member_type=Administrator)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No admins found in {ctx.guild.name}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
    # DONE
    @app_commands.command(name='coords', description='Lists coordinators for a specific voice channel, all, or a member.')
    @app_commands.describe(target='"all", member or channel name/ID/mention')
    async def list_coordinators_app_command(
        self,
        interaction : discord.Interaction,
        target : Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        member_obj = None
        pages = []
        try:
            member_obj = await self.member_service.resolve_member(interaction, target)
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, target)
            except:
                channel_obj = interaction.channel
                await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all coordinators in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            coordinators = await Coordinator.fetch_members_by_guild(guild_snowflake=interaction.guild.id)
            pages = await All.create_pages_to_show_channels_by_guild_and_members(guild_snowflake=interaction.guild.id, members=coordinators, member_type=Coordinator)
        elif member_obj:
            channel_snowflakes = await Coordinator.fetch_channels_by_guild_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_channels_by_guild_and_member(channel_snowflakes=channel_snowflakes, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, member_type=Coordinator)
        elif channel_obj:
            coordinators = await Coordinator.fetch_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            pages = await All.create_pages_to_show_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, members=coordinators, member_type=Coordinator)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        try:
            return await state.end(warning=f'\U000026A0\U0000FE0F No coordinators found for {target}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')

    # DONE
    @commands.command(name='coords', help='Lists coordinators for a specific voice channel, all, or a member.')
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Voice channel name, mention, ID, "all", or member ID')
    ):
        state = State(ctx)
        channel_obj = None
        member_obj = None
        pages = []
        try:
            member_obj = await self.member_service.resolve_member(ctx, target)
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, target)
            except:
                channel_obj = ctx.channel
                await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all coordinators in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            coordinators = await Coordinator.fetch_members_by_guild(guild_snowflake=ctx.guild.id)
            pages = await All.create_pages_to_show_members_by_guild(guild_snowflake=ctx.guild.id, members=coordinators, member_type=Coordinator)
        elif member_obj:
            channel_snowflakes = await Coordinator.fetch_channels_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_channels_by_guild_and_member(channel_snowflakes=channel_snowflakes, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, member_type=Coordinator)
        elif channel_obj:
            coordinators = await Coordinator.fetch_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            pages = await All.create_pages_to_show_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, members=coordinators, member_type=Coordinator)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        try:
            return await state.end(warning=f'\U000026A0\U0000FE0F No coordinators found for {target}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    
    # DONE
    @app_commands.command(name='devs', description='Lists developers.')
    @app_commands.describe(target='"all" or a member snowflake ID/mention')
    async def list_developers_app_command(
        self,
        interaction : discord.Interaction,
        target : Optional[str] = None
    ):
        state = State(interaction)
        member_obj = None
        pages = []
        try:
            member_obj = await self.member_service.resolve_member(interaction, target)
        except:
            pass
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all developers across all guilds.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            developers = await Developer.fetch_all()
            pages = await All.create_pages_to_show_guilds_by_members(members=developers, member_type=Developer)
        elif member_obj:
            guild_snowflakes = await Developer.fetch_guilds_by_member(member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_guilds_by_member(guild_snowflakes=guild_snowflakes, member_snowflake=member_obj.id, member_type=Developer)
        else:
            developers = await Developer.fetch_members_by_guild(guild_snowflake=interaction.guild.id)
            pages = await All.create_pages_to_show_members_by_guild(guild_snowflake=interaction.guild.id, members=developers, member_type=Developer)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No developers found in {interaction.guild.name}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        
    # DONE
    @commands.command(name='devs', help='Lists developers.')
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(default=None, description='"all", or user mention/ID')
    ):
        state = State(ctx)
        member_obj = None
        pages = []
        try:
            member_obj = await self.member_service.resolve_member(ctx, target)
        except:
            pass
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all developers across all guilds.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            developers = await Developer.fetch_all()
            pages = await All.create_pages_to_show_guilds_by_members(members=developers, member_type=Developer)
        elif member_obj:
            guild_snowflakes = await Developer.fetch_guilds_by_member(member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_guilds_by_member(guild_snowflakes=guild_snowflakes, member_snowflake=member_obj.id, member_type=Developer)
        else:
            developers = await Developer.fetch_members_by_guild(guild_snowflake=ctx.guild.id)
            pages = await All.create_pages_to_show_members_by_guild(guild_snowflake=ctx.guild.id, members=developers, member_type=Developer)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No developers found in {ctx.guild.name}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        
    # DONE
    @app_commands.command(name='mods', description='Lists moderator statistics.')
    @app_commands.describe(target='A member or channel snowflake ID/mention')
    async def list_moderators_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        member_obj = None
        pages = []
        try:
            member_obj = await self.member_service.resolve_member(interaction, target)
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, target)
            except:
                channel_obj = interaction.channel
                await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all moderators in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            moderators = await Moderator.fetch_members_by_guild(guild_snowflake=interaction.guild.id)
            pages = await All.create_pages_to_show_members_by_guild(guild_snowflake=interaction.guild.id, members=moderators, member_type=Moderator)
        elif member_obj:
            channel_snowflakes = await Moderator.fetch_channels_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_channels_by_guild_and_member(channel_snowflakes=channel_snowflakes, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, member_type=Moderator)
        elif channel_obj:
            moderators = await Moderator.fetch_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            pages = await All.create_pages_to_show_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, members=moderators, member_type=Moderator)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No moderators found in {interaction.guild.name}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
    # DONE
    @commands.command(name='mods',help='Lists moderator statistics.')
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Voice channel name/mention/ID, "all", or member mention/ID')
    ):
        state = State(ctx)
        channel_obj = None
        member_obj = None
        pages = []
        try:
            member_obj = await self.member_service.resolve_member(ctx, target)
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, target)
            except:
                channel_obj = ctx.channel
                await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all moderators in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            moderators = await Moderator.fetch_members_by_guild(guild_snowflake=ctx.guild.id)
            pages = await All.create_pages_to_show_members_by_guild(guild_snowflake=ctx.guild.id, members=moderators, member_type=Moderator)
        elif member_obj:
            channel_snowflakes = await Moderator.fetch_channels_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_channels_by_guild_and_member(channel_snowflakes=channel_snowflakes, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, member_type=Moderator)
        elif channel_obj:
            moderators = await Moderator.fetch_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            pages = await All.create_pages_to_show_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, members=moderators, member_type=Moderator)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No moderators found in {ctx.guild.name}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
    # DONE
    @app_commands.command(name='owners', description='Show temporary room stats for "all", a channel, or a member.')
    @app_commands.describe(target='"all", a channel mention/ID, or a member mention/ID')
    async def temp_room_stats_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        member_obj = None
        pages = []
        try:
            member_obj = await self.member_service.resolve_member(interaction, target)
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, target)
            except:
                channel_obj = interaction.channel
                await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all temporary room owners in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            rooms = await TemporaryRoom.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not rooms:
                try:
                    return await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F No temporary rooms exist.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            chunk = 12
            for i in range(0, len(rooms), chunk):
                subset = rooms[i:i+chunk]
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Temporary Rooms',
                    color=discord.Color.blurple()
                )
                for room in subset:
                    embed.add_field(
                        name=room.room_name,
                        value=f'• Owner: {room.member_mention}\n• Channel: {room.channel_mention}',
                        inline=False
                    )
                pages.append(embed)
        elif member_obj:
            rooms = await TemporaryRoom.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            if not rooms:
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.display_name} does not own any temporary rooms in {interaction.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            chunk = 12
            for i in range(0, len(rooms), chunk):
                subset = rooms[i:i+chunk]
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Temporary Rooms Owned by {member_obj.display_name}',
                    color=discord.Color.blurple()
                )
                for room in subset:
                    embed.add_field(
                        name=room.room_name,
                        value=f'• Channel: {room.channel.mention}',
                        inline=False
                    )
                pages.append(embed)
        elif channel_obj:
            room = await TemporaryRoom.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not room:
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F {channel_obj.mention} is not a temporary room.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Temporary Room Info for {channel_obj.name}',
                color=discord.Color.blurple()
            )
            member = self.bot.get_user(room.member_snowflake)
            embed.add_field(name='Room Name', value=room.room_name, inline=False)
            embed.add_field(name='Owner', value=f'{member.display_name} ({member.mention})', inline=False)
            embed.add_field(name='Channel', value=channel_obj.mention, inline=False)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F Could not interpret the target. Provide "all", a channel, or a member.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')

    # DONE
    @commands.command(name='owners', help='Show temporary room stats for "all", a channel, or a member.')
    async def temp_room_stats_text_command(
        self,
        ctx,
        target: Optional[str] = commands.parameter(default=None, description='"all", a channel mention/ID, or a member mention/ID')
    ):
        state = State(ctx)
        channel_obj = None
        member_obj = None
        pages = []
        try:
            member_obj = await self.member_service.resolve_member(ctx, target)
        except:
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, target)
            except:
                channel_obj = ctx.channel
                await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F You are not authorized to list all temporary room owners in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            rooms = await TemporaryRoom.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not rooms:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F No temporary rooms exist.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            chunk = 12
            for i in range(0, len(rooms), chunk):
                subset = rooms[i:i+chunk]
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Temporary Rooms',
                    color=discord.Color.blurple()
                )
                for room in subset:
                    embed.add_field(
                        name=room.room_name,
                        value=f'• Owner: {room.member_mention}\n• Channel: {room.channel_mention}',
                        inline=False
                    )
                pages.append(embed)
        elif member_obj:
            rooms = await TemporaryRoom.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if not rooms:
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F {member_obj.display_name} does not own any temporary rooms in {ctx.guild.name}.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            chunk = 12
            for i in range(0, len(rooms), chunk):
                subset = rooms[i:i+chunk]
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Temporary Rooms Owned by {member_obj.display_name}',
                    color=discord.Color.blurple()
                )
                for room in subset:
                    embed.add_field(
                        name=room.room_name,
                        value=f'• Channel: {room.channel.mention}',
                        inline=False
                    )
                pages.append(embed)
        elif channel_obj:
            room = await TemporaryRoom.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not room:
                try:
                   return await state.end(warning=f'\U000026A0\U0000FE0F {channel_obj.mention} is not a temporary room.')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}.')
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Temporary Room Info for {channel_obj.name}',
                color=discord.Color.blurple()
            )
            member = self.bot.get_user(room.member_snowflake)
            embed.add_field(name='Room Name', value=room.room_name, inline=False)
            embed.add_field(name='Owner', value=f'{member.display_name} ({member.mention})', inline=False)
            embed.add_field(name='Channel', value=channel_obj.mention, inline=False)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F Could not interpret the target. Provide "all", a channel, or a member.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
       # DONE
    @app_commands.command(name='ping', description='Ping the bot!')
    async def ping_app_command(
        self,
        interaction: discord.Interaction
    ):
        state = State(interaction)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Pong!')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')

    # DONE
    @commands.command(name='ping', description='Ping the bot!')
    async def ping_text_command(
        self,
        ctx: commands.Context
    ):
        state = State(ctx)
        try:
            return await state.end(success=f"{self.emoji.get_random_emoji()} Pong!")
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')   

    # DONE
    @app_commands.command(name='roleid', description='Get the ID of a role by name in this server.')
    @app_commands.describe(role_name='The name of the role to look up')
    async def get_role_id_app_command(
        self,
        interaction: discord.Interaction,
        role_name: str
    ):
        state = State(interaction)
        role = discord.utils.get(interaction.guild.roles, name=role_name)
        if role:
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Role `{role.name}` has ID `{role.id}`.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No role named "{role_name}" found in this server.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
    # DONE
    @commands.command(name='roleid', help='Get the ID of a role by name in this server.')
    async def get_role_id(self, ctx: commands.Context, *, role_name: str):
        state = State(ctx)
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Role `{role.name}` has ID `{role.id}`.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No role named "{role_name}" found in this server.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
    # DONE
    @app_commands.command(name='survey', description='Survey moderators, developers, owners, and coordinators in the current or specified channel.')
    @app_commands.describe(channel='Tag a voice/stage channel')
    async def stage_survey_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
            await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        owners, developers, administrators, moderators, coordinators = [], [], [], [], []
        for member in channel_obj.members:
            match True:
                case _ if await member_is_owner(guild_snowflake=interaction.guild.id, member_snowflake=member.id):
                    owners.append(member)
                case _ if await member_is_developer(guild_snowflake=interaction.guild.id, member_snowflake=member.id):
                    developers.append(member)
                case _ if await member_is_administrator(guild_snowflake=interaction.guild.id, member_snowflake=member.id):
                    administrators.append(member)
                case _ if await member_is_coordinator(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member.id):
                    coordinators.append(member)
                case _ if await member_is_moderator(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member.id):
                    moderators.append(member)
        def fmt(users):
            return ', '.join(u.mention for u in users) if users else '*None*'
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} Survey results for {channel_obj.name}",
            description=f"Total surveyed: {len(channel_obj.members)}",
            color=discord.Color.blurple()
        )
        embed.add_field(name="Owners", value=fmt(owners), inline=False)
        embed.add_field(name="Developers", value=fmt(developers), inline=False)
        embed.add_field(name="Administrators", value=fmt(administrators), inline=False)
        embed.add_field(name="Coordinators", value=fmt(coordinators), inline=False)
        embed.add_field(name="Moderators", value=fmt(moderators), inline=False)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    # DONE
    @commands.command(name='survey', help='Survey moderators, developers, owners, and coordinators in the current or specified channel.')
    async def stage_survey_text_command(
        self,
        ctx: commands.Context,
        *,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
            await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        owners, developers, administrators, moderators, coordinators = [], [], [], [], []
        for member in channel_obj.members:
            match True:
                case _ if await member_is_owner(guild_snowflake=ctx.guild.id, member_snowflake=member.id):
                    owners.append(member)
                case _ if await member_is_developer(guild_snowflake=ctx.guild.id, member_snowflake=member.id):
                    developers.append(member)
                case _ if await member_is_administrator(guild_snowflake=ctx.guild.id, member_snowflake=member.id):
                    administrators.append(member)
                case _ if await member_is_coordinator(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member.id):
                    coordinators.append(member)
                case _ if await member_is_moderator(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member.id):
                    moderators.append(member)
        def fmt(users):
            return ', '.join(u.mention for u in users) if users else '*None*'
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} Survey results for {channel_obj.name}",
            description=f"Total surveyed: {len(channel_obj.members)}",
            color=discord.Color.blurple()
        )
        embed.add_field(name="Owners", value=fmt(owners), inline=False)
        embed.add_field(name="Developers", value=fmt(developers), inline=False)
        embed.add_field(name="Administrators", value=fmt(administrators), inline=False)
        embed.add_field(name="Coordinators", value=fmt(coordinators), inline=False)
        embed.add_field(name="Moderators", value=fmt(moderators), inline=False)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
        
async def setup(bot: DiscordBot):
    cog = EveryoneCommands(bot)
    await bot.add_cog(cog)

