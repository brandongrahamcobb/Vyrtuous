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
        interaction: discord.Interaction
    ):
        administrators = await Administrator.fetch_members_by_guild(guild_snowflake=interaction.guild.id)
        pages = await All.create_pages_to_show_members(members=administrators, member_type=Administrator)
        if pages:
            paginator = Paginator(self.bot, interaction, pages)
            return await paginator.start()
        else:
            return await self.handler.send_message(interaction, content=f'{self.emoji.get_random_emoji()} No admins found in {interaction.guild.name}.')
            
    # DONE
    @commands.command(name='admins', help='Lists all members with server mute privileges in this guild.')
    async def list_administrators_text_command(
        self,
        ctx: commands.Context
    ) -> None:
        administrators = await Administrator.fetch_members_by_guild(guild_snowflake=ctx.guild.id)
        pages = await All.create_pages_to_show_members(members=administrators, member_type=Administrator)
        if pages:
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        else:
            return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No admins found in {ctx.guild.name}.')
    
    # DONE
    @app_commands.command(name='coords', description='Lists coordinators for a specific voice channel, all, or a member.')
    @app_commands.describe(target='"all", member or channel name/ID/mention')
    async def list_coordinators_app_command(
        self,
        interaction : discord.Interaction,
        target : Optional[str] = None
    ):
        member_obj = await self.member_service.resolve_member(interaction, target)
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(interaction, content='\U0001F6AB You are not authorized to list all coordinators.')
            coordinators = await Coordinator.fetch_members_by_guild(guild_snowflake=interaction.guild.id)
            pages = await All.create_pages_to_show_guilds_by_members(members=coordinators, member_type=Coordinator)
        elif member_obj:
            channel_snowflakes = await Coordinator.fetch_channels_by_guild_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_channels(channel_snowflakes=channel_snowflakes, member_type=Coordinator)
        elif channel_obj:
            coordinators = await Coordinator.fetch_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            pages = await All.create_pages_to_show_members(members=coordinators, member_type=Coordinator)
        if pages:
            paginator = Paginator(self.bot, interaction, pages)
            return await paginator.start()
        return await self.handler.send_message(interaction, content=f'{self.emoji.get_random_emoji()} No coordinators found.')

    # DONE
    @commands.command(name='coords', help='Lists coordinators for a specific voice channel, all, or a member.')
    async def list_coordinators_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Voice channel name, mention, ID, "all", or member ID')
    ) -> None:
        member_obj = await self.member_service.resolve_member(ctx, target)
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all coordinators.')
            coordinators = await Coordinator.fetch_members_by_guild(guild_snowflake=ctx.guild.id)
            pages = await All.create_pages_to_show_guilds_by_members(members=coordinators, member_type=Coordinator)
        elif member_obj:
            channel_snowflakes = await Coordinator.fetch_channels_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_channels(channel_snowflakes=channel_snowflakes, member_type=Coordinator)
        elif channel_obj:
            coordinators = await Coordinator.fetch_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            pages = await All.create_pages_to_show_members(members=coordinators, member_type=Coordinator)
        if pages:
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No coordinators found.')
            
    
    # DONE
    @app_commands.command(name='devs', description='Lists developers.')
    @app_commands.describe(target='"all" or a member snowflake ID/mention')
    async def list_developers_app_command(
        self,
        interaction : discord.Interaction,
        target : Optional[str] = None
    ):
        member_obj = await self.member_service.resolve_member(interaction, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                return await self.handler.send_message(interaction, content='\U0001F6AB You are not authorized to list all developers.')
            developers = await Developer.fetch_all()
            pages = await All.create_pages_to_show_guilds_by_members(members=developers, member_type=Developer)
        elif member_obj:
            guilds = []
            developers = await Developer.fetch_by_member(member_snowflake=member_obj.id)
            for developer in developers:
                guilds.append(developer.guild_snowflake)
            pages = await All.create_pages_to_show_members_by_guild(guilds=guilds, member_snowflake=member_obj.id, member_type=Developer)
        else:
            developers = await Developer.fetch_by_guild(guild_snowflake=interaction.guild.id)
            pages = await All.create_pages_to_show_members(members=developers, member_type=Developer)
        if pages:
            paginator = Paginator(self.bot, interaction, pages)
            return await paginator.start()
        return await self.handler.send_message(interaction, content=f'{self.emoji.get_random_emoji()} No developers found.')
        
    # DONE
    @commands.command(name='devs', help='Lists developers.')
    async def list_developers_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", or user mention/ID')
    ) -> None:
        member_obj = await self.member_service.resolve_member(ctx, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all developers.')
            developers = await Developer.fetch_all()
            pages = await All.create_pages_to_show_guilds_by_members(members=developers, member_type=Developer)
        elif member_obj:
            guilds = []
            developers = await Developer.fetch_by_member(member_snowflake=member_obj.id)
            for developer in developers:
                guilds.append(developer['guild_snowflake'])
            pages = await All.create_pages_to_show_guilds_by_member(guilds=guilds, member_snowflake=member_obj.id, member_type=Developer)
        else:
            developers = await Developer.fetch_by_guild(guild_snowflake=ctx.guild.id)
            pages = await All.create_pages_to_show_members(members=developers, member_type=Developer)
        if pages:
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No developers found.')
        
    # DONE
    @app_commands.command(name='mods', description='Lists moderator statistics.')
    @app_commands.describe(target='A member or channel snowflake ID/mention')
    async def list_moderators_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        member_obj = await self.member_service.resolve_member(interaction, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(interaction, content='\U0001F6AB You are not authorized to list all moderators.')
            moderators = await Moderator.fetch_members_by_guild(guild_snowflake=interaction.guild.id)
            pages = await All.create_pages_to_show_members(members=moderators, member_type=Moderator)
        elif member_obj:
            channel_snowflakes = await Moderator.fetch_channels_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_channels(channel_snowflakes=channel_snowflakes, member_type=Moderator)
        elif channel_obj:
            moderators = await Moderator.fetch_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            pages = await All.create_pages_to_show_members(members=moderators, member_type=Moderator)
        if pages:
            paginator = Paginator(self.bot, interaction, pages)
            return await paginator.start()
        return await self.handler.send_message(interaction, content=f'{self.emoji.get_random_emoji()} No moderators found for {channel_obj.mention}.')
        
    # DONE
    @commands.command(name='mods',help='Lists moderator statistics.')
    async def list_moderators_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Voice channel name/mention/ID, "all", or member mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all moderators.')
            moderators = await Moderator.fetch_members_by_guild(guild_snowflake=ctx.guild.id)
            pages = await All.create_pages_to_show_members(members=moderators, member_type=Moderator)
        elif member_obj:
            channel_snowflakes = await Moderator.fetch_channels_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            pages = await All.create_pages_to_show_channels(channel_snowflakes=channel_snowflakes, member_type=Moderator)
        elif channel_obj:
            moderators = await Moderator.fetch_members_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            pages = await All.create_pages_to_show_members(members=moderators, member_type=Moderator)
        if pages:
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No moderators found.')

    # DONE
    @app_commands.command(name='owners', description='Show temporary room stats for "all", a channel, or a member.')
    @app_commands.describe(target='"all", a channel mention/ID, or a member mention/ID')
    async def temp_room_stats_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await self.handler.send_message(interaction, content='\U0001F6AB This command can only be used in servers.')
        async with self.bot.db_pool.acquire() as conn:
            channel_obj = await self.channel_service.resolve_channel(interaction, target)
            member_obj = await self.member_service.resolve_member(interaction, target)
            highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
            if target and target.lower() == 'all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await self.handler.send_message(interaction, content='\U0001F6AB You are not authorized to list all owners.')
                rooms = await TemporaryRoom.fetch_by_guild(guild_snowflake=interaction.guild.id)
                if not rooms:
                    return await self.handler.send_message(interaction, content=f'{self.emoji.get_random_emoji()} No temporary rooms exist.')
                pages = []
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
                paginator = Paginator(self.bot, interaction, pages)
                return await paginator.start()
            if member_obj:
                rooms = await TemporaryRoom.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
                if not rooms:
                    return await self.handler.send_message(interaction, content=f'{self.emoji.get_random_emoji()}{member_obj.display_name} does not own any temporary rooms.')
                pages = []
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
                paginator = Paginator(self.bot, interaction, pages)
                return await paginator.start()
            if channel_obj:
                room = await TemporaryRoom.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
                if not room:
                    return await self.handler.send_message(interaction, content=f'{self.emoji.get_random_emoji()} {channel_obj.mention} is not a temporary room.')
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Temporary Room Info for {channel_obj.name}',
                    color=discord.Color.blurple()
                )
                member = self.bot.get_user(room.member_snowflake)
                embed.add_field(name='Room Name', value=room.room_name, inline=False)
                embed.add_field(name='Owner', value=f'{member.display_name} ({member.mention})', inline=False)
                embed.add_field(name='Channel', value=channel_obj.mention, inline=False)
                return await self.handler.send_message(interaction, embed=embed)
            return await self.handler.send_message(interaction, content='\U0001F6AB Could not interpret the target. Provide "all", a channel, or a member.')

    # DONE
    @commands.command(name='owners', help='Show temporary room stats for "all", a channel, or a member.')
    async def temp_room_stats_text_command(
        self,
        ctx,
        target: Optional[str] = commands.parameter(default=None, description='"all", a channel mention/ID, or a member mention/ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        async with self.bot.db_pool.acquire() as conn:
            channel_obj = await self.channel_service.resolve_channel(ctx, target)
            member_obj = await self.member_service.resolve_member(ctx, target)
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
            def make_pages(rooms):
                pages = []
                chunk = 12
                for i in range(0, len(rooms), chunk):
                    subset = rooms[i:i+chunk]
                    embed = discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Temporary Rooms',
                        color=discord.Color.blurple()
                    )
                    member = self.bot.get_user(room.member_snowflake)
                    for room in subset:
                        embed.add_field(
                            name=room.room_name,
                            value=f'• Owner: {member.display_name}\n• Channel: {room.channel_mention}',
                            inline=False
                        )
                    pages.append(embed)
                return pages
            if target and target.lower() == 'all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB You are not authorized to list all owners.')
                rooms = await TemporaryRoom.fetch_by_guild(guild_snowflake=ctx.guild.id)
                if not rooms:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No temporary rooms exist.')
                pages = make_pages(rooms)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            if member_obj:
                rooms = await TemporaryRoom.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                if not rooms:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.display_name} does not own any temporary rooms.')
                pages = make_pages(rooms)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            if channel_obj:
                room = await TemporaryRoom.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
                if not room:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {channel_obj.mention} is not a temporary room.')
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Temporary Room Info for {channel_obj.name}',
                    color=discord.Color.blurple()
                )
                embed.add_field(name='Room Name', value=room.room_name, inline=False)
                embed.add_field(name='Owner', value=f"'{room.member_mention}", inline=False)
                embed.add_field(name='Channel', value=channel_obj.mention, inline=False)
                return await self.handler.send_message(ctx, embed=embed)
            return await self.handler.send_message(ctx, content='\U0001F6AB Could not interpret the target. Provide "all", a channel, or a member.')
    
       # DONE
    @app_commands.command(name='ping', description='Ping the bot!')
    async def ping_app_command(
        self,
        interaction: discord.Interaction
    ):
        state = State(interaction)
        return await state.end(success=f'{self.emoji.get_random_emoji()} Pong!')

    # DONE
    @commands.command(name='ping', description='Ping the bot!')
    async def ping_text_command(
        self,
        ctx: commands.Context
    ):
        state = State(ctx)
        return await state.end(success=f"{self.emoji.get_random_emoji()} Pong!")   

    # DONE
    @app_commands.command(name='roleid', description='Get the ID of a role by name in this server.')
    @app_commands.describe(role_name='The name of the role to look up')
    async def get_role_id_app_command(
        self,
        interaction: discord.Interaction,
        role_name: str
    ):
        role = discord.utils.get(interaction.guild.roles, name=role_name)
        if role:
            await self.handler.send_message(interaction, content=f'{self.emoji.get_random_emoji()} Role `{role.name}` has ID `{role.id}`.')
        else:
            await self.handler.send_message(interaction, content=f'{self.emoji.get_random_emoji()} No role named "{role_name}" found in this server.')

    # DONE
    @commands.command(name='roleid', help='Get the ID of a role by name in this server.')
    async def get_role_id(self, ctx: commands.Context, *, role_name: str):
        role = discord.utils.get(ctx.guild.roles, name=role_name)
        if role:
            await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Role `{role.name}` has ID `{role.id}`.')
        else:
            await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No role named "{role_name}" found in this server.')
    
    # DONE
    @app_commands.command(name='survey', description='Survey moderators, developers, owners, and coordinators in the current or specified channel.')
    @app_commands.describe(channel='Tag a voice/stage channel')
    async def stage_survey_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake
    ):
        if not interaction.guild:
            return await self.handler.send_message(interaction, content='This command must be used in a server.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if not isinstance(channel_obj, discord.VoiceChannel):
            return await self.handler.send_message(interaction, content='\U0001F6AB Please specify a valid target.')
        owners, developers, administrators, moderators, coordinators = [], [], [], [], []
        for member in channel_obj.members:
            match True:
                case _ if await member_is_owner(member):
                    owners.append(member)
                case _ if await member_is_developer(member):
                    developers.append(member)
                case _ if await member_is_administrator(member):
                    administrators.append(member)
                case _ if await member_is_coordinator(channel_obj, member):
                    coordinators.append(member)
                case _ if await member_is_moderator(channel_obj, member):
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
        await self.handler.send_message(interaction, embed=embed)
    
    # DONE
    @commands.command(name='survey', help='Survey moderators, developers, owners, and coordinators in the current or specified channel.')
    async def stage_survey_text_command(
        self,
        ctx: commands.Context,
        *,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if not isinstance(channel_obj, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid voice or stage channel.')
        owners, developers, administrators, moderators, coordinators = [], [], [], [], []
        for member in channel_obj.members:
            match True:
                case _ if await member_is_owner(member):
                    owners.append(member)
                case _ if await member_is_developer(member):
                    developers.append(member)
                case _ if await member_is_administrator(member):
                    administrators.append(member)
                case _ if await member_is_coordinator(channel_obj, member):
                    coordinators.append(member)
                case _ if await member_is_moderator(channel_obj, member):
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
        await self.handler.send_message(ctx, embed=embed)
        
async def setup(bot: DiscordBot):
    cog = EveryoneCommands(bot)
    await bot.add_cog(cog)

