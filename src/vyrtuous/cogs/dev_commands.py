''' dev_commands.py A discord.py cog containing developer commands for the Vyrtuous bot.

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
from typing import Literal, Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.check_service import *
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.role_service import RoleService
from vyrtuous.enhanced_members.administrator import Administrator
from vyrtuous.enhanced_members.administrator import AdministratorRole
from vyrtuous.utils.database import Database
from vyrtuous.utils.developer_log import DeveloperLog
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.properties.snowflake import *
from vyrtuous.service.state_service import State
from vyrtuous.rooms.temporary_room import TemporaryRoom


class DevCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.emoji = Emojis()
        self.message_service = MessageService(self.bot, self.bot.db_pool)
        self.member_service = MemberService()
        self.role_service = RoleService()
    
    # DONE
    @app_commands.command(name='backup', description='DB backup.')
    @developer_predicator()
    async def app_backup(
        self,
        interaction: discord.Interaction
    ):
        await interaction.response.defer(ephemeral=True)
        state = State(interaction)
        db = Database(directory='/app/backups')
        try:
            db.create_backup_directory()
            db.execute_backup()
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
        try:
            return await state.end(success=discord.File(db.file_name))
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='backup', help='DB backup.')
    @developer_predicator()
    async def text_backup(
        self,
        ctx: commands.Context
    ):
        state = State(ctx)
        db = Database(directory='/app/backups')
        try:
            db.create_backup_directory()
            db.execute_backup()
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
        try:
            return await state.end(success=discord.File(db.file_name))
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    @app_commands.command(name='cogs', description='Lists cogs.')
    @developer_predicator()
    async def list_cogs_app_command(self, interaction: discord.Interaction):
        state = State(interaction)
        loaded, not_loaded = [], []
        embed = discord.Embed(title=f'{self.emoji.get_random_emoji()} Cogs for {interaction.guild.me.name}', color=discord.Color.blurple())
        for cog in sorted(DISCORD_COGS_CLASSES):
            if cog in self.bot.cogs:
                loaded.append(cog)
            else:
                not_loaded.append(cog)
        if loaded:
            embed.add_field(name='Loaded', value='\n'.join(loaded), inline=False)
        if not_loaded:
            embed.add_field(name='Not Loaded', value='\n'.join(not_loaded), inline=False)
        if not loaded and not not_loaded:
            embed.add_field(name='No cogs available.', inline=False)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    @commands.command(name='cogs', help='Lists cogs.')
    @developer_predicator()
    async def list_cogs_text_command(self, ctx: commands.Context):
        state = State(ctx)
        loaded, not_loaded = [], []
        embed = discord.Embed(title=f'{self.emoji.get_random_emoji()} Cogs for {ctx.guild.me.name}', color=discord.Color.blurple())
        for cog in sorted(DISCORD_COGS_CLASSES):
            if cog in self.bot.cogs:
                loaded.append(cog)
            else:
                not_loaded.append(cog)
        if loaded:
            embed.add_field(name='Loaded', value='\n'.join(loaded), inline=False)
        if not_loaded:
            embed.add_field(name='Not Loaded', value='\n'.join(not_loaded), inline=False)
        if not loaded and not not_loaded:
            embed.add_field(name='No cogs available.', inline=False)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    @app_commands.command(name='dlog', description='Resolve or update the notes on an issue by reference.')
    @app_commands.describe(
        reference='Specify the issue reference ID.',
        action="'resolve' or 'append' or 'overwrite'.",
        notes='Optionally specify notes to append or overwrite.'
    )
    @developer_predicator()
    async def update_developer_logs_app_command(
        self,
        ctx: commands.Context,
        reference: str,
        action: str,
        notes: Optional[str] = None
    ):       
        state = State(ctx)
        developer_log = await DeveloperLog.fetch_unresolved_by_reference(id=reference)
        if not developer_log:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F Issue not found. Received: {reference}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        if action and action.lower() == 'resolve':
            await developer_log.resolve()
            detail = 'resolved the issue. The record will remain in the database for the next 30 days.'
        elif action and action.lower() == 'append':
            await developer_log.append(notes)
            detail = 'appended to the previous notes.'
        elif action and action.lower() == 'overwrite':
            await developer_log.overwrite(notes)
            detail = 'overwrote the previous notes.'
        try:
            return await state.end(success=f'\U000026A0\U0000FE0F You successfully {detail}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    @commands.command(name='dlog', help='Resolve or update the notes on an issue by reference')
    @developer_predicator()
    async def update_developer_logs_text_command(
        self,
        ctx: commands.Context,
        reference: str = commands.parameter(default=None, description='Specify the developer log reference ID.'),
        action: str = commands.parameter(default='append', description="Specify one of: 'resolve' or 'append' or 'overwrite'."),
        *,
        notes: Optional[str] = commands.parameter(default=None, description='Optionally specify notes.')
    ):       
        state = State(ctx)
        developer_log = await DeveloperLog.fetch_unresolved_by_reference(id=reference)
        if not developer_log:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F Issue not found. Received: {reference}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        if action and action.lower() == 'resolve':
            await developer_log.resolve()
            detail = 'resolved the issue'
        elif action and action.lower() == 'append':
            await developer_log.append(notes)
            detail = 'appended to the previous notes'
        elif action and action.lower() == 'overwrite':
            await developer_log.overwrite(notes)
            detail = 'overwrote the previous notes'
        try:
            return await state.end(success=f'\U000026A0\U0000FE0F You successfully {detail}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    
    @app_commands.command(name='dlogs', description="List issues.")
    @app_commands.describe(
        scope="Specify one of: 'all', 'resolved' or 'unresolved'",
        value='Specify one of: channel ID/mention, reference ID and server ID.'
    )
    @developer_predicator()
    async def list_developer_logs_text_command(
        self,
        interaction: discord.Interaction,
        scope: Optional[str],
        value: Optional[str]
    ):
        state = State(interaction)
        is_at_home = False
        channel_obj = None
        developer_logs = []
        guild_obj = None
        chunk_size = 7
        field_count = 0
        lines, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        skipped_message_snowflakes_by_guild_snowflake = {}
        title = f'{self.emoji.get_random_emoji()} Developer Logs'
        if scope and scope.lower() == 'all':
            developer_logs = await DeveloperLog.fetch_all()
        elif scope and scope.lower() == 'resolved':
            try:
                channel_obj = await self.channel_service.search(interaction, value) 
                developer_logs = await DeveloperLog.fetch_resolved_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                guild_obj = self.bot.get_guild(value)
                if guild_obj:
                    developer_logs = await DeveloperLog.fetch_resolved_by_guild(guild_snowflake=value)
                else:
                    developer_log = await DeveloperLog.fetch_resolved_by_reference(id=value)
                    if not developer_log:
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F Value must be one of: channel ID/mention, reference ID, server ID or empty. Received: {value}.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
        elif scope and scope.lower() == 'unresolved':
            try:
                channel_obj = await self.channel_service.search(interaction, value) 
                developer_logs = await DeveloperLog.fetch_unresolved_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            except Exception as e:
                guild_obj = self.bot.get_guild(value)
                if guild_obj:
                    developer_logs = await DeveloperLog.fetch_unresolved_by_guild(guild_snowflake=value)
                else:
                    developer_log = await DeveloperLog.fetch_unresolved_by_reference(id=value)
                    developer_logs = [developer_log]
                    if not developer_logs:
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F Value must be one of: channel ID/mention, reference ID, server ID or empty. Received: {value}.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')                        
        else:
            channel_obj = interaction.channel
            guild_obj = interaction.guild
        
        if not developer_logs:
            if scope:
                msg = f'No developer logs exist for scope: {scope}.'
            else:
                msg = f'No developer logs exist for {channel_obj.mention} in {guild_obj.name}.'
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        
        guild_dictionary = {}
        for developer_log in developer_logs:
            guild_dictionary.setdefault(developer_log.guild_snowflake, {})
            guild_dictionary[developer_log.guild_snowflake].setdefault(developer_log.channel_snowflake, [])
            guild_dictionary[developer_log.guild_snowflake][developer_log.channel_snowflake].append({
                'developer_snowflakes': developer_log.developer_snowflakes,
                'id': developer_log.id,
                'message_snowflake': developer_log.message_snowflake,
                'notes': developer_log.notes,
                'resolved': developer_log.resolved,
            })

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_logs in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                lines = []
                for member_data in channel_logs:
                    try:
                        msg = await channel.fetch_message(member_data['message_snowflake'])
                        lines.append(f'**Message:** {msg.jump_url}')
                    except Exception as e:
                        skipped_message_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_data['message_snowflake'])
                    if member_data['resolved'] == False:
                        resolved = '\u274C'
                    elif member_data['resolved'] == True:
                        resolved = '\u2705'
                    lines.append(f"{resolved}**Reference:** {member_data['id']}")
                    if value:
                        lines.append(f"**Notes:** {member_data['notes']}")
                    if member_data['developer_snowflakes']:
                        lines.append(f"**Assigned to:** {', '.join(member_data['developer_snowflakes'])}")
                    field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(name=f'**Channel:** {channel.mention}', value='\n'.join(lines), inline=False)
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
            if lines:
                embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(lines), inline=False)
            pages.append(embed)
        try:
            is_at_home = at_home(interaction_interaction_or_message=interaction)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guild_snowflakes:
                embed = discord.Embed(title='Skipped Servers', description='\u200b', color=discord.Color.blue())
                lines = []
                for guild_snowflake in skipped_guild_snowflakes:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_channel_snowflakes_by_guild_snowflake:
                for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                    embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in Server ({guild_snowflake})')
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in Server ({guild_snowflake}) continued...')
                            field_count = 0
                            lines = []
                        lines.append(str(channel_snowflake))
                        field_count += 1
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
            if skipped_message_snowflakes_by_guild_snowflake:
                for guild_snowflake, message_list in skipped_message_snowflakes_by_guild_snowflake.items():
                    embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Messages in Server ({guild_snowflake})')
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Messages in Server ({guild_snowflake}) continued...')
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
                    return await state.end(warning=f'\U000026A0\U0000FE0F Embed size is too large. Limit the scope.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No issues found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    @commands.command(name='dlogs', help="List issues.")
    @developer_predicator()
    async def list_developer_logs_text_command(
        self,
        ctx: commands.Context,
        *,
        scope: Optional[str] = commands.parameter(default=None, description="Specify one of: 'all', 'resolved' or 'unresolved."),
        value: Optional[str] = commands.parameter(default=None, description='Specify one of: channel ID/mention, reference ID or server ID.')
    ):
        state = State(ctx)
        is_at_home = False
        channel_obj = None
        developer_logs = []
        guild_obj = None
        chunk_size = 7
        field_count = 0
        lines, pages = [], []
        skipped_channel_snowflakes_by_guild_snowflake = {}
        skipped_guild_snowflakes = set()
        skipped_message_snowflakes_by_guild_snowflake = {}
        title = f'{self.emoji.get_random_emoji()} Developer Logs'
        if scope and scope.lower() == 'all':
            developer_logs = await DeveloperLog.fetch_all()
        elif scope and scope.lower() == 'resolved':
            try:
                channel_obj = await self.channel_service.search(ctx, value) 
                developer_logs = await DeveloperLog.fetch_resolved_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                guild_obj = self.bot.get_guild(value)
                if guild_obj:
                    developer_logs = await DeveloperLog.fetch_resolved_by_guild(guild_snowflake=value)
                else:
                    developer_log = await DeveloperLog.fetch_resolved_by_reference(id=value)
                    if not developer_log:
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F Value must be one of: channel ID/mention, reference ID, server ID or empty. Received: {value}.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')
        elif scope and scope.lower() == 'unresolved':
            try:
                channel_obj = await self.channel_service.search(ctx, value) 
                developer_logs = await DeveloperLog.fetch_unresolved_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            except Exception as e:
                guild_obj = self.bot.get_guild(value)
                if guild_obj:
                    developer_logs = await DeveloperLog.fetch_unresolved_by_guild(guild_snowflake=value)
                else:
                    developer_log = await DeveloperLog.fetch_unresolved_by_reference(id=value)
                    developer_logs = [developer_log]
                    if not developer_logs:
                        try:
                            return await state.end(warning=f'\U000026A0\U0000FE0F Value must be one of: channel ID/mention, reference ID, server ID or empty. Received: {value}.')
                        except Exception as e:
                            return await state.end(error=f'\u274C {str(e).capitalize()}')                        
        else:
            channel_obj = ctx.channel
            guild_obj = ctx.guild
        
        if not developer_logs:
            if scope:
                msg = f'No developer logs exist for scope: {scope}.'
            else:
                msg = f'No developer logs exist for {channel_obj.mention} in {guild_obj.name}.'
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}') 
        
        guild_dictionary = {}
        for developer_log in developer_logs:
            guild_dictionary.setdefault(developer_log.guild_snowflake, {})
            guild_dictionary[developer_log.guild_snowflake].setdefault(developer_log.channel_snowflake, [])
            guild_dictionary[developer_log.guild_snowflake][developer_log.channel_snowflake].append({
                'developer_snowflakes': developer_log.developer_snowflakes,
                'id': developer_log.id,
                'message_snowflake': developer_log.message_snowflake,
                'notes': developer_log.notes,
                'resolved': developer_log.resolved,
            })

        for guild_snowflake in guild_dictionary:
            guild_dictionary[guild_snowflake] = dict(sorted(guild_dictionary[guild_snowflake].items()))

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            if not guild:
                skipped_guild_snowflakes.add(guild_snowflake)
                continue
            embed = discord.Embed(title=title, description=guild.name, color=discord.Color.blue())
            for channel_snowflake, channel_logs in channels.items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    skipped_channel_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(channel_snowflake)
                    continue
                lines = []
                for member_data in channel_logs:
                    try:
                        msg = await channel.fetch_message(member_data['message_snowflake'])
                        lines.append(f'**Message:** {msg.jump_url}')
                    except Exception as e:
                        skipped_message_snowflakes_by_guild_snowflake.setdefault(guild_snowflake, []).append(member_data['message_snowflake'])
                    if member_data['resolved'] == False:
                        resolved = '\u274C'
                    elif member_data['resolved'] == True:
                        resolved = '\u2705'
                    lines.append(f"{resolved}**Reference:** {member_data['id']}")
                    if value:
                        lines.append(f"**Notes:** {member_data['notes']}")
                    if member_data['developer_snowflakes']:
                        lines.append(f"**Assigned to:** {', '.join(member_data['developer_snowflakes'])}")
                    field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(name=f'**Channel:** {channel.mention}', value='\n'.join(lines), inline=False)
                    pages.append(embed)
                    embed = discord.Embed(title=title, description=f'{guild.name} continued...', color=discord.Color.blue())
                    field_count = 0
            if lines:
                embed.add_field(name=f'Channel: {channel.mention}', value='\n'.join(lines), inline=False)
            pages.append(embed)
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            pass
        if is_at_home:
            if skipped_guild_snowflakes:
                embed = discord.Embed(title='Skipped Servers', description='\u200b', color=discord.Color.blue())
                lines = []
                for guild_snowflake in skipped_guild_snowflakes:
                    if field_count >= chunk_size:
                        embed.description = '\n'.join(lines)
                        pages.append(embed)
                        embed = discord.Embed(title='Skipped Servers continued...', color=discord.Color.red())
                        lines = []
                        field_count = 0
                    lines.append(str(guild_snowflake))
                    field_count += 1
                embed.description = '\n'.join(lines)
                pages.append(embed)
            if skipped_channel_snowflakes_by_guild_snowflake:
                for guild_snowflake, channel_list in skipped_channel_snowflakes_by_guild_snowflake.items():
                    embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in Server ({guild_snowflake})')
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Channels in Server ({guild_snowflake}) continued...')
                            field_count = 0
                            lines = []
                        lines.append(str(channel_snowflake))
                        field_count += 1
                    embed.description = '\n'.join(lines)
                    pages.append(embed)
            if skipped_message_snowflakes_by_guild_snowflake:
                for guild_snowflake, message_list in skipped_message_snowflakes_by_guild_snowflake.items():
                    embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Messages in Server ({guild_snowflake})')
                    field_count = 0
                    lines = []
                    for channel_snowflake in channel_list:
                        if field_count >= chunk_size:
                            embed.description = '\n'.join(lines)
                            pages.append(embed)
                            embed = discord.Embed(color=discord.Color.red(), title=f'Skipped Messages in Server ({guild_snowflake}) continued...')
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
                    return await state.end(warning=f'\U000026A0\U0000FE0F Embed size is too large. Limit the scope.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No issues found.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @app_commands.command(name='load', description="Loads a cog by name 'vyrtuous.cog.<cog_name>.'")
    @developer_predicator()
    async def load_app_command(self, interaction: discord.Interaction, module: str):
        await interaction.response.defer(ephemeral=True)
        state = State(interaction)
        try:
            await interaction.client.load_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully loaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE    
    @commands.command(name='load', help="Loads a cog by name 'vyrtuous.cog.<cog_name>.'")
    @developer_predicator()
    async def load_text_command(self, ctx: commands.Context, *, module: str):
        state = State(ctx)
        try:
            await self.bot.load_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully loaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @app_commands.command(name='ping', description='Ping me!')
    @developer_predicator()
    async def ping_app_command(
        self,
        interaction: discord.Interaction
    ):
        state = State(interaction)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Pong!')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @commands.command(name='ping', help='Ping me!')
    @developer_predicator()
    async def ping_text_command(
        self,
        ctx: commands.Context
    ):
        state = State(ctx)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Pong!')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
               
    # DONE
    @app_commands.command(name='reload', description="Reloads a cog by name 'vyrtuous.cog.<cog_name>'.")
    @app_commands.check(at_home)
    @developer_predicator()
    async def reload_app_command(self, interaction: discord.Interaction, module: str):
        await interaction.response.defer(ephemeral=True)
        state = State(interaction)
        try:
            await interaction.client.reload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully reloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(name='reload', help="Reloads a cog by name 'vyrtuous.cog.<cog_name>'.")
    @developer_predicator()
    async def reload_text_command(self, ctx: commands.Context, *, module: str):
        state = State(ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully reloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(name='sync', description="Sync app commands.")
    @developer_predicator()
    async def sync_app_command(
        self,
        interaction: discord.Interaction,
        spec: Optional[Literal['~', '*', '^']] = None
    ):
        await interaction.response.defer(ephemeral=True)
        state = State(interaction)
        guilds = interaction.client.guilds
        synced = []
        if not guilds:
            if spec == '~':
                synced = await interaction.client.tree.sync(guild=interaction.guild)
            elif spec == '*':
                interaction.client.tree.copy_global_to(guild=interaction.guild)
                synced = await interaction.client.tree.sync(guild=interaction.guild)
            elif spec == '^':
                interaction.client.tree.clear_commands(guild=interaction.guild)
                await interaction.client.tree.sync(guild=interaction.guild)
            else:
                synced = await interaction.client.tree.sync()
            try:
                if spec is None:
                    msg = f'Synced {len(synced)} commands globally.'
                else:
                    msg = f'Synced {len(synced)} commands to the current server.'
                return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
            except Exception as e:   
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        ret = 0
        for guild in guilds:
            try:
                await interaction.client.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Synced the tree to {ret}/{len(guilds)}.')
        except Exception as e:   
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='sync', help="Sync app commands.")
    @developer_predicator()
    async def sync_text_command(
        self,
        ctx: commands.Context,
        guilds: commands.Greedy[discord.Object],
        spec: Optional[Literal['~', '*', '^']] = None
    ):
        state = State(ctx)
        synced = []
        if not guilds:
            if spec == '~':
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '*':
                ctx.bot.tree.copy_global_to(guild=ctx.guild)
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '^':
                ctx.bot.tree.clear_commands(guild=ctx.guild)
                await ctx.bot.tree.sync(guild=ctx.guild)
            else:
                synced = await ctx.bot.tree.sync()
            try:
                if spec is None:
                    msg = f'Synced {len(synced)} commands globally.'
                else:
                    msg = f'Synced {len(synced)} commands to the current server.'
                return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Synced the tree to {ret}/{len(guilds)}.')
        except Exception as e:   
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(name='unload', description="Unloads a cog by name 'vyrtuous.cog.<cog_name>'.")
    @developer_predicator()
    async def unload_app_command(self, interaction: discord.Interaction, module: str):
        await interaction.response.defer(ephemeral=True)
        state = State(interaction)
        try:
            await interaction.client.unload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully unloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @commands.command(name='unload', help="Unloads a cog by name 'vyrtuous.cog.<cog_name>'.")
    @developer_predicator()
    async def unload_text_command(self, ctx: commands.Context, *, module: str):
        state = State(ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully unloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}') 

async def setup(bot: DiscordBot):
    await bot.add_cog(DevCommands(bot))
