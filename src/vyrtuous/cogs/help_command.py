''' help_command.py

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
import datetime
import inspect
import os
import re
import types
from typing import Any, Coroutine, Optional

from discord.ext.commands import Command

from vyrtuous.inc.helpers import *
from vyrtuous.service.check_service import *
from vyrtuous.service.discord_message_service import DiscordMessageService, Paginator, UserPaginator
from vyrtuous.utils.setup_logging import logger
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.permission import Permission

class Help(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.aliases_cog = bot.get_cog("Aliases")
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
        self.permission_page_title_pairs = [
            ('Owner', '`Owner` inherits `developer`.'),
            ('Developer', '`Developer` inherits `administrator`.'),
            ('Administrator', '`Administrator` inherits `coordinator`.'),
            ('Coordinator', '`Coordinator` inherits `moderator`.'),
            ('Moderator', 'Moderators can use these commands.'),
            ('Everyone', 'Commands available to everyone.')
        ]
    
    async def get_alias_help(self, channel: discord.abc.GuildChannel):
        contextual_aliases = []
        for alias_type, help_lines in self.aliases_cog.alias_help.items():
            for line in help_lines:
                contextual_aliases.append(f"**{alias_type}** – {line}")
        return contextual_aliases
        
    async def get_available_commands(self, bot, ctx_or_interaction) -> list[commands.Command]:
        available = []
        user_highest = await is_owner_developer_administrator_coordinator_moderator(ctx_or_interaction)
        for command in bot.commands:
            try:
                perm_level = await self.get_command_permission_level(bot, command)
                if Permission.PERMISSION_TYPES.index(user_highest) <= Permission.PERMISSION_TYPES.index(perm_level):
                    available.append(command)
            except Exception as e:
                logger.warning(f'\U0001F6AB Exception while evaluating command {command}: {e}')
        return available
    
    async def get_command_permission_level(self, bot, command):
        if not hasattr(command,'checks') or not command.checks:
            return 'Everyone'
        for check in command.checks:
            if hasattr(check,'__wrapped__'):
                func = check.__wrapped__
            else:
                func = check
            if hasattr(func,'_permission_level'):
                return func._permission_level
        return 'Everyone'

    def get_permission_color(self, perm_level):
        colors = {
            'Owner': discord.Color.red(),
            'Developer': discord.Color.purple(),
            'Administrator': discord.Color.blurple(),
            'Coordinator': discord.Color.orange(),
            'Moderator': discord.Color.blue(),
            'Everyone': discord.Color.green()
        }
        return colors.get(perm_level, discord.Color.greyple())
    
    async def group_commands_by_permission(self, bot, ctx_or_interaction, commands_list):
        permission_groups = {level:[] for level in Permission.PERMISSION_TYPES}
        for command in commands_list:
            perm_level = await self.get_command_permission_level(bot, command)
            if perm_level in permission_groups:
                permission_groups[perm_level].append(command)
            else:
                permission_groups['Everyone'].append(command)
        return permission_groups


    def split_command_list(self, commands_list, max_length=1024):
        chunks = []
        current_chunk = []
        current_length = 0
        for cmd in commands_list:
            cmd_line = f'**{self.config['discord_command_prefix']}{cmd.name}** – {cmd.help or 'No description'}'
            cmd_length = len(cmd_line)
            if current_length + cmd_length > max_length and current_chunk:
                chunks.append('\n'.join(current_chunk))
                current_chunk = [cmd_line.rstrip()]
                current_length = cmd_length
            else:
                current_chunk.append(cmd_line.rstrip())
                current_length += cmd_length
        if current_chunk:
            chunks.append('\n'.join(current_chunk))
        return chunks
    
    @app_commands.command(name='help', description='Show command information or your available commands.')
    @app_commands.describe(command_name='The command to view details for.')
    async def help_app_command(self, interaction: discord.Interaction, command_name: Optional[str] = None):
        bot = interaction.client

        if command_name:
            # Detailed view for a specific command
            cmd = bot.get_command(command_name.lower())
            if not cmd:
                return await interaction.response.send_message(f'\U0001F6AB Command `{command_name}` not found.', ephemeral=True)

            embed = discord.Embed(
                title=f'{self.bot.config["discord_command_prefix"]}{cmd.name}',
                description=cmd.help or 'No description provided.',
                color=discord.Color.blue()
            )

            sig = inspect.signature(cmd.callback)
            parameters = list(sig.parameters.items())
            if parameters and parameters[0][0]=='self':
                parameters.pop(0)
            if parameters and parameters[0][0]=='ctx':
                parameters.pop(0)

            if parameters:
                usage_parts = [f'{self.bot.config["discord_command_prefix"]}{cmd.name}']
                param_details = []
                for name, param in parameters:
                    is_optional = param.kind == inspect.Parameter.KEYWORD_ONLY
                    default = param.default
                    description = getattr(default, 'description', None) if isinstance(default, commands.Parameter) else None
                    annotation = param.annotation.__name__ if hasattr(param.annotation, '__name__') else str(param.annotation)
                    usage_parts.append(f'[{name}]' if is_optional else f'<{name}>')
                    detail = f'**{name}** ({annotation})'
                    if description:
                        detail += f': {description}'
                    param_details.append(detail)

                embed.add_field(name='Usage', value=f'`{" ".join(usage_parts)}`', inline=False)
                if param_details:
                    embed.add_field(name='Parameter Details', value='\n'.join(param_details), inline=False)

            # Include alias help only for this command
            alias_help_lines = await self.get_alias_help_for_command(cmd.name, interaction.channel)
            if alias_help_lines:
                embed.add_field(name='Aliases', value=''.join(alias_help_lines), inline=False)

            return await interaction.response.send_message(embed=embed, ephemeral=True)

        # General help: show only commands grouped by permission
        all_commands = await self.get_available_commands(bot, interaction)
        permission_groups = await self.group_commands_by_permission(bot, interaction, all_commands)
        pages = []

        user_highest = await is_owner_developer_administrator_coordinator_moderator(interaction)
        user_index = Permission.PERMISSION_TYPES.index(user_highest)

        for perm_level, description in self.permission_page_title_pairs:
            if Permission.PERMISSION_TYPES.index(perm_level) < user_index:
                continue
            commands_in_level = sorted(permission_groups.get(perm_level, []), key=lambda c: c.name)
            if not commands_in_level:
                continue
            embed = discord.Embed(title=f"{perm_level} Commands", description=description, color=self.get_permission_color(perm_level))
            command_lines = []
            for cmd in commands_in_level:
                cmd_name = cmd.name
                cmd_help = cmd.help or "No description"
                command_lines.append(f'**{self.config["discord_command_prefix"]}{cmd_name}** – {cmd_help}')
            
            command_text = '\n'.join(command_lines)
            if len(command_text) > 1024:
                chunks = self.split_command_list(commands_in_level)
                for i, chunk in enumerate(chunks):
                    field_name = f'{perm_level} Commands' if i == 0 else f'{perm_level} Commands (cont.)'
                    embed.add_field(name='', value=chunk, inline=False)
            else:
                embed.add_field(name='', value=command_text, inline=False)
            pages.append(embed)

        if not pages:
            return await interaction.response.send_message('\U0001F6AB No commands available to you.', ephemeral=True)

        paginator = UserPaginator(bot, interaction, pages)
        await paginator.start()
        
    @commands.command(name='help')
    async def help_text_command(self, ctx, *, command_name: str = None):
        bot = ctx.bot
        if command_name:
            # Detailed view for a single command
            cmd = bot.get_command(command_name.lower())
            if not cmd:
                return await self.handler.send_message(ctx, f'\U0001F6AB Command `{command_name}` not found.')
            embed = discord.Embed(
                title=f'{self.config["discord_command_prefix"]}{cmd.name}',
                description=cmd.help or 'No description provided.',
                color=discord.Color.blue()
            )
    
            # Show parameters for the command
            sig = inspect.signature(cmd.callback)
            parameters = list(sig.parameters.items())
            if parameters and parameters[0][0] == 'self':
                parameters.pop(0)
            if parameters and parameters[0][0] == 'ctx':
                parameters.pop(0)
            if parameters:
                usage_parts = [f'{self.config["discord_command_prefix"]}{cmd.name}']
                param_details = []
                for name, param in parameters:
                    is_optional = param.kind == inspect.Parameter.KEYWORD_ONLY
                    default = param.default
                    description = getattr(default, 'description', None) if isinstance(default, commands.Parameter) else None
                    annotation = param.annotation.__name__ if hasattr(param.annotation, '__name__') else str(param.annotation)
                    usage_parts.append(f'[{name}]' if is_optional else f'<{name}>')
                    detail = f'**{name}** ({annotation})'
                    if description:
                        detail += f': {description}'
                    param_details.append(detail)
                embed.add_field(name='Usage', value=f'`{" ".join(usage_parts)}`', inline=False)
                if param_details:
                    embed.add_field(name='Parameter Details', value='\n'.join(param_details), inline=False)
    
            # Include alias help for this command
            alias_help_lines = []
            for alias in self.aliases_cog.get_all_aliases():
                if alias.alias_name and alias.alias_name.lower() == cmd.name.lower():
                    for line in alias.help_lines:
                        alias_help_lines.append(f'**{alias.alias_type}** – {line}')
            if alias_help_lines:
                embed.add_field(name='Aliases', value=''.join(alias_help_lines), inline=False)
    
            return await self.handler.send_message(ctx, embed=embed)
        all_commands = await self.get_available_commands(bot, ctx)
        permission_groups = await self.group_commands_by_permission(bot, ctx, all_commands)
        pages = []
        user_highest = await is_owner_developer_administrator_coordinator_moderator(ctx)
        user_index = Permission.PERMISSION_TYPES.index(user_highest)
    
        for perm_level, description in self.permission_page_title_pairs:
            if Permission.PERMISSION_TYPES.index(perm_level) < user_index:
                continue
            commands_in_level = sorted(permission_groups.get(perm_level, []), key=lambda c: c.name)
            if not commands_in_level:
                continue
            embed = discord.Embed(title=f"{perm_level} Commands", description=description, color=self.get_permission_color(perm_level))
            command_lines = []
            for cmd in commands_in_level:
                cmd_name = cmd.name
                cmd_help = cmd.help or "No description"
                command_lines.append(f'**{self.config["discord_command_prefix"]}{cmd_name}** – {cmd_help}')
            
            command_text = '\n'.join(command_lines)
            if len(command_text) > 1024:
                chunks = self.split_command_list(commands_in_level)
                for i, chunk in enumerate(chunks):
                    field_name = f'{perm_level} Commands' if i == 0 else f'{perm_level} Commands (cont.)'
                    embed.add_field(name='', value=chunk, inline=False)
            else:
                embed.add_field(name='', value=command_text, inline=False)
            pages.append(embed)
    
        if not pages:
            return await self.handler.send_message(ctx, '\U0001F6AB No commands available to you.')
        paginator = Paginator(bot, ctx, pages)
        await paginator.start()

async def setup(bot: DiscordBot):
    cog = Help(bot)
    await bot.add_cog(cog)
