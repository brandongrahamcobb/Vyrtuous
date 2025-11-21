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
from typing import Any, Coroutine, Optional

from discord.ext.commands import Command

from vyrtuous.inc.helpers import *
from vyrtuous.service.check_service import *
from vyrtuous.service.discord_message_service import DiscordMessageService, Paginator
from vyrtuous.utils.setup_logging import logger
from vyrtuous.bot.discord_bot import DiscordBot

PERMISSION_ORDER = ['Owner', 'Developer', 'Coordinator', 'Moderator', 'Everyone']

class Help(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)

    async def get_available_commands(self, bot, ctx) -> list[commands.Command]:
        available_commands = []
        for command in bot.commands:
            try:
                if await command.can_run(ctx):
                    available_commands.append(command)
            except commands.CheckFailure:
                continue
            except Exception as e:
                logger.warning(f'\U0001F6AB Exception while checking command \'{command}\': {e}')
        return available_commands

    async def get_command_permission_level(self, bot, ctx, command):
        if not hasattr(command, 'checks') or not command.checks:
            return 'Everyone'
        for check in command.checks:
            func = check
            if hasattr(func, '__wrapped__'):
                func = func.__wrapped__
            if hasattr(func, '_permission_level'):
                return func._permission_level
        return 'Everyone'
        
    def get_permission_color(self, perm_level):
        colors = {
            'Owner': discord.Color.red(),
            'Developer': discord.Color.purple(),
            'Coordinator': discord.Color.orange(),
            'Moderator': discord.Color.blue(),
            'Everyone': discord.Color.green()
        }
        return colors.get(perm_level, discord.Color.greyple())
    
    async def get_user_highest_permission(self, bot, ctx):
        try:
            if await is_system_owner(ctx): return 'Owner'
        except commands.CheckFailure: pass
        try:
            if await is_guild_owner(ctx): return 'Owner'
        except commands.CheckFailure: pass
        try:
            if await is_developer(ctx): return 'Developer'
        except commands.CheckFailure: pass
        if ctx.guild:
            room_name = getattr(ctx.channel, 'name', None)
            try:
                async with ctx.bot.db_pool.acquire() as conn:
                    user_row = await conn.fetchrow(
                        'SELECT coordinator_channel_ids, moderator_channel_ids, coordinator_room_names, moderator_room_names '
                        'FROM users WHERE discord_snowflake=$1',
                        ctx.author.id
                    )
                    if user_row:
                        c_chan = user_row.get('coordinator_channel_ids') or []
                        m_chan = user_row.get('moderator_channel_ids') or []
                        c_room = user_row.get('coordinator_room_names') or []
                        m_room = user_row.get('moderator_room_names') or []
                        if c_chan or c_room: return 'Coordinator'
                        if m_chan or m_room: return 'Moderator'
            except Exception as e:
                logger.warning(f'Error checking coordinator/moderator permissions: {e}')
        return 'Everyone'
    
    async def group_commands_by_permission(self, bot, ctx, commands_list):
        permission_groups = {level: [] for level in PERMISSION_ORDER}
        for command in commands_list:
            perm_level = await self.get_command_permission_level(bot, ctx, command)
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
            cmd_line = f'**{config['discord_command_prefix']}{cmd.name}** – {cmd.help or 'No description'}\n'
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
        
    @commands.command(name='help')
    async def help(self, ctx, *, command_name: str = None):
        bot = ctx.bot
        if command_name:
            cmd = bot.get_command(command_name.lower())
            if not cmd:
                return await self.handler.send_message(ctx, f'\U0001F6AB Command `{command_name}` not found.')
            if cmd.hidden:
                return await self.handler.send_message(ctx, f'\U0001F6AB Command `{command_name}` is hidden.')
            embed = discord.Embed(
                title=f'{config['discord_command_prefix']}{cmd.name}',
                description=cmd.help or 'No description provided.',
                color=discord.Color.blue()
            )
            sig = inspect.signature(cmd.callback)
            parameters = list(sig.parameters.items())
            if parameters and parameters[0][0] == 'self':
                parameters.pop(0)
            if parameters and parameters[0][0] == 'ctx':
                parameters.pop(0)
            if parameters:
                usage_parts = [f'{config['discord_command_prefix']}{cmd.name}']
                param_details = []
                for name, param in parameters:
                    is_optional = param.kind == inspect.Parameter.KEYWORD_ONLY
                    default = param.default
                    description = (
                        getattr(default, 'description', None)
                        if isinstance(default, commands.Parameter)
                        else None
                    )
                    annotation = (
                        param.annotation.__name__
                        if hasattr(param.annotation, '__name__')
                        else str(param.annotation)
                    )
                    if is_optional:
                        usage_parts.append(f'[{name}]')
                    else:
                        usage_parts.append(f'<{name}>')
                    detail = f'**{name}** ({annotation})'
                    if description:
                        detail += f': {description}'
                    param_details.append(detail)
                embed.add_field(
                    name='Usage',
                    value=f'`{' '.join(usage_parts)}`',
                    inline=False
                )
                if param_details:
                    embed.add_field(
                        name='Parameter Details',
                        value='\n'.join(param_details),
                        inline=False
                    )
            return await self.handler.send_message(ctx, embed=embed)
        all_commands = await self.get_available_commands(bot, ctx)
        if not all_commands:
            return await self.handler.send_message(ctx, '\U0001F6AB No commands available to you.')
        current_text_channel_id = ctx.channel.id
        current_room_name = getattr(ctx.channel, 'name', None)
        current_guild_id = ctx.guild.id if ctx.guild else None
        guild_aliases = self.bot.command_aliases.get(current_guild_id, {})
        guild_channel_aliases = guild_aliases.get('channel_aliases', {})
        guild_role_aliases = guild_aliases.get('role_aliases', {})
        guild_temp_aliases = guild_aliases.get('temp_room_aliases', {})
        current_guild_alias_commands = set()
        for type_map in guild_channel_aliases.values():
            current_guild_alias_commands.update(type_map.keys())
        for type_map in guild_role_aliases.values():
            current_guild_alias_commands.update(type_map.keys())
        for type_map in guild_temp_aliases.values():
            current_guild_alias_commands.update(type_map.keys())
        all_alias_commands = set()
        for guild_id, guild_data in self.bot.command_aliases.items():
            for type_map in guild_data.get('channel_aliases', {}).values():
                all_alias_commands.update(type_map.keys())
            for type_map in guild_data.get('role_aliases', {}).values():
                all_alias_commands.update(type_map.keys())
            for type_map in guild_data.get('temp_room_aliases', {}).values():
                all_alias_commands.update(type_map.keys())
        contextual_commands = []
        for command in all_commands:
            in_channel = any(
                command.name in type_map and type_map[command.name] == current_text_channel_id
                for type_map in guild_channel_aliases.values()
            )
            in_role = any(
                command.name in type_map and isinstance(type_map[command.name], dict) and type_map[command.name].get('channel_id') == current_text_channel_id
                for type_map in guild_role_aliases.values()
            )
            in_temp_room = any(
                command.name in type_map and isinstance(type_map[command.name], dict) and type_map[command.name].get('room_name') == current_room_name
                for type_map in guild_temp_aliases.values()
            )
            if in_channel or in_role or in_temp_room or command.name not in all_alias_commands:
                contextual_commands.append(command)
        if not contextual_commands:
            return await self.handler.send_message(ctx, '\U0001F6AB No commands available to you.')
        permission_groups = await self.group_commands_by_permission(bot, ctx, contextual_commands)
        pages = []
        user_highest = await self.get_user_highest_permission(bot, ctx)
        user_index = PERMISSION_ORDER.index(user_highest)
        permission_order = [
            ('Owner', '`Owner` inherits `developer`.'),
            ('Developer', '`Developer` inherits `coordinator`.'),
            ('Coordinator', '`Coordinator` inherits `moderator`.'),
            ('Moderator', 'Moderators can use these commands.'),
            ('Everyone', 'Commands available to everyone.')
        ]
        for i, (perm_level, description) in enumerate(permission_order):
            if i < user_index:
                continue
            commands_in_level = sorted(permission_groups.get(perm_level, []), key=lambda c: c.name)
            if not commands_in_level:
                continue
            embed = discord.Embed(
                title=f'{perm_level} Commands',
                description=description,
                color=self.get_permission_color(perm_level)
            )
            cog_map = {}
            for command in commands_in_level:
                if command.hidden:
                    continue
                cog_name = command.cog_name or 'Aliases'
                cog_map.setdefault(cog_name, []).append(command)
            for cog_name in sorted(cog_map):
                commands_in_cog = sorted(cog_map[cog_name], key=lambda c: c.name)
                command_list = '\n'.join(
                    f'**{config["discord_command_prefix"]}{cmd.name}** – {cmd.help or "No description"}'
                    for cmd in commands_in_cog
                )
                if len(command_list) > 1024:
                    chunks = self.split_command_list(commands_in_cog)
                    for j, chunk in enumerate(chunks):
                        field_name = f'{cog_name}' if j == 0 else f'{cog_name} (cont.)'
                        embed.add_field(name=field_name, value=chunk, inline=False)
                else:
                    embed.add_field(name=cog_name, value=command_list, inline=False)
            pages.append(embed)
        async with self.bot.db_pool.acquire() as conn:
            rows = await self.bot.db_pool.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
        team_role_ids = [r['role_id'] for r in rows]
        author_role_ids = [role.id for role in ctx.author.roles]
        is_team_member = any(rid in team_role_ids for rid in author_role_ids)
        if is_team_member:
            team_member_commands = [cmd for cmd in all_commands if getattr(cmd.callback, "_team_command", False) and not cmd.hidden]
            if team_member_commands:
                embed = discord.Embed(
                    title="Team Member Commands",
                    description="Commands available to users in team member roles:",
                    color=discord.Color.yellow()
                )
                cog_map = {}
                for cmd in team_member_commands:
                    cog_map.setdefault(cmd.cog_name or "No Cog", []).append(cmd)
                for cog_name, cmds in cog_map.items():
                    cmd_list = '\n'.join(f'**{config["discord_command_prefix"]}{c.name}** – {c.help or "No description"}' for c in cmds)
                    embed.add_field(name=cog_name, value=cmd_list, inline=False)
                pages.insert(-1, embed)
        if not pages:
            return await self.handler.send_message(ctx, '\U0001F6AB No commands available to you.')
        paginator = Paginator(bot, ctx, pages)
        return await paginator.start()

async def setup(bot: DiscordBot):
    cog = Help(bot)
    await bot.add_cog(cog)
