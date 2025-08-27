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
                logger.warning(f"❌ Exception while checking command '{command}': {e}")
        return available_commands

    async def get_command_permission_level(self, bot, ctx, command):
        if not hasattr(command, 'checks') or not command.checks:
            return 'Everyone'
        check_names = []
        for check in command.checks:
            func = check
            if hasattr(func, '__wrapped__'):  # unwrap decorators
                func = func.__wrapped__
            if hasattr(func, '__name__'):
                check_names.append(func.__name__)
        if any(name in ['is_owner', 'is_guild_owner', 'is_system_owner'] for name in check_names):
            return 'Owner'
        if 'is_owner_developer' in check_names:
            return 'Developer'
        if 'is_owner_developer_coordinator' in check_names:
            return 'Coordinator'
        if 'is_owner_developer_coordinator_moderator' in check_names:
            return 'Moderator'
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
        permission_checks = [
            ('Owner', is_owner),
            ('Developer', is_developer),
            ('Coordinator', is_coordinator),
            ('Moderator', is_moderator)
        ]
        for level, check in permission_checks:
            try:
                if await check(ctx):
                    return level
            except commands.CheckFailure:
                continue
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
            cmd_line = f'**{config['discord_command_prefix']}{cmd.name}** – {cmd.help or "No description"}\n'
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
                return await self.handler.send_message(ctx, f'❌ Command `{command_name}` not found.')
            if cmd.hidden:
                return await self.handler.send_message(ctx, f'❌ Command `{command_name}` is hidden.')
            if not await cmd.can_run(ctx):
                return await self.handler.send_message(ctx, f'❌ You do not have permission to run `{command_name}`.')
            embed = discord.Embed(
                title=f'{config['discord_command_prefix']}{cmd.name}',
                description=cmd.help or 'No description provided.',
                color=discord.Color.blue()
            )
            sig = inspect.signature(cmd.callback)
            parameters = list(sig.parameters.items())

            # Remove leading 'self' and 'ctx' if they exist
            if parameters and parameters[0][0] == 'self':
                parameters.pop(0)
            if parameters and parameters[0][0] == 'ctx':
                parameters.pop(0)
            if parameters:
                usage_parts = [f"{config['discord_command_prefix']}{cmd.name}"]
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
                        usage_parts.append(f"[{name}]")
                    else:
                        usage_parts.append(f"<{name}>")
                    detail = f"**{name}** ({annotation})"
                    if description:
                        detail += f": {description}"
                    param_details.append(detail)
                embed.add_field(
                    name="Usage",
                    value=f"`{' '.join(usage_parts)}`",
                    inline=False
                )
                if param_details:
                    embed.add_field(
                        name="Parameter Details",
                        value="\n".join(param_details),
                        inline=False
                    )
            return await self.handler.send_message(ctx, embed=embed)
        all_commands = await self.get_available_commands(bot, ctx)
        if not all_commands:
            return await self.handler.send_message(ctx, '❌ No commands available to you.')
        current_text_channel_id = ctx.channel.id
        guild_aliases = self.bot.command_aliases.get(ctx.guild.id, {})
        alias_to_channel_map = {}
        for alias_type_map in guild_aliases.values():
            for alias_name, channel_id in alias_type_map.items():
                alias_to_channel_map[alias_name] = channel_id
        contextual_commands = []
        for command in all_commands:
            alias_channel_id = alias_to_channel_map.get(command.name)

            # Command is aliased in this guild and must match the current channel
            if alias_channel_id is not None:
                if alias_channel_id == current_text_channel_id:
                    contextual_commands.append(command)
            # Command is not aliased in this guild at all — skip it
            elif command.name not in alias_to_channel_map:
                # Only include global commands (no alias binding anywhere)
                is_global = True
                for guild_id, alias_type_map in self.bot.command_aliases.items():
                    for channel_map in alias_type_map.values():
                        if command.name in channel_map:
                            is_global = False
                            break
                    if not is_global:
                        break
                if is_global:
                    contextual_commands.append(command)
        if not contextual_commands:
            await self.handler.send_message(ctx, '❌ No commands available to you.')
            return
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
                    f'**{config['discord_command_prefix']}{cmd.name}** – {cmd.help or "No description"}'
                    for cmd in commands_in_cog
                )
                if len(command_list) > 1024:
                    chunks = self.split_command_list(commands_in_cog)
                    for j, chunk in enumerate(chunks):
                        field_name = f'{cog_name}' if j == 0 else f'{cog_name} (cont.)'
                        embed.add_field(
                            name=field_name,
                            value=chunk,
                            inline=False
                        )
                else:
                    embed.add_field(
                        name=cog_name,
                        value=command_list,
                        inline=False
                    )
            pages.append(embed)
        if not pages:
            return await self.handler.send_message(ctx, '❌ No commands available to you.')
        paginator = Paginator(bot, ctx, pages)
        return await paginator.start()

async def setup(bot: DiscordBot):
    cog = Help(bot)
    await bot.add_cog(cog)
