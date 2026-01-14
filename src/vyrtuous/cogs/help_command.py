"""help_command.py A discord.py cog containing a custom help command for the Vyrtuous bot.

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
"""

from collections import defaultdict
from discord import app_commands
from discord.ext import commands
from typing import Optional
from vyrtuous.service.check_service import role_check_without_specifics
from vyrtuous.service.message_service import MessageService
from vyrtuous.utils.setup_logging import logger
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.permission import PERMISSION_TYPES
from vyrtuous.database.actions.alias import Alias
from vyrtuous.service.state_service import StateService

import discord
import inspect

class HelpCommand(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.aliases_cog = bot.get_cog("Aliases")
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.message_service = MessageService(self.bot, self.bot.db_pool)
        self.permission_page_title_pairs = [
            ("System Owner", "`System Owner` inherits `Developer`."),
            ("Developer", "`Developer` inherits `Guild Owner`."),
            ("Guild Owner", "`Guild Owner` inherits `Administrator`."),
            ("Administrator", "`Administrator` inherits `Coordinator`."),
            ("Coordinator", "`Coordinator` inherits `Moderator`."),
            ("Moderator", "`Moderator` inherits `Everyone`."),
            ("Everyone", "Commands available to everyone."),
        ]

    async def get_channel_alias_help(
        self, channel_snowflake: Optional[int], guild_snowflake: Optional[int]
    ) -> list[str]:
        lines = []
        aliases = await Alias.select(guild_snowflake=guild_snowflake)
        if aliases:
            grouped = defaultdict(list)
            for alias in aliases:
                grouped[alias.alias_type].append(alias.alias_name)
            for alias_type, alias_names in grouped.items():
                help_lines = self.aliases_cog.alias_help.get(alias_type)
                if not help_lines:
                    continue
                for alias_name in alias_names:
                    lines.append(f"**{alias_name}**")
                    for line in help_lines:
                        lines.append(f"• {line}")
            return lines

    async def get_available_commands(
        self, bot, ctx_or_interaction
    ) -> list[commands.Command]:
        available = []
        user_highest = await role_check_without_specifics(ctx_or_interaction)
        for command in bot.commands:
            try:
                perm_level = await self.get_command_permission_level(bot, command)
                if PERMISSION_TYPES.index(user_highest) >= PERMISSION_TYPES.index(
                    perm_level
                ):
                    available.append(command)
            except Exception as e:
                logger.warning(
                    f"Exception while evaluating command {command}: {str(e).capitalize()}"
                )
        return available

    async def get_command_permission_level(self, bot, command):
        if not hasattr(command, "checks") or not command.checks:
            return "Everyone"
        for check in command.checks:
            if hasattr(check, "__wrapped__"):
                func = check.__wrapped__
            else:
                func = check
            if hasattr(func, "_permission_level"):
                return func._permission_level
        return "Everyone"

    def get_permission_color(self, perm_level):
        colors = {
            "System Owner": discord.Color.dark_red(),
            "Developer": discord.Color.red(),
            "Guild Owner": discord.Color.purple(),
            "Administrator": discord.Color.blue(),
            "Coordinator": discord.Color.orange(),
            "Moderator": discord.Color.green(),
            "Everyone": discord.Color.gold(),
        }
        return colors.get(perm_level, discord.Color.greyple())

    async def group_commands_by_permission(
        self, bot, ctx_or_interaction, commands_list
    ):
        permission_groups = {level: [] for level in PERMISSION_TYPES}
        for command in commands_list:
            perm_level = await self.get_command_permission_level(bot, command)
            if perm_level in permission_groups:
                permission_groups[perm_level].append(command)
            else:
                permission_groups["Everyone"].append(command)
        return permission_groups

    async def resolve_command_or_alias(self, ctx_or_interaction, name: str):
        cmd = self.bot.get_command(name.lower())
        if cmd:
            return ("command", cmd)
        alias = await Alias.select(
            alias_name=name.lower(), guild_snowflake=ctx_or_interaction.guild.id
        )
        if alias and alias.guild_snowflake == ctx_or_interaction.guild.id:
            return ("alias", alias)
        return (None, None)

    def split_command_list(self, commands_list, max_length=1024):
        current_chunk, chunks = [], []
        current_length = 0
        for cmd in commands_list:
            cmd_line = f"**{self.config['discord_command_prefix']}{cmd.name}** – {cmd.help or 'No description'}"
            cmd_length = len(cmd_line)
            if current_length + cmd_length > max_length and current_chunk:
                chunks.append("\n".join(current_chunk))
                current_chunk = [cmd_line.rstrip()]
                current_length = cmd_length
            else:
                current_chunk.append(cmd_line.rstrip())
                current_length += cmd_length
        if current_chunk:
            chunks.append("\n".join(current_chunk))
        return chunks

    def unwrap_callback(self, func):
        while hasattr(func, "__wrapped__"):
            func = func.__wrapped__
        return func

    async def get_permission_filtered_aliases(self, ctx_or_interaction):
        aliases = await Alias.select(
            channel_snowflake=ctx_or_interaction.channel.id,
            guild_snowflake=ctx_or_interaction.guild.id,
        )
        if aliases:
            grouped = defaultdict(list)
            for alias in aliases:
                grouped[alias.alias_type].append(alias)
            perm_alias_map = defaultdict(list)
            for alias_type, alias_list in grouped.items():
                perm_level = self.aliases_cog.alias_type_to_permission_level.get(
                    alias_type, "Everyone"
                )
                for a in alias_list:
                    help_lines = self.aliases_cog.alias_help.get(alias_type, [])
                    perm_alias_map[perm_level].append(
                        f"**{a.alias_name}**\n"
                        + "\n".join(f"• {line}" for line in help_lines)
                    )
            return perm_alias_map

    @app_commands.command(
        name="help", description="Show command information or your available commands."
    )
    @app_commands.describe(command_name="The command to view details for.")
    async def help_app_command(
        self, interaction: discord.Interaction, command_name: Optional[str] = None
    ):
        state = StateService(interaction)
        bot = interaction.client
        pages, param_details, parameters = [], [], []
        if command_name:
            kind, obj = await self.resolve_command_or_alias(interaction, command_name)
            if not kind:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f Command or alias `{command_name}` not found."
                    )
                except Exception as e:
                    try:
                        return await state.end(
                            warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                        )
                    except Exception as e:
                        return await state.end(error=f"\u274c {str(e).capitalize()}")
            if kind == "command":
                cmd = obj
                embed = discord.Embed(
                    title=f"{self.config['discord_command_prefix']}{cmd.name}",
                    description=cmd.help or "No description provided.",
                    color=discord.Color.blue(),
                )
                callback = self.unwrap_callback(cmd.callback)
                sig = inspect.signature(callback)
                for name, param in sig.parameters.items():
                    if param.annotation == discord.Interaction:
                        continue
                    parameters.append((name, param))
                if parameters and parameters[0][0] == "self":
                    parameters.pop(0)
                if parameters and parameters[0][0] == "ctx":
                    parameters.pop(0)
                if parameters:
                    usage_parts = [f"{self.config['discord_command_prefix']}{cmd.name}"]
                    for name, param in parameters:
                        param_desc = None
                        if isinstance(param.default, commands.Parameter):
                            param_desc = param.default.description
                        annotation = (
                            param.annotation.__name__
                            if hasattr(param.annotation, "__name__")
                            else str(param.annotation)
                        )
                        if annotation == "_empty":
                            annotation = "any"
                        is_optional = param.default != inspect.Parameter.empty
                        usage_parts.append(f"[{name}]" if is_optional else f"<{name}>")
                        if param_desc:
                            param_details.append(
                                f"**{name}** ({annotation}): {param_desc}"
                            )
                        else:
                            param_details.append(f"**{name}** ({annotation})")
                    embed.add_field(
                        name="Usage", value=f"`{' '.join(usage_parts)}`", inline=False
                    )
                    if param_details:
                        embed.add_field(
                            name="Parameters",
                            value="\n".join(param_details),
                            inline=False,
                        )
                        try:
                            return await state.end(success=embed)
                        except Exception as e:
                            try:
                                return await state.end(
                                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                                )
                            except Exception as e:
                                return await state.end(
                                    error=f"\u274c {str(e).capitalize()}"
                                )
            if kind == "alias":
                alias = obj
                help_lines = self.aliases_cog.alias_help.get(alias.alias_type)
                if not help_lines:
                    try:
                        return await state.end(
                            warning=f"\U000026a0\U0000fe0f No help available for `{alias.alias_name}`."
                        )
                    except Exception as e:
                        try:
                            return await state.end(
                                warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                            )
                        except Exception as e:
                            return await state.end(
                                error=f"\u274c {str(e).capitalize()}"
                            )
                embed = discord.Embed(
                    title=f"{self.config['discord_command_prefix']}{alias.alias_name}",
                    description=f"Alias for **{alias.alias_type}**",
                    color=discord.Color.green(),
                )
                embed.add_field(
                    name="Usage",
                    value="\n".join(f"• {line}" for line in help_lines),
                    inline=False,
                )
                try:
                    return await state.end(success=embed)
                except Exception as e:
                    try:
                        return await state.end(
                            warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                        )
                    except Exception as e:
                        return await state.end(error=f"\u274c {str(e).capitalize()}")
        all_commands = await self.get_available_commands(bot, interaction)
        permission_groups = await self.group_commands_by_permission(
            bot, interaction, all_commands
        )
        aliases = await Alias.select(
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
        )
        perm_alias_map = defaultdict(list)
        if aliases:
            for alias in aliases:
                short_desc = self.aliases_cog.alias_type_to_description.get(
                    alias.alias_type, "No description"
                )
                perm_level_for_alias = (
                    self.aliases_cog.alias_type_to_permission_level.get(
                        alias.alias_type, "Everyone"
                    )
                )
                perm_alias_map[perm_level_for_alias].append(
                    f"**{alias.alias_name}** – {short_desc}"
                )
        user_highest = await role_check_without_specifics(interaction)
        user_index = PERMISSION_TYPES.index(user_highest)
        for perm_level, description in self.permission_page_title_pairs:
            if PERMISSION_TYPES.index(perm_level) > user_index:
                continue
            commands_in_level = sorted(
                permission_groups.get(perm_level, []), key=lambda c: c.name
            )
            embed = discord.Embed(
                title=f"{perm_level} Commands",
                description=description,
                color=self.get_permission_color(perm_level),
            )
            if commands_in_level:
                command_lines = [
                    f"**{self.config['discord_command_prefix']}{cmd.name}** – {cmd.help or 'No description'}"
                    for cmd in commands_in_level
                ]
                command_text = "\n".join(command_lines)
                if len(command_text) > 1024:
                    chunks = self.split_command_list(commands_in_level)
                    for i, chunk in enumerate(chunks):
                        field_name = (
                            f"{perm_level} Commands"
                            if i == 0
                            else f"{perm_level} Commands (cont.)"
                        )
                        embed.add_field(name=field_name, value=chunk, inline=False)
                else:
                    embed.add_field(name="", value=command_text, inline=False)
            if perm_level in perm_alias_map:
                aliases_text = "\n".join(perm_alias_map[perm_level])
                embed.add_field(name="Aliases", value=aliases_text, inline=False)
            pages.append(embed)
        if not pages:
            try:
                return await state.end(
                    warning="\U000026a0\U0000fe0f No commands available to you."
                )
            except Exception as e:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(success=pages)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

    @commands.command(name="help")
    async def help_text_command(self, ctx, *, command_name: str = None):
        state = StateService(ctx)
        bot = ctx.bot
        pages, param_details, parameters = [], [], []
        if command_name:
            kind, obj = await self.resolve_command_or_alias(ctx, command_name)
            if not kind:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f Command or alias `{command_name}` not found."
                    )
                except Exception as e:
                    try:
                        return await state.end(
                            warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                        )
                    except Exception as e:
                        return await state.end(error=f"\u274c {str(e).capitalize()}")
            if kind == "command":
                cmd = obj
                embed = discord.Embed(
                    title=f"{self.config['discord_command_prefix']}{cmd.name}",
                    description=cmd.help or "No description provided.",
                    color=discord.Color.blue(),
                )
                callback = self.unwrap_callback(cmd.callback)
                sig = inspect.signature(callback)
                for name, param in sig.parameters.items():
                    if param.annotation == commands.Context:
                        continue
                    parameters.append((name, param))
                if parameters and parameters[0][0] == "self":
                    parameters.pop(0)
                if parameters and parameters[0][0] == "ctx":
                    parameters.pop(0)
                if parameters:
                    usage_parts = [f"{self.config['discord_command_prefix']}{cmd.name}"]
                    for name, param in parameters:
                        param_desc = None
                        if isinstance(param.default, commands.Parameter):
                            param_desc = param.default.description
                        annotation = (
                            param.annotation.__name__
                            if hasattr(param.annotation, "__name__")
                            else str(param.annotation)
                        )
                        if annotation == "_empty":
                            annotation = "any"
                        is_optional = param.default != inspect.Parameter.empty
                        usage_parts.append(f"[{name}]" if is_optional else f"<{name}>")
                        if param_desc:
                            param_details.append(
                                f"**{name}** ({annotation}): {param_desc}"
                            )
                        else:
                            param_details.append(f"**{name}** ({annotation})")
                    embed.add_field(
                        name="Usage", value=f"`{' '.join(usage_parts)}`", inline=False
                    )
                    if param_details:
                        embed.add_field(
                            name="Parameters",
                            value="\n".join(param_details),
                            inline=False,
                        )
                        try:
                            return await state.end(success=embed)
                        except Exception as e:
                            try:
                                return await state.end(
                                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                                )
                            except Exception as e:
                                return await state.end(
                                    error=f"\u274c {str(e).capitalize()}"
                                )
            if kind == "alias":
                alias = obj
                help_lines = self.aliases_cog.alias_help.get(alias.alias_type)
                if not help_lines:
                    try:
                        return await state.end(
                            warning=f"\U000026a0\U0000fe0f No help available for `{alias.alias_name}`."
                        )
                    except Exception as e:
                        try:
                            return await state.end(
                                warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                            )
                        except Exception as e:
                            return await state.end(
                                error=f"\u274c {str(e).capitalize()}"
                            )
                embed = discord.Embed(
                    title=f"{self.config['discord_command_prefix']}{alias.alias_name}",
                    description=f"Alias for **{alias.alias_type}**",
                    color=discord.Color.green(),
                )
                embed.add_field(
                    name="Usage",
                    value="\n".join(f"• {line}" for line in help_lines),
                    inline=False,
                )
                try:
                    return await state.end(success=embed)
                except Exception as e:
                    try:
                        return await state.end(
                            warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                        )
                    except Exception as e:
                        return await state.end(error=f"\u274c {str(e).capitalize()}")
        all_commands = await self.get_available_commands(bot, ctx)
        permission_groups = await self.group_commands_by_permission(
            bot, ctx, all_commands
        )
        aliases = await Alias.select(
            channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id
        )
        perm_alias_map = defaultdict(list)
        if aliases:
            for alias in aliases:
                short_desc = self.aliases_cog.alias_type_to_description.get(
                    alias.alias_type, "No description"
                )
                perm_level_for_alias = (
                    self.aliases_cog.alias_type_to_permission_level.get(
                        alias.alias_type, "Everyone"
                    )
                )
                perm_alias_map[perm_level_for_alias].append(
                    f"**{alias.alias_name}** – {short_desc}"
                )
        user_highest = await role_check_without_specifics(ctx)
        user_index = PERMISSION_TYPES.index(user_highest)
        for perm_level, description in self.permission_page_title_pairs:
            if PERMISSION_TYPES.index(perm_level) > user_index:
                continue
            commands_in_level = sorted(
                permission_groups.get(perm_level, []), key=lambda c: c.name
            )
            embed = discord.Embed(
                title=f"{perm_level} Commands",
                description=description,
                color=self.get_permission_color(perm_level),
            )
            if commands_in_level:
                command_lines = [
                    f"**{self.config['discord_command_prefix']}{cmd.name}** – {cmd.help or 'No description'}"
                    for cmd in commands_in_level
                ]
                command_text = "\n".join(command_lines)
                if len(command_text) > 1024:
                    chunks = self.split_command_list(commands_in_level)
                    for i, chunk in enumerate(chunks):
                        field_name = (
                            f"{perm_level} Commands"
                            if i == 0
                            else f"{perm_level} Commands (cont.)"
                        )
                        embed.add_field(name=field_name, value=chunk, inline=False)
                else:
                    embed.add_field(name="", value=command_text, inline=False)
            if perm_level in perm_alias_map:
                aliases_text = "\n".join(perm_alias_map[perm_level])
                embed.add_field(name="Aliases", value=aliases_text, inline=False)
            pages.append(embed)
        if not pages:
            try:
                return await state.end(
                    warning="\U000026a0\U0000fe0f No commands available to you."
                )
            except Exception as e:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(success=pages)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

async def setup(bot: DiscordBot):
    cog = HelpCommand(bot)
    await bot.add_cog(cog)
