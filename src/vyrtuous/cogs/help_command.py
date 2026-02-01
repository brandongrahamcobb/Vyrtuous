"""help_command.py A discord.py cog containing a custom help command for the Vyrtuous bot.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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

import inspect
from collections import defaultdict

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.inc.helpers import PERMISSION_TYPES
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.roles.moderator_service import moderator_predicator
from vyrtuous.service.state_service import StateService
from vyrtuous.utils.highest_role import resolve_highest_role
from vyrtuous.utils.logger import logger


def skip_help_discovery():
    async def predicate(ctx):
        return True

    predicate._skip_help_discovery = True
    return commands.check(predicate)


class HelpCommand(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.message_service = MessageService(self.bot)
        self.permission_page_title_pairs = [
            ("Sysadmin", "`Sysadmin` inherits `Developer`."),
            ("Developer", "`Developer` inherits `Guild Owner`."),
            ("Guild Owner", "`Guild Owner` inherits `Administrator`."),
            ("Administrator", "`Administrator` inherits `Coordinator`."),
            ("Coordinator", "`Coordinator` inherits `Moderator`."),
            ("Moderator", "`Moderator` inherits `Everyone`."),
            ("Everyone", "Commands available to everyone."),
        ]

    async def get_channel_alias_help(
        self, channel_snowflake: int, guild_snowflake: int
    ) -> list[str]:
        lines = []
        aliases = await Alias.select(guild_snowflake=int(guild_snowflake))
        if aliases:
            grouped = defaultdict(list)
            for alias in aliases:
                grouped[alias.category].append(alias.alias_name)
            for category, alias_names in grouped.items():
                help_lines = Alias.alias_help.get(category, None)
                if not help_lines:
                    continue
                for alias_name in alias_names:
                    lines.append(f"**{alias_name}**")
                    for line in help_lines:
                        lines.append(f"• {line}")
            return lines

    async def get_available_commands(
        self, all, bot, user_highest
    ) -> list[commands.Command]:
        available = []
        for command in bot.commands:
            try:
                perm_level = await self.get_command_permission_level(bot, command)
                if PERMISSION_TYPES.index(user_highest) >= PERMISSION_TYPES.index(
                    perm_level
                ):
                    fail = False
                    for check in command.checks:
                        if hasattr(check, "_skip_help_discovery") and not all:
                            fail = True
                    if not fail:
                        available.append(command)
            except Exception as e:
                logger.warning(
                    f"Exception while evaluating command {command}: {str(e).capitalize()}"
                )
        return available

    async def get_command_permission_level(self, bot, command):
        if not hasattr(command, "checks") or not command.checks:
            return "Everyone"
        for verify in command.checks:
            if hasattr(verify, "__wrapped__"):
                func = verify.__wrapped__
            else:
                func = verify
            if hasattr(func, "_permission_level"):
                return func._permission_level
        return "Everyone"

    def get_permission_color(self, perm_level):
        colors = {
            "Sysadmin": discord.Color.dark_red(),
            "Developer": discord.Color.red(),
            "Guild Owner": discord.Color.purple(),
            "Administrator": discord.Color.blue(),
            "Coordinator": discord.Color.orange(),
            "Moderator": discord.Color.green(),
            "Everyone": discord.Color.gold(),
        }
        return colors.get(perm_level, discord.Color.greyple())

    async def group_commands_by_permission(self, bot, source, commands_list):
        permission_groups = {level: [] for level in PERMISSION_TYPES}
        for command in commands_list:
            perm_level = await self.get_command_permission_level(bot, command)
            if perm_level in permission_groups:
                permission_groups[perm_level].append(command)
            else:
                permission_groups["Everyone"].append(command)
        return permission_groups

    async def resolve_command_or_alias(self, source, name: str):
        cmd = self.bot.get_command(name.lower())
        if cmd:
            return ("command", cmd)
        alias = await Alias.select(
            alias_name=name.lower(), guild_snowflake=source.guild.id, singular=True
        )
        if alias and alias.guild_snowflake == source.guild.id:
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

    async def get_permission_filtered_aliases(self, source):
        aliases = await Alias.select(
            channel_snowflake=source.channel.id,
            guild_snowflake=source.guild.id,
        )
        if aliases:
            grouped = defaultdict(list)
            for alias in aliases:
                grouped[alias.category].append(alias)
            perm_alias_map = defaultdict(list)
            for category, alias_list in grouped.items():
                perm_level = Alias.category_to_permission_level.get(
                    category, "Everyone"
                )
                for a in alias_list:
                    help_lines = Alias.alias_help.get(category, [])
                    perm_alias_map[perm_level].append(
                        f"**{self.config['discord_command_prefix']}{a.alias_name}**\n"
                        + "\n".join(f"• {line}" for line in help_lines)
                    )
            return perm_alias_map

    @app_commands.command(
        name="help", description="Show command information or your available commands."
    )
    @app_commands.describe(command_name="The command to view details for.")
    @moderator_predicator()
    async def help_app_command(
        self, interaction: discord.Interaction, command_name: str
    ):
        state = StateService(interaction=interaction)
        bot = interaction.client
        pages, param_details, parameters = [], [], []
        if command_name and command_name != "all":
            kind, obj = await self.resolve_command_or_alias(interaction, command_name)
            if not kind:
                return await state.end(
                    warning=f"Command or alias `{command_name}` not found."
                )
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
                        is_optional = param.default != inspect.Parameter.empty
                        usage_parts.append(f"[{name}]" if is_optional else f"<{name}>")
                        if param_desc:
                            param_details.append(f"**{name}**: {param_desc}")
                        else:
                            param_details.append(f"**{name}**")
                    embed.add_field(
                        name="Usage", value=f"`{' '.join(usage_parts)}`", inline=False
                    )
                    if param_details:
                        embed.add_field(
                            name="Parameters",
                            value="\n".join(param_details),
                            inline=False,
                        )
                        return await state.end(success=embed)
                alias = obj
                help_lines = Alias.alias_help.get(alias.category, None)
                if not help_lines:
                    return await state.end(
                        warning=f"No help available for `{alias.alias_name}`."
                    )
                embed = discord.Embed(
                    title=f"{self.config['discord_command_prefix']}{alias.alias_name}",
                    description=f"Alias for **{alias.category}**",
                    color=discord.Color.green(),
                )
                embed.add_field(
                    name="Usage",
                    value="\n".join(f"• {line}" for line in help_lines),
                    inline=False,
                )
                return await state.end(success=embed)
        all = False
        if command_name and command_name == "all":
            all = True
        user_highest = await resolve_highest_role(
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
        )
        all_commands = await self.get_available_commands(
            all=all, bot=bot, user_highest=user_highest
        )
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
                short_desc = Alias.category_to_description.get(
                    alias.category, "No description"
                )
                perm_level_for_alias = Alias.category_to_permission_level.get(
                    alias.category, "Everyone"
                )
                perm_alias_map[perm_level_for_alias].append(
                    f"**{self.config['discord_command_prefix']}{alias.alias_name}** – {short_desc}"
                )
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
            return await state.end(
                warning="\U000026a0\U0000fe0f No commands available to you."
            )
        return await state.end(success=pages)

    @commands.command(name="help")
    @moderator_predicator()
    async def help_text_command(self, ctx, *, command_name: str):
        state = StateService(ctx=ctx)
        bot = ctx.bot
        pages, param_details, parameters = [], [], []
        if command_name and command_name != "all":
            kind, obj = await self.resolve_command_or_alias(ctx, command_name)
            if not kind:
                return await state.end(
                    warning=f"Command or alias `{command_name}` not found."
                )
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
                        is_optional = param.default != inspect.Parameter.empty
                        usage_parts.append(f"[{name}]" if is_optional else f"<{name}>")
                        if param_desc:
                            param_details.append(f"**{name}**: {param_desc}")
                        else:
                            param_details.append(f"**{name}**")
                    embed.add_field(
                        name="Usage", value=f"`{' '.join(usage_parts)}`", inline=False
                    )
                    if param_details:
                        embed.add_field(
                            name="Parameters",
                            value="\n".join(param_details),
                            inline=False,
                        )
                        return await state.end(success=embed)
            if kind == "alias":
                alias = obj
                help_lines = Alias.alias_help.get(alias.category, None)
                if not help_lines:
                    return await state.end(
                        warning=f"No help available for `{alias.alias_name}`."
                    )
                embed = discord.Embed(
                    title=f"{self.config['discord_command_prefix']}{alias.alias_name}",
                    description=f"Alias for **{alias.category}**",
                    color=discord.Color.green(),
                )
                embed.add_field(
                    name="Usage",
                    value="\n".join(f"• {line}" for line in help_lines),
                    inline=False,
                )
                return await state.end(success=embed)
        all = False
        if command_name and command_name == "all":
            all = True
        user_highest = await resolve_highest_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
        )
        all_commands = await self.get_available_commands(
            all=all, bot=bot, user_highest=user_highest
        )
        permission_groups = await self.group_commands_by_permission(
            bot, ctx, all_commands
        )
        aliases = await Alias.select(
            channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id
        )
        perm_alias_map = defaultdict(list)
        if aliases:
            for alias in aliases:
                short_desc = Alias.category_to_description.get(
                    alias.category, "No description"
                )
                perm_level_for_alias = Alias.category_to_permission_level.get(
                    alias.category, "Everyone"
                )
                perm_alias_map[perm_level_for_alias].append(
                    f"**{self.config['discord_command_prefix']}{alias.alias_name}** – {short_desc}"
                )
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
            return await state.end(warning="No commands available to you.")
        return await state.end(success=pages)


async def setup(bot: DiscordBot):
    cog = HelpCommand(bot)
    await bot.add_cog(cog)
