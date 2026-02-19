"""!/bin/python3
help_command.py A discord.py cog containing a custom help command for the Vyrtuous bot.

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
from typing import Any, Coroutine, Union

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.administrator.administrator_service import AdministratorService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.moderator.moderator_service import ModeratorService, NotModerator
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.state_service import StateService


def skip_app_command_help_discovery():
    async def predicate(interaction):
        return True

    predicate._skip_app_command_help_discovery = True
    return app_commands.check(predicate)


class HelpCommand(commands.Cog):
    def __init__(self, *, bot: DiscordBot | None = None):
        self.__bot = bot
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__author_service = AuthorService()
        self.config = self.__bot.config
        self.permission_page_title_pairs = [
            ("Sysadmin", "`Sysadmin` inherits `Developer`."),
            ("Developer", "`Developer` inherits `Guild Owner`."),
            ("Guild Owner", "`Guild Owner` inherits `Administrator`."),
            ("Administrator", "`Administrator` inherits `Coordinator`."),
            ("Coordinator", "`Coordinator` inherits `Moderator`."),
            ("Moderator", "`Moderator` inherits `Everyone`."),
            ("Everyone", "Commands available to everyone."),
        ]
        self.__sysadmin_service = SysadminService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__bug_service = BugService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__moderator_service = ModeratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__developer_service = DeveloperService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            database_factory=self.__database_factory,
            emoji=self.__emoji,
        )
        self.__coordinator_service = CoordinatorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__administrator_service = AdministratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__guild_owner_service = GuildOwnerService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )

    async def cog_check(self, interaction) -> Coroutine[Any, Any, bool]:
        async def predicate(
            source: Union[commands.Context, discord.Interaction, discord.Message],
        ):
            for verify in (
                self.__sysadmin_service.is_sysadmin_wrapper,
                self.__developer_service.is_developer_wrapper,
                self.__guild_owner_service.is_guild_owner_wrapper,
                self.__administrator_service.is_administrator_wrapper,
                self.__coordinator_service.is_coordinator_at_all_wrapper,
                self.__moderator_service.is_moderator_at_all_wrapper,
            ):
                try:
                    if await verify(source):
                        return True
                except app_commands.CheckFailure:
                    continue
            raise NotModerator

        predicate._permission_level = "Moderator"
        return await predicate(interaction)

    async def get_available_commands(self, all, bot, user_highest):
        available = []
        skipped = []
        for command in bot.commands:
            try:
                perm_level = await self.get_command_permission_level(bot, command)
                if self.__moderator_service.PERMISSION_TYPES.index(
                    user_highest
                ) >= self.__moderator_service.PERMISSION_TYPES.index(perm_level):
                    has_skip = any(
                        hasattr(check, "_skip_app_command_help_discovery")
                        for check in command.checks
                    )
                    if has_skip:
                        skipped.append(command)
                    else:
                        available.append(command)
            except Exception as e:
                self.__bot.logger.warning(
                    f"Exception while evaluating command {command}: {str(e).capitalize()}"
                )
        return available, skipped

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
        permission_groups = {
            level: [] for level in self.__moderator_service.PERMISSION_TYPES
        }
        for command in commands_list:
            perm_level = await self.get_command_permission_level(bot, command)
            if perm_level in permission_groups:
                permission_groups[perm_level].append(command)
            else:
                permission_groups["Everyone"].append(command)
        return permission_groups

    async def resolve_app_command(self, source, name: str):
        command = self.bot.tree.get_command(name.lower(), guild=source.guild)
        if command:
            return ("command", command)
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

    @app_commands.command(name="help", description="List commands")
    @app_commands.describe(command_name="The name of the command.")
    async def help_app_command(
        self, interaction: discord.Interaction, *, command_name: str | None = None
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            interaction=interaction,
        )
        bot = interaction.client
        pages, param_details, parameters = [], [], []
        if command_name and command_name != "all":
            kind, obj = await self.resolve_app_command(interaction, command_name)
            if not kind:
                return await state.end(warning=f"Command `{command_name}` not found.")
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
                return await state.end(success=embed)
        all = False
        if command_name and command_name == "all":
            all = True
        user_highest = await self.__moderator_service.resolve_highest_role_at_all(
            member_snowflake=interaction.user.id,
        )
        all_commands, skipped_commands = await self.get_available_commands(
            all=all, bot=bot, user_highest=user_highest
        )
        permission_groups = await self.group_commands_by_permission(
            bot, interaction, all_commands
        )
        skipped_permission_groups = await self.group_commands_by_permission(
            bot, interaction, skipped_commands
        )
        user_index = self.__moderator_service.PERMISSION_TYPES.index(user_highest)
        for perm_level, description in self.permission_page_title_pairs:
            if self.__moderator_service.PERMISSION_TYPES.index(perm_level) > user_index:
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
            if all:
                skipped_in_level = sorted(
                    skipped_permission_groups.get(perm_level, []), key=lambda c: c.name
                )
                if skipped_in_level:
                    skipped_lines = [
                        f"**{self.config['discord_command_prefix']}{cmd.name}** – {cmd.help or 'No description'}"
                        for cmd in skipped_in_level
                    ]
                    skipped_text = "\n".join(skipped_lines)
                    if len(skipped_text) > 1024:
                        chunks = self.split_command_list(skipped_in_level)
                        for i, chunk in enumerate(chunks):
                            embed.add_field(
                                name="Additional" if i == 0 else "Additional (cont.)",
                                value=chunk,
                                inline=False,
                            )
                    else:
                        embed.add_field(
                            name="Additional", value=skipped_text, inline=False
                        )
            pages.append(embed)
        if not pages:
            return await state.end(
                warning="\U000026a0\U0000fe0f No commands available to you."
            )
        return await state.end(success=pages)


async def setup(bot: DiscordBot):
    cog = HelpCommand(bot)
    await bot.add_cog(cog)
