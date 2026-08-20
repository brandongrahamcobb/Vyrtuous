"""!/bin/python3
help_text_command.py A discord.py cog containing a custom help text command for the Vyrtuous bot.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import inspect
from collections import defaultdict
from typing import Callable

import discord
from discord.ext import commands

from vyrtuous.aliases import alias_service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.alias import Alias
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.metadata import metadata
from vyrtuous.permissions import permission_service
from vyrtuous.utils.messaging.tick import Tick


def skip_text_command_help_discovery() -> Callable:
    async def predicate(ctx) -> bool:
        return True

    return commands.check(predicate)


class HelpTextCommand(commands.Cog):

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def get_channel_alias_help(
        self, channel_snowflake: int, guild_snowflake: int
    ) -> list[str]:
        lines = []
        database_factory: DatabaseFactory = DatabaseFactory(Alias)
        aliases = await database_factory.select(
            channel_snowflake=channel_snowflake,
            guild_snowflake=int(guild_snowflake),
            singular=False,
        )
        if aliases:
            grouped = defaultdict(list)
            for alias in aliases:
                grouped[alias.category].append(alias.alias_name)
            for category, alias_names in grouped.items():
                help_lines = alias_service.CATEGORY_TO_HELP.get(category, None)
                if not help_lines:
                    continue
                for alias_name in alias_names:
                    lines.append(f"**{alias_name}**")
                    for line in help_lines:
                        lines.append(f"• {line}")
        return lines

    async def get_available_commands(
        self,
    ) -> tuple[list[commands.Command], list[commands.Command]]:
        available = []
        skipped = []
        for command in self.__bot.commands:
            try:
                has_skip = False
                if command.checks:
                    has_skip = any(
                        getattr(
                            getattr(check, "__wrapped__", check),
                            "_skip_text_command_help_discovery",
                            False,
                        )
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

    def get_permission_color(self) -> discord.Color:
        return discord.Color.random()

    def build_alias_pages(
        self,
        permission_state: PermissionState,
        aliases: list[Alias],
        group_order: list[str],
        default_group_alias: str,
    ) -> dict[str, list[Alias]]:
        pages: dict[str, list[Alias]] = {alias: [] for alias in group_order}
        for a in aliases:
            node = alias_service.CATEGORY_TO_PERMISSION.get(a.category)
            if not node:
                pages[default_group_alias].append(a)
                continue
            lowest_match = None
            for alias_key in group_order:
                if permission_service.group_defines_permission(
                    permission_state, alias_key, node
                ):
                    lowest_match = alias_key
            if lowest_match is not None:
                pages[lowest_match].append(a)
        return pages

    def compute_group_page_order(
        self,
        permission_state: PermissionState,
    ) -> list[str]:
        groups = list(permission_state.groups.values())
        groups.sort(key=lambda g: (-len(g.ancestors), g.name))
        return [g.alias for g in groups]

    async def compute_visible_page_order(
        self,
        permission_state: PermissionState,
        member_snowflake: int,
        guild_snowflake: int | None,
        group_order: list[str],
    ) -> list[str]:
        assigned_groups = await permission_service.resolve_all_assigned_groups(
            permission_state=permission_state,
            member_snowflake=member_snowflake,
        )
        visible_groups = [
            group
            for group, assigned_guild_snowflake, _ in assigned_groups
            if assigned_guild_snowflake is None
            or assigned_guild_snowflake == guild_snowflake
        ]
        highest_group = max(
            visible_groups,
            key=lambda group: len(group.ancestors),
            default=None,
        )
        if highest_group is None:
            return []
        member_chain = {highest_group.alias, *highest_group.ancestors}
        return [alias for alias in group_order if alias in member_chain]

    def build_command_pages(
        self,
        permission_state: PermissionState,
        commands_list: list,
        group_order: list[str],
        default_group_alias: str,
    ) -> dict[str, list]:
        pages: dict[str, list] = {alias: [] for alias in group_order}
        for command in commands_list:
            node = getattr(command.callback, "permission", None)
            if not node:
                pages[default_group_alias].append(command)
                continue
            lowest_match = None
            for alias in group_order:
                if permission_service.group_defines_permission(
                    permission_state, alias, node
                ):
                    lowest_match = alias
            if lowest_match is not None:
                pages[lowest_match].append(command)
        return pages

    async def resolve_command_or_alias(
        self, source, name: str
    ) -> tuple[str | None, commands.Command | Alias | None]:
        database_factory: DatabaseFactory = DatabaseFactory(Alias)
        cmd = self.__bot.get_command(name.lower())
        if cmd:
            return ("command", cmd)
        alias = await database_factory.select(
            alias_name=name.lower(), guild_snowflake=source.guild.id, singular=True
        )
        if alias and alias.guild_snowflake == source.guild.id:
            return ("alias", alias)
        return (None, None)

    def split_command_list(self, commands_list, max_length=1024) -> list[str]:
        current_chunk: list[str] = []
        chunks: list[str] = []
        current_length = 0
        for cmd in commands_list:
            cmd_line = f"**{self.__bot.config['discord_command_prefix']}{cmd.name}** – {cmd.help or 'No description'}"
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

    def unwrap_callback(self, func) -> Callable:
        while hasattr(func, "__wrapped__"):
            func = func.__wrapped__
        return func

    async def get_permission_filtered_aliases(
        self, source
    ) -> defaultdict[str, list[str]]:
        bot = DiscordBot.get_instance()
        database_factory: DatabaseFactory = DatabaseFactory(Alias)
        permission_state: PermissionState = bot.registry.get(PermissionState)
        aliases = await database_factory.select(
            channel_snowflake=source.channel.id,
            guild_snowflake=source.guild.id,
            singular=False,
        )
        grouped = defaultdict(list)
        perm_alias_map = defaultdict(list)
        if aliases:
            for alias in aliases:
                grouped[alias.category].append(alias)
            for category, alias_list in grouped.items():
                group = permission_service.get_default_group(
                    permission_state=permission_state
                )
                perm_level = alias_service.CATEGORY_TO_PERMISSION.get(
                    category, group.name
                )
                for a in alias_list:
                    help_lines = alias_service.CATEGORY_TO_HELP.get(category, [])
                    perm_alias_map[perm_level].append(
                        f"**{self.__bot.config['discord_command_prefix']}{a.alias_name}**\n"
                        + "\n".join(f"• {line}" for line in help_lines)
                    )
        return perm_alias_map

    @commands.command(name="help", help="List commands.")
    @metadata(permission="command.info.help")
    async def help_text_command(
        self, ctx, *, command_name: str | None = None
    ) -> discord.Message | None:
        tick = Tick(bot=self.__bot, ctx=ctx)
        bot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        guild = ctx.guild
        if guild is None:
            return await tick.end(warning="This command must be used in a server.")
        channel = ctx.channel
        if not isinstance(
            channel, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
        ):
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        await permission_service.has_permissions_at_all(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            requested=["command.info.help"],
        )
        pages, param_details, parameters = [], [], []
        if command_name and command_name != "all":
            kind, obj = await self.resolve_command_or_alias(ctx, command_name)
            if not kind:
                return await tick.end(
                    warning=f"Command or alias `{command_name}` not found."
                )
            if kind == "command":
                cmd = obj
                if kind == "command":
                    cmd = obj
                    if not isinstance(cmd, commands.Command):
                        return None
                    embed = discord.Embed(
                        title=f"{self.__bot.config['discord_command_prefix']}{cmd.name}",
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
                        usage_parts = [
                            f"{self.__bot.config['discord_command_prefix']}{cmd.name}"
                        ]
                        for name, param in parameters:
                            param_desc = None
                            if isinstance(param.default, commands.Parameter):
                                param_desc = param.default.description
                            is_optional = param.default != inspect.Parameter.empty
                            usage_parts.append(
                                f"[{name}]" if is_optional else f"<{name}>"
                            )
                            if param_desc:
                                param_details.append(
                                    f"**{name}**{' (Optional)' if not param.required else ''}: {param_desc}"
                                )
                            else:
                                param_details.append(
                                    f"**{name}**{' (Optional)' if not param.required else ''}"
                                )
                        embed.add_field(
                            name="Usage",
                            value=f"`{' '.join(usage_parts)}`",
                            inline=False,
                        )
                        if param_details:
                            embed.add_field(
                                name="Parameters",
                                value="\n".join(param_details),
                                inline=False,
                            )
                    return await tick.end(success=embed)
            elif kind == "alias":
                alias = obj
                if not isinstance(alias, Alias):
                    return None
                help_lines = alias_service.CATEGORY_TO_HELP.get(alias.category, None)
                if not help_lines:
                    return await tick.end(
                        warning=f"No help available for `{alias.alias_name}`."
                    )
                embed = discord.Embed(
                    title=f"{self.__bot.config['discord_command_prefix']}{alias.alias_name}",
                    description=f"Alias for **{alias.category}**",
                    color=discord.Color.green(),
                )
                embed.add_field(
                    name="Usage",
                    value="\n".join(f"• {line}" for line in help_lines),
                    inline=False,
                )
                return await tick.end(success=embed)
        all = False
        if command_name and command_name == "all":
            all = True
        all_commands, skipped_commands = await self.get_available_commands()
        group_order = self.compute_group_page_order(permission_state)
        visible_order = await self.compute_visible_page_order(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            guild_snowflake=guild.id,
            group_order=group_order,
        )
        if not visible_order:
            return await tick.end(warning="No commands available to you.")
        default_group = next(
            group for group in permission_state.groups.values() if group.default
        )
        pages_by_group = self.build_command_pages(
            permission_state=permission_state,
            commands_list=all_commands,
            group_order=visible_order,
            default_group_alias=default_group.alias,
        )
        skipped_pages_by_group = self.build_command_pages(
            permission_state=permission_state,
            commands_list=skipped_commands,
            group_order=visible_order,
            default_group_alias=default_group.alias,
        )
        database_factory: DatabaseFactory = DatabaseFactory(Alias)
        aliases = await database_factory.select(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            singular=False,
        )
        alias_pages_by_group = self.build_alias_pages(
            permission_state=permission_state,
            aliases=aliases or [],
            group_order=visible_order,
            default_group_alias=default_group.alias,
        )
        for alias_key in visible_order:
            group = permission_state.groups[alias_key]
            if group.ancestors:
                parent = permission_state.groups.get(group.ancestors[0])
                description = (
                    f"`{group.name}` inherits `{parent.name}`."
                    if parent
                    else "Commands available to everyone."
                )
            else:
                description = "Commands available to everyone."
            commands_in_level = sorted(
                pages_by_group.get(alias_key, []), key=lambda c: c.name
            )
            embed = discord.Embed(
                title=f"{group.name} Commands",
                description=description,
                color=self.get_permission_color(),
            )
            if commands_in_level:
                command_lines = [
                    f"**{self.__bot.config['discord_command_prefix']}{cmd.name}** – {cmd.help or 'No description'}"
                    for cmd in commands_in_level
                ]
                command_text = "\n".join(command_lines)
                if len(command_text) > 1024:
                    chunks = self.split_command_list(commands_in_level)
                    for i, chunk in enumerate(chunks):
                        field_name = (
                            f"{group.name} Commands"
                            if i == 0
                            else f"{group.name} Commands (cont.)"
                        )
                        embed.add_field(name=field_name, value=chunk, inline=False)
                else:
                    embed.add_field(name="", value=command_text, inline=False)
            aliases_in_level = alias_pages_by_group.get(alias_key, [])
            if aliases_in_level:
                alias_lines = []
                for a in aliases_in_level:
                    short_desc = alias_service.CATEGORY_TO_DESCRIPTION.get(
                        a.category, "No description"
                    )
                    alias_lines.append(
                        f"**{self.__bot.config['discord_command_prefix']}{a.alias_name}** – {short_desc}"
                    )
                embed.add_field(
                    name="Aliases", value="\n".join(alias_lines), inline=False
                )
            if all:
                skipped_in_level = sorted(
                    skipped_pages_by_group.get(alias_key, []), key=lambda c: c.name
                )
                if skipped_in_level:
                    skipped_lines = [
                        f"**{self.__bot.config['discord_command_prefix']}{cmd.name}** – {cmd.help or 'No description'}"
                        for cmd in skipped_in_level
                    ]
                    skipped_text = "\n".join(skipped_lines)
                    if len(skipped_text) > 1024:
                        chunks = self.split_command_list(skipped_in_level)
                        for i, chunk in enumerate(chunks):
                            embed.add_field(
                                name=("Additional" if i == 0 else "Additional (cont.)"),
                                value=chunk,
                                inline=False,
                            )
                    else:
                        embed.add_field(
                            name="Additional", value=skipped_text, inline=False
                        )
            pages.append(embed)
        if not pages:
            return await tick.end(warning="No commands available to you.")
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HelpTextCommand(bot=bot))
