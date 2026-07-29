"""!/bin/python3
help_command.py A discord.py cog containing a custom help command for the Vyrtuous bot.

Copyright (C) 2026  https://{ithub.com/brandongrahamcobb/Vyrtuous.git

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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.models.metadata import metadata
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.permissions import permission_service


def skip_app_command_help_discovery():
    async def predicate(interaction):
        return True

    return app_commands.check(predicate)


class HelpAppCommand(commands.Cog):

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def get_available_commands(self, interaction):
        available = []
        skipped = []
        for command in self.__bot.tree.get_commands():
            if not isinstance(command, discord.app_commands.Command):
                continue
            try:
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

    def get_permission_color(self):
        return discord.Color.random()

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
        channel_snowflake: int | None,
        guild_snowflake: int | None,
        group_order: list[str],
    ) -> list[str]:
        effective_group = await permission_service.resolve_effective_group(
            permission_state=permission_state,
            member_snowflake=member_snowflake,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
        )
        if effective_group is None:
            return []
        member_chain = {effective_group.alias, *effective_group.ancestors}
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
            node = getattr(command, "permission", None)
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

    async def resolve_app_command(self, name: str):
        command = self.__bot.tree.get_command(name.lower(), guild=None)
        if command:
            return ("command", command)
        return (None, None)

    def split_command_list(self, commands_list, max_length=1024):
        current_chunk, chunks = [], []
        current_length = 0
        for cmd in commands_list:
            cmd_line = f"**/{cmd.name}** – {cmd.description or 'No description'}"
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

    @metadata(permission="command.info.help")
    @app_commands.command(name="help", description="List commands")
    @app_commands.describe(command_name="The name of the command.")
    async def help_app_command(
        self, interaction: discord.Interaction, *, command_name: str | None = None
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        guild = interaction.guild
        if guild is None:
            return await tick.end(warning="This command must be used in a server.")
        channel = interaction.channel
        if not isinstance(
            channel, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
        ):
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        await permission_service.any_group_has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            requested=["command.info.help"],
        )
        pages, param_details, parameters = [], [], []
        if command_name and command_name != "all":
            kind, obj = await self.resolve_app_command(command_name)
            if not kind:
                return await tick.end(warning=f"Command `{command_name}` not found.")
            if kind == "command":
                cmd = obj
                if isinstance(cmd, app_commands.Command):
                    embed = discord.Embed(
                        title=f"/{cmd.name}",
                        description=cmd.description or "No description provided.",
                        color=discord.Color.blue(),
                    )
                    if cmd.parameters:
                        usage_parts = [f"/{cmd.name}"]
                        param_details = []
                        for param in cmd.parameters:
                            usage_parts.append(
                                f"[{param.name}]"
                                if not param.required
                                else f"<{param.name}>"
                            )
                            if param.description and param.description != "…":
                                param_details.append(
                                    f"**{param.name}**: {param.description}"
                                )
                            else:
                                param_details.append(f"**{param.name}**")
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
        all = False
        if command_name and command_name == "all":
            all = True
        all_commands, skipped_commands = await self.get_available_commands(
            interaction=interaction
        )
        group_order = self.compute_group_page_order(permission_state)
        visible_order = await self.compute_visible_page_order(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel.id,
            guild_snowflake=guild.id,
            group_order=group_order,
        )
        if not visible_order:
            return await tick.end(
                warning="\U000026a0\U0000fe0f No commands available to you."
            )
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
        for alias in visible_order:
            group = permission_state.groups[alias]
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
                pages_by_group.get(alias, []), key=lambda c: c.name
            )
            embed = discord.Embed(
                title=f"{group.name} Commands",
                description=description,
                color=self.get_permission_color(),
            )
            if commands_in_level:
                command_lines = [
                    f"**/{cmd.name}** – {cmd.description or 'No description'}"
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
            if all:
                skipped_in_level = sorted(
                    skipped_pages_by_group.get(alias, []), key=lambda c: c.name
                )
                if skipped_in_level:
                    skipped_lines = [
                        f"**/{cmd.name}** – {cmd.description or 'No description'}"
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
            return await tick.end(
                warning="\U000026a0\U0000fe0f No commands available to you."
            )
        return await tick.end(success=pages)

    async def cog_app_command_error(self, interaction, error):
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))


async def setup(bot: DiscordBot):
    await bot.add_cog(HelpAppCommand(bot=bot))
