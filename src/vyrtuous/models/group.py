"""!/bin/python3
target.py The purpose of this program is to provide the Target properties class.

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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.permissions import PermissionGroup
from vyrtuous.cache.registry import PermissionState
from vyrtuous.utils.errors.error import BadArgument


class GroupObject:

    def __init__(self, group: PermissionGroup):
        self.__group = group

    @property
    def group(self) -> PermissionGroup:
        return self.__group

    @group.setter
    def group(self, new_group: PermissionGroup) -> None:
        self.__group = new_group


class Converter(commands.Converter):

    def __init__(self, group_cls=GroupObject):
        self.group_cls = group_cls

    async def convert(self, ctx: commands.Context, argument: str) -> GroupObject | None:
        bot: DiscordBot = DiscordBot.get_instance()
        groups = bot.registry.get(PermissionState).groups
        for group in groups.values():
            if (
                argument.lower() == group.alias
                or argument == group.name
                or argument.lower() == group.name
            ):
                return self.group_cls(group)
        else:
            raise BadArgument("Group alias/name is invalid.")


class Transformer(app_commands.Transformer):

    def __init__(self, group_cls=GroupObject):
        self.group_cls = group_cls

    async def transform(
        self, interaction: discord.Interaction, argument: str
    ) -> GroupObject | None:
        bot: DiscordBot = DiscordBot.get_instance()
        groups = bot.registry.get(PermissionState).groups
        for group in groups.values():
            if (
                argument.lower() == group.alias
                or argument == group.name
                or argument.lower() == group.name
            ):
                return self.group_cls(group)
        else:
            raise app_commands.TransformerError(
                argument, discord.AppCommandOptionType.string, self
            )


class Group(Converter):
    def __init__(self):
        super().__init__(GroupObject)


class AppGroup(Transformer):
    def __init__(self):
        super().__init__(GroupObject)
