"""!/bin/python3
module.py The purpose of this program is to provide the Module properties class.

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
from rapidfuzz import fuzz

from vyrtuous.inc.helpers import DISCORD_COGS


class ModuleObject:

    def __init__(self, module: str):
        self.__module = module

    @property
    def module(self) -> str:
        return self.__module

    @module.setter
    def module(self, new_module: str) -> None:
        self.__module = find_best_cog(new_module, DISCORD_COGS)


def find_best_cog(query: str, cogs: list[str]) -> str:
    return max(
        cogs,
        key=lambda cog: fuzz.token_set_ratio(
            query, cog.rsplit(".", 1)[-1].replace("_", " ")
        ),
    )


class Converter(commands.Converter):

    def __init__(self, module_cls=ModuleObject):
        self.module_cls = module_cls

    async def convert(self, ctx: commands.Context, argument) -> ModuleObject:
        return self.module_cls(argument)


class Transformer(app_commands.Transformer):

    def __init__(self, module_cls=ModuleObject):
        self.module_cls = module_cls

    async def transform(self, interaction: discord.Interaction, arg) -> ModuleObject:
        return self.module_cls(arg)


class Module(Converter):
    def __init__(self):
        super().__init__(ModuleObject)


class AppModule(Transformer):
    def __init__(self):
        super().__init__(ModuleObject)
