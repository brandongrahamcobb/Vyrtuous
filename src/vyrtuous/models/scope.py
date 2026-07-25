"""!/bin/python3
scope.py The purpose of this program is to provide the Scope properties class.

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

SCOPES = ["all", "auto", "click", "command", "server"]


class ScopeObject:

    def __init__(self, scope: str):
        self.__scope = scope

    @property
    def scope(self) -> str:
        return self.__scope


class Converter(commands.Converter):

    def __init__(self, scope_cls=ScopeObject):
        self.scope_cls = scope_cls

    async def convert(self, ctx: commands.Context, argument: str) -> ScopeObject:
        if argument not in SCOPES:
            raise commands.CheckFailure(f"Invalid scope specified: ({argument}).")
        return self.scope_cls(argument)


class Transformer(app_commands.Transformer):

    def __init__(self, scope_cls=ScopeObject):
        self.scope_cls = scope_cls

    async def transform(
        self, interaction: discord.Interaction, arg: str
    ) -> ScopeObject:
        if arg not in SCOPES:
            raise app_commands.CheckFailure(f"Invalid scope specified: ({arg}).")
        return self.scope_cls(arg)


class Scope(Converter):
    def __init__(self):
        super().__init__(ScopeObject)


class AppScope(Transformer):
    def __init__(self):
        super().__init__(ScopeObject)
