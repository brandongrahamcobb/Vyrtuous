"""category.py The purpose of this program is to provide the Category properties class.

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

from pathlib import Path

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.utils.dir_to_classes import dir_to_classes
from vyrtuous.utils.logger import logger


class CategoryObject:

    EXTRA_CATEGORIES = ["all"]

    def __init__(self, category=None):
        self.__category = category

    @property
    def category(self):
        return self.__category

    @category.setter
    def category(self, new_cat):
        dir_paths = []
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/infractions")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/mgmt")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/roles")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/rooms")
        classes = dir_to_classes(dir_paths=dir_paths)
        categories = [obj.category for obj in classes if obj.category is not None]
        for extra in self.EXTRA_CATEGORIES:
            categories.append(extra)
        if new_cat not in categories:
            logger.warning(f"Invalid category type ({str(new_cat)}).")
        self.__category = new_cat


class Converter(commands.Converter):

    def __init__(self, category_cls=CategoryObject):
        self.category_cls = category_cls

    async def convert(self, ctx: commands.Context, arg):
        return self.category_cls(arg).category


class Transformer(app_commands.Transformer):

    def __init__(self, category_cls=CategoryObject):
        self.category_cls = category_cls

    async def transform(self, interaction: discord.Interaction, arg):
        return self.category_cls(arg).category


class Category(Converter):
    def __init__(self):
        super().__init__(CategoryObject)


class AppCategory(Transformer):
    def __init__(self):
        super().__init__(CategoryObject)
