"""!/bin/python3
category.py The purpose of this program is to provide the Category properties class.

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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.db.active_member import ActiveMember
from vyrtuous.db.administrator import Administrator, AdministratorRole
from vyrtuous.db.alias import Alias
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.ban import Ban
from vyrtuous.db.bug import Bug
from vyrtuous.db.cap import Cap
from vyrtuous.db.coordinator import Coordinator
from vyrtuous.db.data import Data
from vyrtuous.db.developer import Developer
from vyrtuous.db.flag import Flag
from vyrtuous.db.moderator import Moderator
from vyrtuous.db.stream import Stream
from vyrtuous.db.text_mute import TextMute
from vyrtuous.db.vegan import Vegan
from vyrtuous.db.video_channel import VideoChannel
from vyrtuous.db.voice_mute import VoiceMute

CLASSES = [
    ActiveMember,
    Administrator,
    AdministratorRole,
    Alias,
    AutoMute,
    Ban,
    Bug,
    Cap,
    Coordinator,
    Data,
    Developer,
    Flag,
    Moderator,
    Stream,
    TextMute,
    Vegan,
    VideoChannel,
    VoiceMute,
]
EXTRA_CATEGORIES = ["all"]


class CategoryObject:

    def __init__(self, category: str):
        self.__category = category

    @property
    def category(self) -> str:
        return self.__category


class Converter(commands.Converter):

    def __init__(self, category_cls=CategoryObject):
        self.category_cls = category_cls
        self.__identifiers = [obj.identifier for obj in CLASSES]

    async def convert(self, ctx: commands.Context, argument: str) -> CategoryObject:
        categories = self.__identifiers
        for extra in EXTRA_CATEGORIES:
            categories.append(extra)
        if argument not in categories:
            raise commands.CheckFailure(f"Invalid category specified: ({argument}).")
        return self.category_cls(argument)


class Transformer(app_commands.Transformer):

    def __init__(self, category_cls=CategoryObject):
        self.category_cls = category_cls
        self.__identifiers = [obj.identifier for obj in CLASSES]

    async def transform(
        self, interaction: discord.Interaction, arg: str
    ) -> CategoryObject:
        categories = self.__identifiers
        for extra in EXTRA_CATEGORIES:
            categories.append(extra)
        if arg not in categories:
            raise app_commands.CheckFailure(f"Invalid category specified: ({arg}).")
        return self.category_cls(arg)


class Category(Converter):
    def __init__(self):
        super().__init__(CategoryObject)


class AppCategory(Transformer):
    def __init__(self):
        super().__init__(CategoryObject)
