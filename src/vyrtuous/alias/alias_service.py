"""!/bin/python3
alias_service.py The purpose of this program is to extend Service to service aliases.

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

from copy import copy
from dataclasses import dataclass, field
import importlib.util
import inspect
from pathlib import Path
from typing import Dict, List, Union

import discord
from discord.ext import commands

from vyrtuous.alias.alias import Alias


class NotAlias(commands.CheckFailure):
    def __init__(self):
        super().__init__(
            message=("Invalid alias."),
        )


@dataclass(frozen=True)
class AliasDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[str, Union[List, str]]]]]] = (
        field(default_factory=dict)
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


class AliasService:
    __CHUNK_SIZE = 7
    MODEL = Alias
    CATEGORY_TO_HELP = {
        "ban": [
            "**member**: Tag a member or include their ID",
            "**duration**: m/h/d",
            "**reason**: Reason for ban",
        ],
        "flag": [
            "**member**: Tag a member or include their ID",
            "**reason**: Reason for flag",
        ],
        "role": [
            "**member**: Tag a member or include their ID",
            "**role**: Tag a role or include its ID",
        ],
        "tmute": [
            "**member**: Tag a member or include their ID",
            "**duration**: m/h/d",
            "**reason**: Reason for text-mute",
        ],
        "vegan": ["**member**: Tag a member or include their ID"],
        "vmute": [
            "**member**: Tag a member or include their ID",
            "**duration**: m/h/d",
            "**reason**: Reason for voice-mute",
        ],
    }
    CATEGORY_TO_DESCRIPTION = {
        "ban": "Toggles a ban.",
        "flag": "Toggles a moderation flag.",
        "role": "Toggles a role to a user.",
        "tmute": "Toggles a mute in text channels.",
        "vegan": "Toggles a going vegan flag.",
        "vmute": "Toggles a mute in voice channels.",
    }
    CATEGORY_TO_PERMISSION_LEVEL = {
        "ban": "Moderator",
        "flag": "Moderator",
        "role": "Coordinator",
        "tmute": "Moderator",
        "vegan": "Moderator",
        "vmute": "Moderator",
    }

    def __init__(
        self, *, bot=None, database_factory=None, dictionary_service=None, emoji=None
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji

    async def build_dictionary(self, where_kwargs):
        dictionary = {}
        aliases = await self.__database_factory.select(singular=False, **where_kwargs)
        for alias in aliases:
            dictionary.setdefault(alias.guild_snowflake, {"channels": {}})
            dictionary[alias.guild_snowflake]["channels"].setdefault(
                alias.channel_snowflake, {"aliases": {}}
            )
            dictionary[alias.guild_snowflake]["channels"][alias.channel_snowflake][
                "aliases"
            ].setdefault(alias.category, {})[alias.alias_name] = []
            if alias.category == "role":
                guild = self.__bot.get_guild(alias.guild_snowflake)
                if guild:
                    role = guild.get_role(alias.role_snowflake)
                    dictionary[alias.guild_snowflake]["channels"][
                        alias.channel_snowflake
                    ]["aliases"][alias.category][alias.alias_name] = role.mention
        return dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Command Aliases"

        where_kwargs = object_dict.get("columns", None)

        dictionary = await self.build_dictionary(where_kwargs=where_kwargs)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=AliasDictionary, dictionary=dictionary
        )
        if is_at_home:
            pages.extend(processed_dictionary.skipped_channels)
            pages.extend(processed_dictionary.skipped_guilds)

        alias_n = 0
        for guild_snowflake, guild_data in processed_dictionary.data.items():
            pages = []
            field_count = 0
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_dictionary in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                for category, alias_data in channel_dictionary["aliases"].items():
                    lines.append(f"{category}")
                    for name, role_mention in alias_data.items():
                        if category == "role":
                            lines.append(f"  ↳ `{name}` -> {role_mention}")
                        else:
                            lines.append(f"  ↳ `{name}`")
                        alias_n += 1
                    field_count += 1
                if field_count >= self.__CHUNK_SIZE:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = self.__dictionary_service.flush_page(
                        embed, pages, title, guild.name
                    )
                    lines = []
                    field_count = 0
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)
        if pages:
            pages[0].description = f"**({alias_n})**"
        return pages

    async def delete_alias(self, alias_name, context):
        guild_snowflake = context.guild.id
        where_kwargs = {
            "alias_name": alias_name,
            "guild_snowflake": int(guild_snowflake),
        }
        alias = await self.__database_factory.select(singular=True, **where_kwargs)
        if not alias:
            return f"No aliases found for `{alias_name}`."
        guild = self.__bot.get_guild(guild_snowflake)
        channel = guild.get_channel(alias.channel_snowflake)
        if getattr(alias, "role_snowflake"):
            role = guild.get_role(alias.role_snowflake)
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel.mention if channel else alias.channel_snowflake} "
                f" and role {role.mention if role else alias.role_snowflake} deleted successfully."
            )
        else:
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel.mention} "
                f"deleted successfully."
            )
        await self.__database_factory.delete(**where_kwargs)
        return msg

    async def create_alias(self, alias_name, category, channel, *, role=None):
        kwargs = {"channel_snowflake": channel.id, "guild_snowflake": channel.guild.id}
        msg = (
            f"Alias `{alias_name}` of type `{category}` "
            f"created successfully for channel {channel.mention}."
        )
        alias = await self.__database_factory.select(
            category=category,
            alias_name=alias_name,
            guild_snowflake=channel.guild.id,
            singular=True,
        )
        if alias and alias.category != "role":
            return (
                f"Alias of type `{category}` "
                f"with the name `{alias_name} already exists in this server channel ({channel.guild.name})."
            )
        if role:
            kwargs.update({"role_snowflake": role.id})
            msg = (
                f"Alias `{alias_name}` of type `{category}` "
                f"created successfully for channel {channel.mention} with role {role.mention}."
            )
        alias = self.MODEL(alias_name=alias_name, category=str(category), **kwargs)
        await self.__database_factory.create(alias)
        return msg

    def alias_category_to_alias(self, alias_category):
        dir_paths = []
        dir_paths.append(Path("/app/vyrtuous"))
        typed_aliases = self.dir_to_classes(dir_paths=dir_paths)
        for a in typed_aliases:
            if a.category == alias_category:
                return a()
        raise NotAlias()

    def dir_to_classes(self, dir_paths, *, attr="ARGS_MAP"):
        classes = []
        for dir_path in dir_paths:
            for py_file in dir_path.rglob("*.py"):
                if py_file.name == "__init__.py":
                    continue
                module_name = py_file.stem
                spec = importlib.util.spec_from_file_location(module_name, str(py_file))
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                for _, cls in inspect.getmembers(module, inspect.isclass):
                    if cls.__module__ != module.__name__:
                        continue
                    if attr in getattr(cls, "__annotations__", {}):
                        classes.append(cls)
        return classes
