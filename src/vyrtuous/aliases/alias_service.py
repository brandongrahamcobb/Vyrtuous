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

import importlib.util
import inspect
from pathlib import Path

from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.alias import Alias, NotAlias
from vyrtuous.db.database_factory import DatabaseFactory

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


async def delete_alias(alias_name: str, context):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    guild_snowflake = context.guild.id
    where_kwargs = {
        "alias_name": alias_name,
        "guild_snowflake": int(guild_snowflake),
    }
    alias = await database_factory.select(singular=True, **where_kwargs)
    if not alias:
        return f"No aliases found for `{alias_name}`."
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(guild_snowflake)
    channel = guild.get_channel(alias.channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(alias.channel_snowflake)
    if getattr(alias, "role_snowflake"):
        role = guild.get_role(alias.role_snowflake)
        if not role:
            raise commands.RoleNotFound(alias.role_snowflake)
        msg = (
            f"Alias `{alias.alias_name}` of type "
            f"`{alias.category}` for channel {channel.mention if channel else alias.channel_snowflake} "
            f" and role {role.mention} deleted successfully."
        )
    else:
        msg = (
            f"Alias `{alias.alias_name}` of type "
            f"`{alias.category}` for channel {channel.mention} "
            f"deleted successfully."
        )
    await database_factory.delete(**where_kwargs)
    return msg


async def create_alias(alias_name: str, category: str, channel, *, role=None):
    database_factory = DatabaseFactory(MODEL)
    kwargs = {"channel_snowflake": channel.id, "guild_snowflake": channel.guild.id}
    msg = (
        f"Alias `{alias_name}` of type `{category}` "
        f"created successfully for channel {channel.mention}."
    )
    alias = await database_factory.select(
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
    alias = MODEL(alias_name=alias_name, category=str(category), **kwargs)
    await database_factory.create(alias)
    return msg


def alias_category_to_alias(category: str):
    dir_paths = []
    dir_paths.append(Path("/app/vyrtuous"))
    typed_aliases = dir_to_classes(dir_paths=dir_paths)
    for a in typed_aliases:
        if a.category == category:
            return a()
    raise NotAlias()


def dir_to_classes(dir_paths, *, attr="ARGS_MAP"):
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


async def migrate(kwargs):
    database_factory = DatabaseFactory(MODEL)
    await database_factory.update(**kwargs)
